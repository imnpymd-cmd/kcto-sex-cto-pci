#===============================================================================
# 01_clean_raw.R
# K-CTO Registry 2.0 — Sex/Gender Difference Analysis
# Phase 1: raw → cleaned (superset 원칙, 환자 flag 방식)
#
# 원본: data/raw/Data_for_release_v1.xlsx (시트 Data: 9,316행x255열, 시트 codebook)
# 산출:
#   - data/processed/kcto_sex_analysis.rds / .csv   (superset, 9,316행 유지)
#   - data/processed/variable_dictionary.csv        (cleaned 파일 기준)
#   - interim_report/phase1_data/*.csv              (진단용 보조 표, 보고서 근거자료)
#
# 규약: R(tidyverse) 전용. 원본 파일 절대 수정 안 함(읽기 전용). 결측 대치 없음.
#       변수 의미는 반드시 codebook 시트를 근거로 함 (CLAUDE.md §4 변수 해석 규약).
#       set.seed() 미사용 (randomness 도입 안 함 — Phase 1 규칙).
#===============================================================================

suppressMessages({
  library(tidyverse)
  library(readxl)
  library(lubridate)
  library(janitor)
})

proj_root   <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
raw_path    <- file.path(proj_root, "data/raw/Data_for_release_v1.xlsx")
out_dir     <- file.path(proj_root, "data/processed")
report_dir  <- file.path(proj_root, "interim_report/phase1_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# 0. 유틸 함수
#-------------------------------------------------------------------------------

# 문자열 "NA"/빈문자를 실제 NA로 치환 (readxl이 혼합타입 컬럼을 문자형으로 guess하는
# 문제 대응. 과제 지시사항: "문자열 NA가 숫자 컬럼에 섞여 type coercion 경고 발생 →
# 명시적으로 처리")
clean_na_chr <- function(x) {
  x <- as.character(x)
  x[x == "NA" | x == ""] <- NA_character_
  x
}

# "NA" 문자열 처리 후 숫자로 변환. 과학적 표기(E-2 등)는 as.numeric이 정상 처리.
# 숫자로 변환 불가능한 잔존 문자(자유서술 오입력 등)는 NA 처리하고 별도 로그에 기록.
# 주의: log_env는 반드시 environment (list가 아님) — R의 list는 값 전달이라 함수 내부
# 수정이 호출부에 반영되지 않으므로, 참조 의미론을 갖는 environment를 사용함.
clean_num <- function(x, colname = NA_character_, log_env = NULL) {
  x_chr <- clean_na_chr(x)
  x_num <- suppressWarnings(as.numeric(x_chr))
  bad_idx <- which(!is.na(x_chr) & is.na(x_num))
  if (length(bad_idx) > 0 && !is.null(log_env)) {
    log_env$logs[[length(log_env$logs) + 1]] <- tibble(
      variable = colname, row = bad_idx, raw_value = x_chr[bad_idx]
    )
  }
  x_num
}

# 혼합 날짜 파서: raw에는 (1) 텍스트 ISO 형식("YYYY-MM-DD", 일부는 달력상 불가능한
# 날짜 포함: 예 "2022-13-01") 과 (2) 엑셀 시리얼 넘버(문자로 저장됨) 가 혼재.
# 잘못된 달력 날짜는 파싱 불가 → NA 처리하고 로그에 기록 (원본은 수정하지 않음).
parse_mixed_date <- function(x, colname = NA_character_, log_env = NULL) {
  x_chr <- clean_na_chr(x)
  out <- as.Date(rep(NA, length(x_chr)))
  is_iso <- grepl("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", x_chr) & !is.na(x_chr)
  parsed_iso <- suppressWarnings(lubridate::ymd(x_chr[is_iso]))
  out[is_iso] <- parsed_iso
  invalid_idx <- which(is_iso)[is.na(parsed_iso)]
  if (length(invalid_idx) > 0 && !is.null(log_env)) {
    log_env$logs[[length(log_env$logs) + 1]] <- tibble(
      variable = colname, row = invalid_idx, raw_value = x_chr[invalid_idx],
      issue = "달력상 불가능한 날짜(파싱 불가, NA 처리)"
    )
  }
  is_serial <- grepl("^[0-9]+$", x_chr) & !is.na(x_chr) & !is_iso
  out[is_serial] <- as.Date(as.numeric(x_chr[is_serial]), origin = "1899-12-30")
  out
}

#-------------------------------------------------------------------------------
# 1. 원본 로딩 (Data, codebook 시트)
#-------------------------------------------------------------------------------

raw <- suppressWarnings(read_excel(raw_path, sheet = "Data"))
codebook <- read_excel(raw_path, sheet = "codebook", col_types = "text") %>%
  filter(!is.na(`Variable ID`)) %>%
  transmute(var_name = `Variable names`, codebook_desc = `Code descriptions`)

stopifnot(nrow(raw) == 9316, ncol(raw) == 255, nrow(codebook) == 255)
stopifnot(length(setdiff(names(raw), codebook$var_name)) == 0)

cat("원본 로딩 완료: ", nrow(raw), "행 x", ncol(raw), "열\n")

#-------------------------------------------------------------------------------
# 2. 컬럼 분류 (codebook 근거) — 타입 변환 전략 결정
#-------------------------------------------------------------------------------

# (a) 식별자/자유텍스트 (숫자 변환 불가, 그대로 문자형 보존)
id_text_cols <- c("idx", "id_code", "Center_name", "id")

# (b) 날짜 컬럼 (codebook: "Date for...", "The date when...")
date_cols <- c("CTO_intervention_date", "last_FU_date",
               "TVR_date", "TLR_date", "NLR_date", "CABG_date", "ST_date",
               "Sponmi_date", "TV_Sponmi_date", "CVA_date", "bleeding_date",
               "Death_date")

# (c) 라벨형 문자 컬럼 (codebook에 숫자코드 없이 텍스트 그룹명으로 정의됨)
label_text_cols <- c("hos_location", "hos_location_2", "hos_location_3",
                      "hos_location_city", "index_yr_cat", "index_yr_cat5",
                      "approach_direction", "Wire.name", "Stent_name")

# (d) 코드북 설명이 없거나(공란) / 코드북 정의(단일값)와 실제 raw 값(콤마결합·순번
#     표기 등)이 불일치하여 의미가 확정되지 않는 컬럼 — 추측 변환 금지, 원문자 그대로
#     보존. Phase 1 완료 보고에 연구자 질의 항목으로 포함.
ambiguous_multivalue_cols <- c("Approach_2", "Wiring_tech", "Type_monotherapy",
                                "non_CTO.PCI_vessel", "Main_donor_vessel_2",
                                "dc_antiplt", "dc_antiplt_2")

# (e) 나머지 전체 → 숫자형 변환 (codebook상 모두 0/1, 1/2, 순서형 코드 등 숫자 코딩)
all_cols <- names(raw)
numeric_cols <- setdiff(all_cols,
                         c(id_text_cols, date_cols, label_text_cols, ambiguous_multivalue_cols))

cat("컬럼 분류: id/text=", length(id_text_cols),
    " date=", length(date_cols),
    " label_text=", length(label_text_cols),
    " ambiguous_multivalue=", length(ambiguous_multivalue_cols),
    " numeric=", length(numeric_cols), " (합계=", length(all_cols), ")\n", sep = "")

#-------------------------------------------------------------------------------
# 3. 타입 변환 실행 (+ 변환 실패 로그)
#-------------------------------------------------------------------------------

numeric_parse_fail_log <- new.env(); numeric_parse_fail_log$logs <- list()
date_parse_fail_log <- new.env(); date_parse_fail_log$logs <- list()

cleaned <- raw

for (cn in numeric_cols) {
  cleaned[[cn]] <- clean_num(raw[[cn]], colname = cn, log_env = numeric_parse_fail_log)
}
for (cn in date_cols) {
  cleaned[[cn]] <- parse_mixed_date(raw[[cn]], colname = cn, log_env = date_parse_fail_log)
}
for (cn in c(id_text_cols, label_text_cols, ambiguous_multivalue_cols)) {
  cleaned[[cn]] <- clean_na_chr(raw[[cn]])
}

numeric_parse_fail_df <- if (length(numeric_parse_fail_log$logs) > 0) bind_rows(numeric_parse_fail_log$logs) else
  tibble(variable = character(), row = integer(), raw_value = character())
date_parse_fail_df <- if (length(date_parse_fail_log$logs) > 0) bind_rows(date_parse_fail_log$logs) else
  tibble(variable = character(), row = integer(), raw_value = character(), issue = character())

write_csv(numeric_parse_fail_df, file.path(report_dir, "log_numeric_parse_failures.csv"))
write_csv(date_parse_fail_df, file.path(report_dir, "log_date_parse_failures.csv"))

cat("숫자 변환 실패(비숫자 잔존값) 건수:", nrow(numeric_parse_fail_df), "\n")
cat("날짜 파싱 실패(달력상 불가능) 건수:", nrow(date_parse_fail_df), "\n")

#-------------------------------------------------------------------------------
# 4. 파생 변수 — 성별, 환자 키, 코호트 플래그 (행 삭제 없음, superset 원칙)
#-------------------------------------------------------------------------------

# 4-1. 성별: codebook 근거 "Sex; male=1, female=2" → female_bi 파생 (raw에는
#      female_bi 변수가 존재하지 않음 — CLAUDE.md는 legacy 관행을 전제로 존재를
#      가정했으나 실제 raw 255개 컬럼에는 없음. Sex 코드북 정의를 그대로 적용하여
#      Phase 1에서 새로 파생함).
cleaned <- cleaned %>%
  mutate(female_bi = case_when(Sex == 2 ~ 1, Sex == 1 ~ 0, TRUE ~ NA_real_))

# 4-2. 환자 키: CLAUDE.md 지정 "id_code + id" 조합 (row 자체는 idx로 이미 고유).
cleaned <- cleaned %>%
  mutate(patient_key = paste(id_code, id, sep = "__"),
         row_id = row_number())

# 4-3. 등록기간 window flag (parent paper 등록기간 2014-2022)
cleaned <- cleaned %>%
  mutate(flag_in_window = index_yr >= 2014 & index_yr <= 2022)

# 4-4. CTO PCI 시행 flag — codebook 근거 ctopci_bi("CTO PCI binary variable; No=0,
#      Yes=1") 를 그대로 사용 (CTO_manage_3g 파생 불필요, 이미 이진 변수 존재).
cleaned <- cleaned %>%
  mutate(flag_cto_pci = ctopci_bi == 1)

# 4-5. 환자당 index(최초) 시술 flag — patient_key 내 CTO_intervention_date가 가장
#      이른 행을 index로 채택. 동일 최소일자가 여러 행이면 row_id(원본 행 순서)가
#      가장 작은 행을 채택 (tie-break 규칙 명시적 선언).
cleaned <- cleaned %>%
  group_by(patient_key) %>%
  mutate(
    n_in_key = n(),
    flag_index_procedure = row_id == row_id[order(CTO_intervention_date, row_id)][1]
  ) %>%
  ungroup() %>%
  select(-n_in_key)

# 4-6. 1차 분석후보 flag (모든 조건 충족: window AND CTO PCI AND index procedure).
#      ※ 최종 분석 코호트 범위는 CLAUDE.md상 미확정 — 이 flag는 "현재까지 알려진
#      제약 조건을 모두 만족하는 candidate" 이며 Phase 2에서 재확정 필요.
cleaned <- cleaned %>%
  mutate(flag_analysis_primary = flag_in_window & flag_cto_pci & flag_index_procedure)

#-------------------------------------------------------------------------------
# 5. patient_key 중복 진단 (환자 단위 검증 — 동일 환자 재시술 vs id 재사용 판별)
#-------------------------------------------------------------------------------

dup_diag <- cleaned %>%
  select(row_id, patient_key, Sex, Age, CTO_intervention_date) %>%
  group_by(patient_key) %>%
  arrange(CTO_intervention_date, .by_group = TRUE) %>%
  mutate(
    n_in_key = n(),
    n_unique_date = n_distinct(CTO_intervention_date, na.rm = TRUE),
    n_unique_sex = n_distinct(Sex, na.rm = TRUE),
    sex_mismatch = n_unique_sex > 1,
    # 나이 궤적 타당성: 성별 일치 시, 시술 간 경과연수만큼 나이가 증가하는지(오차 ±2세)
    yrs_from_first = as.numeric(difftime(CTO_intervention_date, CTO_intervention_date[1], units = "days")) / 365.25,
    expected_age = Age[1] + yrs_from_first,
    age_diff = abs(Age - expected_age),
    dup_full_row = n_in_key > 1 & duplicated(tibble(CTO_intervention_date, Age, Sex)) |
                   n_in_key > 1 & duplicated(tibble(CTO_intervention_date, Age, Sex), fromLast = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    flag_dupkey_multi = n_in_key > 1,
    flag_dupkey_sex_mismatch = n_in_key > 1 & sex_mismatch,
    flag_dupkey_age_implausible = n_in_key > 1 & !sex_mismatch & !is.na(age_diff) & age_diff > 2,
    flag_dupkey_exact_dup_row = n_in_key > 1 & dup_full_row
  ) %>%
  select(row_id, flag_dupkey_multi, flag_dupkey_sex_mismatch,
         flag_dupkey_age_implausible, flag_dupkey_exact_dup_row)

cleaned <- cleaned %>% left_join(dup_diag, by = "row_id")

#-------------------------------------------------------------------------------
# 6. 로직 정합성 flag (행 삭제 없이 표시만)
#-------------------------------------------------------------------------------

event_fudur_pairs <- list(
  c("Death", "death_fudur"), c("MACE", "MACE_fudur"), c("NACE", "NACE_fudur"),
  c("TVR_re2", "TVR_fudur2"), c("TLR_re2", "TLR_fudur2"), c("NLR_re2", "NLR_fudur2"),
  c("CABG_re", "CABG_fudur"), c("ST", "ST_fudur"), c("TVST", "TVST_fudur"),
  c("Spons_MI_re", "Sponmi_fudur"), c("MI", "MI_fudur"), c("CVA_re", "CVA_fudur"),
  c("bleeding_re", "bleeding_fudur"), c("Bleeding_major", "Bleeding_major_fudur")
)
event_fudur_zero_summary <- map_dfr(event_fudur_pairs, function(p) {
  ev <- cleaned[[p[1]]]; fu <- cleaned[[p[2]]]
  tibble(event_var = p[1], fudur_var = p[2],
         n_event1_fudur_na = sum(ev == 1 & is.na(fu), na.rm = TRUE),
         n_event1_fudur0   = sum(ev == 1 & fu == 0, na.rm = TRUE))
})
write_csv(event_fudur_zero_summary, file.path(report_dir, "log_event_fudur_zero_or_na.csv"))

cleaned <- cleaned %>%
  mutate(
    flag_logic_fudur_exceeds_total =
      (death_fudur > fudur2) | (MACE_fudur > fudur2) | (NACE_fudur > fudur2),
    flag_logic_invalid_last_FU_text = row_id %in%
      which(grepl("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", clean_na_chr(raw$last_FU_date)) &
              is.na(parse_mixed_date(raw$last_FU_date)))
  )

#-------------------------------------------------------------------------------
# 7. 최종 정렬 및 저장 (superset, 9,316행 그대로 유지)
#-------------------------------------------------------------------------------

stopifnot(nrow(cleaned) == 9316)
stopifnot(n_distinct(cleaned$idx) == 9316)

saveRDS(cleaned, file.path(out_dir, "kcto_sex_analysis.rds"))
write_csv(cleaned, file.path(out_dir, "kcto_sex_analysis.csv"))
cat("저장 완료:", file.path(out_dir, "kcto_sex_analysis.rds"), "/",
    file.path(out_dir, "kcto_sex_analysis.csv"), "\n")

#-------------------------------------------------------------------------------
# 8. 변수 카탈로그 (cleaned 파일 기준으로 생성 — legacy와 동일 실수 방지)
#-------------------------------------------------------------------------------

derived_notes <- tribble(
  ~var_name, ~codebook_desc,
  "female_bi", "[Phase1 파생] Sex 코드북 정의(male=1, female=2)를 그대로 적용: Sex==2 -> 1(여성), Sex==1 -> 0(남성). raw 255개 컬럼에는 female_bi가 존재하지 않아 신규 파생함.",
  "patient_key", "[Phase1 파생] id_code + '__' + id 결합 문자열. CLAUDE.md 지정 환자 단위 키. 주의: id 필드 재사용/중복 이슈 있음(02_cleaning_report.md 참조).",
  "row_id", "[Phase1 파생] 원본 행 순서 번호(1~9316). tie-break 및 재현성 목적.",
  "flag_in_window", "[Phase1 파생] index_yr 2014~2022(parent paper 등록기간) 포함 여부.",
  "flag_cto_pci", "[Phase1 파생] ctopci_bi==1 (codebook: CTO PCI binary variable) 과 동일.",
  "flag_index_procedure", "[Phase1 파생] patient_key 내 CTO_intervention_date가 가장 이른 행(동률 시 row_id 최솟값)을 TRUE로 표시.",
  "flag_analysis_primary", "[Phase1 파생] flag_in_window & flag_cto_pci & flag_index_procedure. 최종 분석 코호트 범위는 미확정 — candidate flag.",
  "flag_dupkey_multi", "[Phase1 진단] 동일 patient_key(id_code+id)를 가진 행이 2개 이상 존재.",
  "flag_dupkey_sex_mismatch", "[Phase1 진단] 동일 patient_key 그룹 내 Sex 값이 서로 다름 -> 서로 다른 실제 환자에게 id가 재사용되었을 가능성이 높음(확정적 판단 아님).",
  "flag_dupkey_age_implausible", "[Phase1 진단] 동일 patient_key/동일 Sex 그룹 내에서, 시술일자 간 경과연수 대비 나이 증가가 ±2세 범위를 벗어남 -> id 재사용 의심.",
  "flag_dupkey_exact_dup_row", "[Phase1 진단] 동일 patient_key 그룹 내 (날짜,나이,성별)이 완전히 동일한 행 존재 -> 동일 시술의 중복 입력 의심.",
  "flag_logic_fudur_exceeds_total", "[Phase1 진단] death_fudur/MACE_fudur/NACE_fudur 중 하나라도 총 추적기간(fudur2)을 초과.",
  "flag_logic_invalid_last_FU_text", "[Phase1 진단] last_FU_date가 텍스트 ISO 형식으로 입력되었으나 달력상 존재하지 않는 날짜(예: 11월 31일, 13월)여서 파싱 불가."
)

classify_type <- function(x) {
  if (inherits(x, "Date")) return("date")
  if (is.logical(x)) return("flag(logical)")
  if (is.character(x)) return("character/text")
  if (is.numeric(x)) {
    ux <- unique(x[!is.na(x)])
    if (length(ux) == 0) return("numeric(all-missing)")
    if (all(ux %in% c(0, 1))) return("binary(0/1)")
    if (length(ux) <= 12 && all(ux == round(ux))) return("categorical(numeric code)")
    return("continuous")
  }
  class(x)[1]
}

# PHI(환자 식별정보) 위험 컬럼 — id 필드에는 한글 성명 등 실제 식별정보가 섞여 있음
# (02_cleaning_report.md 참조). 어떤 산출물에도 실제 값 예시를 노출하지 않는다.
phi_risk_cols <- c("id", "patient_key", "id_code", "Center_name")

range_or_levels <- function(x, var_name = NA_character_) {
  xv <- x[!is.na(x)]
  if (length(xv) == 0) return("all missing")
  if (!is.na(var_name) && var_name %in% phi_risk_cols) {
    return(paste0(length(unique(xv)), " distinct values — 예시값 비노출(PHI 위험, id 필드에 ",
                   "환자 성명 등 식별정보 포함 가능; 02_cleaning_report.md 참조)"))
  }
  if (inherits(x, "Date")) return(paste0(min(xv), " ~ ", max(xv)))
  if (is.character(x)) {
    u <- sort(unique(xv))
    if (length(u) <= 15) return(paste(u, collapse = " | "))
    return(paste0(length(u), " distinct values (top5: ",
                   paste(names(sort(table(xv), decreasing = TRUE))[1:5], collapse=", "), ")"))
  }
  if (is.logical(x)) return(paste0("TRUE=", sum(xv), ", FALSE=", sum(!xv)))
  if (is.numeric(x)) {
    ux <- sort(unique(xv))
    if (length(ux) <= 12 && all(ux == round(ux))) return(paste(ux, collapse = ", "))
    return(paste0(round(min(xv), 2), " ~ ", round(max(xv), 2),
                   " (median ", round(median(xv), 2), ")"))
  }
  NA_character_
}

# category / Table1 역할 / model 역할 분류를 위한 그룹 매핑 (변수명 패턴 기반, 문서화)
categorize_var <- function(v) {
  case_when(
    v %in% c("idx","id_code","Center_name","id","patient_key","row_id") ~ "identifier",
    v %in% c("Center_volume","hos_location","hos_location_bi","hos_location_2",
             "hos_location_3","hos_location_city") ~ "site/administrative",
    v %in% c("Age","Sex","female_bi","Height","Weight","BMI_re","bsa","Current_sm",
             "Clindx","clindx_stable","HTN","DM","DM_Insulin","Dyslipidemia","CKD_2",
             "prePCI","preCABG","preMI","preCVA","prePVD","AF") ~ "baseline_demographic_comorbidity",
    v %in% c("Cr","Cr_2d","cr_change","cr_increase","CIN_2","egfr","egfr_2d","esrd_eq",
             "CKMB_pre","CKMB_peak","CKMB_pre_5g","CKMB_pre_cat","CKMB_peak_5g",
             "CKMB_peak_cat","preEF","FU_EF") ~ "baseline_lab_imaging",
    v %in% c("CTO_intervention_date","last_FU_date","fudur2","index_yr","index_yr_cat",
             "index_yr_cat5") ~ "index_date/follow_up_admin",
    v %in% c("CTO_manage_3g","Vessel.Disease","CTO_number","CTO_Vessel","CTO_Vessel_3g",
             "CTO_location","CTO_location_pm","LM.or.LAD.OS","ISR_CTO",
             "number_donor_vessel","Main_donor_vessel_2","lad_donor","lcx_donor",
             "rca_donor","Rentrop","Rentrop3","JCTO_re","JCTO_bi","JCTO_re.g4","Blunt",
             "Calcific","Tortuocity","Length20","Retry") ~ "cto_lesion_characteristics",
    v %in% c("pci_bi","non_CTO.PCI","number_non_CTO.PCI","non_CTO.PCI_vessel",
             "noncto_lad","noncto_lcx","noncto_rca","noncto_lm","noncto_ri","ctopci_bi",
             "CR","cto_success","Bilati_inj","IVUS","approach_direction",
             "ante_approach_only","Approach_2","AWE","ADR","RWE","rCART","CART",
             "Fail_tech","Wiring_tech","Colla_approach_re","Retro_attempt_re",
             "Retro_success","Retro_first","TV_DEB","TV_POBA","TV_STENT","Stent_name",
             "Stent_number_re","n.stent.5g","Final_TIMI","fluorotime","dye") ~ "procedure_technique",
    str_detect(v, "^Stent_(len|dia)_[1-5]$") | v %in% c("stent_dia_mean","stent_len_total",
             "stent_dia_min","stent_dia_max") ~ "procedure_stent_detail",
    v %in% c("n.wire","n.wire.5g","Wire.name","total_Microcath_number","n.micro.3g") ~ "procedure_wire_microcatheter",
    v %in% c("inhos_death","inhos_CVA","urgent_CABG","urgent_PCI",
             "urgent_tampon_intervention","PeriMI") ~ "periprocedural_complication",
    str_detect(v, "^dc_") | v %in% c("dapt","dapt_tapt","cdapt_dc","tdapt_dc","pdapt_dc") ~ "discharge_medication",
    str_detect(v, "_(1m|3m|6m|1Y|2Y|1Y\\.1)$") | str_detect(v, "^(asp|clopid|ticagrelor|prasugrel|cilostazol|cdapt|tdapt|pdapt)_") |
      v %in% c("dapt_time","Type_monotherapy") ~ "longitudinal_medication_adherence",
    str_detect(v, "_date$") ~ "outcome_event_date",
    str_detect(v, "_fudur") ~ "outcome_followup_duration",
    v %in% c("Death","NCVdeath","CVdeath","TVR_re2","TVR_re","TLR_re2","TLR_re",
             "NLR_re2","NLR_re","CABG_re","ST","ST_def","TVST","Spons_MI_re",
             "TV_Sponmi_re","MI","TVMI","CVA_re","bleeding_re","Bleeding_re",
             "Bleeding_major","Bleeding_minor","MACE","NACE") ~ "outcome_event",
    str_detect(v, "^flag_") ~ "phase1_flag",
    TRUE ~ "other/unclassified"
  )
}

role_table1 <- function(v, cat) {
  if (cat %in% c("baseline_demographic_comorbidity","baseline_lab_imaging",
                 "cto_lesion_characteristics") && !str_detect(v, "_5g$|_cat$")) return("candidate")
  if (cat == "outcome_event") return("outcome_summary_candidate")
  if (v %in% c("Age","Sex","female_bi")) return("candidate")
  "no"
}
role_model <- function(v, cat) {
  case_when(
    v == "female_bi" ~ "exposure(sex)",
    cat == "outcome_event" ~ "outcome_candidate",
    cat == "outcome_followup_duration" ~ "time_to_event",
    cat %in% c("baseline_demographic_comorbidity","baseline_lab_imaging",
               "cto_lesion_characteristics") ~ "covariate_candidate",
    cat == "identifier" ~ "id_not_for_model",
    cat == "phase1_flag" ~ "cohort_definition_flag",
    TRUE ~ "undetermined"
  )
}

var_dict <- map_dfr(names(cleaned), function(v) {
  x <- cleaned[[v]]
  n_miss <- sum(is.na(x))
  cat_v <- categorize_var(v)
  tibble(
    var_name = v,
    r_type = class(x)[1],
    type_label = classify_type(x),
    n_missing = n_miss,
    pct_missing = round(n_miss / nrow(cleaned) * 100, 2),
    range_or_levels = range_or_levels(x, var_name = v),
    category = cat_v,
    role_in_table1 = role_table1(v, cat_v),
    role_in_model = role_model(v, cat_v)
  )
})

var_dict <- var_dict %>%
  left_join(bind_rows(codebook, derived_notes), by = "var_name") %>%
  mutate(codebook_desc = str_replace_all(coalesce(codebook_desc, "[정의 없음 — 연구자 확인 필요]"), "[\r\n]+", " ")) %>%
  rename(notes = codebook_desc) %>%
  arrange(match(var_name, names(cleaned)))

write_csv(var_dict, file.path(out_dir, "variable_dictionary.csv"))
cat("변수 카탈로그 저장 완료:", file.path(out_dir, "variable_dictionary.csv"),
    " (", nrow(var_dict), "개 변수)\n")

#-------------------------------------------------------------------------------
# 9. 보조 진단 산출물 (보고서 근거자료) — interim_report/phase1_data/ 에 저장
#-------------------------------------------------------------------------------

# 9-1. STROBE flow 수치
flow <- tibble(
  step = c("raw_total", "flag_in_window_TRUE", "flag_cto_pci_TRUE (전체)",
           "flag_in_window & flag_cto_pci", "flag_index_procedure_TRUE (전체)",
           "flag_analysis_primary_TRUE (candidate)"),
  n = c(nrow(cleaned),
        sum(cleaned$flag_in_window),
        sum(cleaned$flag_cto_pci),
        sum(cleaned$flag_in_window & cleaned$flag_cto_pci),
        sum(cleaned$flag_index_procedure),
        sum(cleaned$flag_analysis_primary))
)
write_csv(flow, file.path(report_dir, "strobe_flow_counts.csv"))

# 9-2. parent paper Table 1 수치 재현
tbl1_repro <- tibble(
  metric = c("N", "Age mean(SD)", "Female %", "HTN %", "DM %", "CKD %",
             "Stable angina %", "Median FU (IQR)"),
  parent_paper = c("9316", "65.9 (11.6)", "19.7", "69.7", "39.9", "24.1", "47.4",
                    "651 [156-1390]"),
  raw_recalculated = c(
    as.character(nrow(cleaned)),
    sprintf("%.1f (%.1f)", mean(cleaned$Age, na.rm=TRUE), sd(cleaned$Age, na.rm=TRUE)),
    sprintf("%.1f", mean(cleaned$female_bi, na.rm=TRUE) * 100),
    sprintf("%.1f", mean(cleaned$HTN, na.rm=TRUE) * 100),
    sprintf("%.1f", mean(cleaned$DM, na.rm=TRUE) * 100),
    sprintf("%.1f", mean(cleaned$CKD_2, na.rm=TRUE) * 100),
    sprintf("%.1f", mean(cleaned$Clindx == 1, na.rm=TRUE) * 100),
    sprintf("%.0f [%.0f-%.0f]", median(cleaned$fudur2, na.rm=TRUE),
            quantile(cleaned$fudur2, .25, na.rm=TRUE), quantile(cleaned$fudur2, .75, na.rm=TRUE))
  ),
  denominator_note = c("N=9316",
    sprintf("결측 제외(n=%d)", sum(!is.na(cleaned$Age))),
    "결측 0 (Sex 완전)",
    sprintf("결측 제외(n=%d, 결측 %d)", sum(!is.na(cleaned$HTN)), sum(is.na(cleaned$HTN))),
    sprintf("결측 제외(n=%d, 결측 %d)", sum(!is.na(cleaned$DM)), sum(is.na(cleaned$DM))),
    sprintf("결측 제외(n=%d, 결측 %d)", sum(!is.na(cleaned$CKD_2)), sum(is.na(cleaned$CKD_2))),
    sprintf("결측 제외(n=%d, 결측 %d)", sum(!is.na(cleaned$Clindx)), sum(is.na(cleaned$Clindx))),
    "결측 0 (fudur2 완전)")
)
write_csv(tbl1_repro, file.path(report_dir, "table1_reproduction_check.csv"))

# 9-3. patient_key 중복 분류 요약
# 주의: flag_dupkey_age_implausible / flag_dupkey_exact_dup_row 는 "행" 단위 플래그라
# 동일 그룹 내에서 행마다 값이 다를 수 있음(예: 3행 중 2행만 완전동일). 그룹 단위 요약을
# 위해서는 group_by 후 우선순위(성별불일치 > 나이궤적비타당 > 완전동일행중복 > 진짜재시술)
# 규칙으로 그룹당 1개 범주만 배정해야 함 (단순 distinct()는 그룹을 중복 집계하는 버그 유발).
# PHI 보호: patient_key(=id_code+id)는 raw id 필드에 환자 성명 등 식별정보를 포함할 수
# 있으므로, interim_report/ (git 보호 미적용 디렉토리) 로 나가는 어떤 산출물에도 절대
# 포함하지 않는다. 그룹 참조는 idx(신규 생성 익명 식별자, "병원코드-일련번호")만 사용.
dup_group_class <- cleaned %>%
  filter(flag_dupkey_multi) %>%
  group_by(patient_key) %>%
  summarise(
    idx_list = paste(idx, collapse = ";"),
    any_sex_mismatch = any(flag_dupkey_sex_mismatch),
    any_age_implausible = any(flag_dupkey_age_implausible),
    any_exact_dup = any(flag_dupkey_exact_dup_row),
    .groups = "drop"
  ) %>%
  mutate(category = case_when(
    any_sex_mismatch ~ "id_reuse_sex_mismatch",
    any_age_implausible ~ "id_reuse_age_implausible",
    any_exact_dup ~ "exact_duplicate_row",
    TRUE ~ "genuine_repeat_procedure"
  )) %>%
  select(-patient_key)  # PHI 위험 컬럼 제거 (idx_list로 대체 참조)

dup_summary <- dup_group_class %>%
  count(category, name = "n_groups") %>%
  mutate(n_dup_groups_total = nrow(dup_group_class))

write_csv(dup_summary, file.path(report_dir, "patient_key_duplicate_summary.csv"))
write_csv(dup_group_class, file.path(report_dir, "patient_key_duplicate_group_detail.csv"))

# 9-4. 결측률 상위 변수 (>20%)
high_missing <- var_dict %>% filter(pct_missing > 20, !str_detect(category, "phase1_flag")) %>%
  arrange(desc(pct_missing)) %>% select(var_name, pct_missing, category, notes)
write_csv(high_missing, file.path(report_dir, "high_missing_variables.csv"))

# 9-5. 성별 간 차등 결측 검정 (금지된 Table1 비교통계 아님 — 결측률 자체의 성별차 검정은
#      과제 지시상 예외로 허용됨)
key_vars_for_missingness <- c("HTN","DM","CKD_2","Dyslipidemia","AF","preEF","egfr",
  "FU_EF","Retro_success","fluorotime","n.wire","JCTO_re","IVUS","Death","MACE","NACE",
  "TVR_re2","death_fudur","MACE_fudur")
key_vars_for_missingness <- intersect(key_vars_for_missingness, names(cleaned))

diff_missing <- map_dfr(key_vars_for_missingness, function(v) {
  x <- cleaned[[v]]
  tab <- table(is.na(x), cleaned$Sex)
  if (nrow(tab) < 2) {
    return(tibble(variable = v, pct_missing_male = 0,
                   pct_missing_female = 0, p_value = NA_real_, test = "no_missing"))
  }
  ct <- suppressWarnings(chisq.test(tab))
  n_male <- sum(cleaned$Sex == 1); n_female <- sum(cleaned$Sex == 2)
  miss_male <- sum(is.na(x) & cleaned$Sex == 1)
  miss_female <- sum(is.na(x) & cleaned$Sex == 2)
  tibble(variable = v,
         pct_missing_male = round(miss_male / n_male * 100, 2),
         pct_missing_female = round(miss_female / n_female * 100, 2),
         p_value = round(ct$p.value, 4),
         test = "chisq.test(is.na(x) ~ Sex)")
})
write_csv(diff_missing, file.path(report_dir, "differential_missingness_by_sex.csv"))

cat("\n=== Phase 1 cleaning script 완료 ===\n")
cat("최종 행수:", nrow(cleaned), " / 변수 수(파생 포함):", ncol(cleaned), "\n")
