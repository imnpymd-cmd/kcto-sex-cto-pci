###############################################################################
# 00_cohort_prep.R  —  공통 코호트 구축 · 재코딩 모듈
#
# 목적: Phase 2 의 모든 스크립트(02 기술통계, 04 주분석)가 **완전히 동일한**
#       코호트/파생변수를 쓰도록 단일 진입점을 제공한다.
#
# 코호트 (CLAUDE.md §1 연구자 확정 2026-08-08 최종):
#   flag_cto_pci == TRUE & index_yr >= 2014 & index_yr <= 2023  ->  N = 6,777
#   ⛔ flag_analysis_primary / flag_index_procedure 사용 금지 (CLAUDE.md §2-1)
#
# Endpoint (CLAUDE.md §1-1, §1-2):
#   scripts/03_derive_endpoints.R 의 derive_endpoints() 로 파생한
#   DMS(primary) / MACCE_pub / MACE_pub / NACE_pub / AnyRevasc 만 사용.
#   ⛔ 데이터셋 기존 MACE / NACE 컬럼 직접 사용 금지.
#
# 본 모듈은 코호트 외에 **공변량 세트 · outcome 목록 · 모델 helper** 도 함께
# 제공한다 (04_primary_analysis.R / 05_qa_sensitivity.R 공용, 중복 정의 금지).
#
# ⛔ PHI: id / patient_key / id_code / Center_name / idx 의 실제 값은 어떤
#         산출물에도 출력하지 않는다. 로드 직후 컬럼을 제거하고 익명 정수
#         center_anon 만 남긴다.
#
# 사용법:
#   source("scripts/00_cohort_prep.R")
#   co <- build_cohort()                       # 주 코호트 (N = 6,777)
#   full <- build_cohort(filter_cohort = FALSE) # 전체 9,316 (출판본 대조용)
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")

source(file.path(PROJ, "scripts", "03_derive_endpoints.R"))

# 코호트 상수 (assertion 용) --------------------------------------------------
COHORT_N        <- 6777L
COHORT_N_MALE   <- 5491L
COHORT_N_FEMALE <- 1286L
COHORT_LABEL    <- "CTO PCI, index year 2014-2023 (N=6,777)"

###############################################################################
# 재코딩 · 파생
###############################################################################

recode_common <- function(d) {

  # --- 노출 -----------------------------------------------------------------
  d$sex_f <- factor(d$Sex, levels = c(1, 2), labels = c("Male", "Female"))
  d$female <- as.integer(d$Sex == 2)
  stopifnot(all(d$female == d$female_bi, na.rm = TRUE))   # female_bi 와 정합성

  # --- CR: codebook 상 Yes = 1 / No = 2 (다른 이진변수와 **역방향**) --------
  #     명시적으로 재코딩하고 아래에서 assertion 으로 방향 역전을 차단한다.
  d$CR_yes <- dplyr::case_when(d$CR == 1 ~ 1, d$CR == 2 ~ 0, TRUE ~ NA_real_)

  # --- 임상 표현형 ----------------------------------------------------------
  d$Clindx_f <- factor(d$Clindx, levels = 1:4,
                       labels = c("Stable angina", "Unstable angina",
                                  "MI (culprit CTO)", "MI (non-culprit CTO)"))
  d$acute_MI <- ifelse(is.na(d$Clindx), NA_real_, as.numeric(d$Clindx %in% c(3, 4)))
  d$acs      <- ifelse(is.na(d$Clindx), NA_real_, as.numeric(d$Clindx %in% c(2, 3, 4)))
  d$stable_angina <- ifelse(is.na(d$Clindx), NA_real_, as.numeric(d$Clindx == 1))

  # --- 시설 ----------------------------------------------------------------
  d$Center_volume_f <- factor(d$Center_volume, levels = 1:4,
                              labels = c("<100", "100-299", "300-599", ">=600"))
  d$hos_location_f <- factor(d$hos_location, levels = c("seoul", "kyeonggi", "others"),
                             labels = c("Seoul", "Kyeonggi", "Others"))

  # --- index year (신 코호트 2014-2023 에 맞춘 2년 구간) --------------------
  d$index_yr <- suppressWarnings(as.numeric(as.character(d$index_yr)))
  d$index_period_f <- cut(d$index_yr,
                          breaks = c(-Inf, 2015, 2017, 2019, 2021, Inf),
                          labels = c("2014-2015", "2016-2017", "2018-2019",
                                     "2020-2021", "2022-2023"))
  d$index_half_f <- factor(ifelse(d$index_yr <= 2018, "2014-2018", "2019-2023"),
                           levels = c("2014-2018", "2019-2023"))

  # --- 병변 ----------------------------------------------------------------
  d$Vessel.Disease_f <- factor(d$Vessel.Disease, levels = 1:3,
                               labels = c("1-vessel", "2-vessel", "3-vessel"))
  d$multivessel <- ifelse(is.na(d$Vessel.Disease), NA_real_,
                          as.numeric(d$Vessel.Disease >= 2))
  d$CTO_number_f <- factor(d$CTO_number, levels = 1:3, labels = c("1", "2", "3"))
  d$CTO_Vessel_f <- factor(d$CTO_Vessel, levels = 1:4,
                           labels = c("LAD", "LCX", "RCA", "LM"))
  d$CTO_location_f <- factor(d$CTO_location, levels = 1:3,
                             labels = c("Proximal", "Mid", "Distal"))
  d$number_donor_vessel_f <- factor(d$number_donor_vessel, levels = 0:3,
                                    labels = c("0", "1", "2", "3"))
  d$Rentrop_f <- factor(d$Rentrop, levels = 0:3,
                        labels = c("Grade 0", "Grade 1", "Grade 2", "Grade 3"))
  d$JCTO_g4_f <- factor(d$JCTO_re.g4, levels = 0:3, labels = c("0", "1", "2", ">=3"))

  # --- 시술 ----------------------------------------------------------------
  d$approach_direction_f <- factor(d$approach_direction,
                                   levels = c("antegrade", "both", "retrograde"),
                                   labels = c("Antegrade only", "Both", "Retrograde"))
  d$Colla_approach_f <- factor(d$Colla_approach_re, levels = 0:2,
                               labels = c("Not done", "Septal channel", "Epicardial channel"))
  d$Retro_attempt_f <- factor(d$Retro_attempt_re, levels = 0:2,
                              labels = c("Not attempted", "Successful collateral crossing",
                                         "Failed collateral crossing"))
  # 역행성 시도 여부 (성공/실패 모두 '시도'로 간주)
  d$retro_attempt_any <- ifelse(is.na(d$Retro_attempt_re), NA_real_,
                                as.numeric(d$Retro_attempt_re %in% c(1, 2)))
  d$Final_TIMI_f <- factor(d$Final_TIMI, levels = 0:3,
                           labels = c("TIMI 0", "TIMI 1", "TIMI 2", "TIMI 3"))
  d$Final_TIMI3 <- ifelse(is.na(d$Final_TIMI), NA_real_, as.numeric(d$Final_TIMI == 3))
  d$n_stent_f <- factor(d$n.stent.5g, levels = 0:4, labels = c("0", "1", "2", "3", ">=4"))
  d$n_wire_f  <- factor(d$n.wire.5g,  levels = 1:5, labels = c("1", "2", "3", "4", ">=5"))
  d$n_micro_f <- factor(d$n.micro.3g, levels = 0:3, labels = c("0", "1", "2", ">=3"))

  # --- 퇴원 항혈소판 (CLAUDE.md §2-2 확정 정의) ----------------------------
  d$dc_antiplt2_f <- factor(d$dc_antiplt_2,
                            levels = c("ac", "ap", "at", "single", "triple"),
                            labels = c("Aspirin + clopidogrel",
                                       "Aspirin + prasugrel",
                                       "Aspirin + ticagrelor",
                                       "P2Y12 inhibitor monotherapy",
                                       "Triple antiplatelet therapy"))
  d$dapt_potent <- dplyr::case_when(          # 분모: 2제요법(DAPT) 환자만
    d$dc_antiplt_2 %in% c("ap", "at") ~ 1,
    d$dc_antiplt_2 %in% c("ac")       ~ 0,
    TRUE ~ NA_real_)

  # --- 기타 파생 -----------------------------------------------------------
  d$age75  <- ifelse(is.na(d$Age),   NA_real_, as.numeric(d$Age >= 75))
  d$lvef40 <- ifelse(is.na(d$preEF), NA_real_, as.numeric(d$preEF < 40))
  d$fu_zero <- as.numeric(d$fudur2 == 0)

  # 원내 합병증 복합 (6개 구성요소 모두 결측 0)
  d$inhos_any <- as.numeric(
    d$inhos_death == 1 | d$inhos_CVA == 1 | d$urgent_CABG == 1 |
      d$urgent_PCI == 1 | d$urgent_tampon_intervention == 1 | d$PeriMI == 1)

  d
}

###############################################################################
# 코호트 빌더
###############################################################################

build_cohort <- function(filter_cohort = TRUE, verbose = TRUE) {

  raw <- readRDS(file.path(PROJ, "data", "processed", "kcto_sex_analysis.rds"))
  stopifnot(nrow(raw) == 9316)

  # 익명 센터 코드 (실제 병원명/코드 노출 금지 — clustering 전용)
  raw$center_anon <- as.integer(factor(raw$id_code))
  PHI_COLS <- c("id", "patient_key", "id_code", "Center_name", "idx")
  raw <- raw[, setdiff(names(raw), PHI_COLS), drop = FALSE]

  # 출판본 endpoint 파생 (CLAUDE.md §1-1)
  raw <- derive_endpoints(raw)

  raw$index_yr <- suppressWarnings(as.numeric(as.character(raw$index_yr)))

  if (!filter_cohort) return(recode_common(raw))

  d <- raw[raw$flag_cto_pci %in% TRUE &
             !is.na(raw$index_yr) &
             raw$index_yr >= 2014 & raw$index_yr <= 2023, , drop = FALSE]

  d <- recode_common(d)

  # --- 코호트 assertion ----------------------------------------------------
  stopifnot(nrow(d) == COHORT_N)
  stopifnot(sum(d$sex_f == "Male")   == COHORT_N_MALE)
  stopifnot(sum(d$sex_f == "Female") == COHORT_N_FEMALE)

  # --- CR 역방향 재코딩 assertion (Yes=1 / No=2) ---------------------------
  stopifnot(sum(d$CR_yes == 1, na.rm = TRUE) == sum(d$CR == 1, na.rm = TRUE))
  stopifnot(sum(d$CR_yes == 0, na.rm = TRUE) == sum(d$CR == 2, na.rm = TRUE))
  stopifnot(all(is.na(d$CR) == is.na(d$CR_yes)))

  # --- CVA_re 재코딩 assertion (값 2 -> 0) ---------------------------------
  stopifnot(all(d$CVA_re %in% c(0, 1) | is.na(d$CVA_re)))

  # --- 파생 endpoint 정합성 assertion --------------------------------------
  stopifnot(all(d$DMS[d$Death == 1] == 1))
  stopifnot(all(d$MACCE_pub[d$DMS == 1] == 1))          # DMS ⊂ MACCE_pub
  stopifnot(all(d$MACCE_pub[d$AnyRevasc == 1] == 1))    # AnyRevasc ⊂ MACCE_pub
  stopifnot(all(d$NACE_pub[d$MACCE_pub == 1] == 1))     # MACCE_pub ⊂ NACE_pub

  if (verbose) {
    message(sprintf("[cohort] %s | Male %d / Female %d (%.1f%% female) | centers %d",
                    COHORT_LABEL, sum(d$sex_f == "Male"), sum(d$sex_f == "Female"),
                    100 * mean(d$sex_f == "Female"), length(unique(d$center_anon))))
  }
  d
}

###############################################################################
# day-0 event 처리 (주 규칙: 0.5일 대체 / 민감도: 제외)
#   근거: event==1 & fudur==0 은 시술 당일·재원 중 발생한 실제 사건이다.
#         제외하면 가장 중증인 사건만 선택적으로 소실되어 누적발생률이
#         하향 편향되고, 성별 원내 사건률이 다르면 노출과 상관된 차등 편향이 된다.
###############################################################################

apply_day0 <- function(time, event, rule = c("half_day", "exclude")) {
  rule <- match.arg(rule)
  keep <- !is.na(time) & !is.na(event)
  if (rule == "half_day") {
    time[keep & event == 1 & time == 0] <- 0.5
  } else {
    keep <- keep & !(event == 1 & time == 0)
  }
  list(time = time, event = event, keep = keep)
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

###############################################################################
# 공변량 세트 (DAG 기반, parent 출판본 공변량을 출발점으로 함)
#
#  parent 출판본 Methods 에 기술된 공변량 세트:
#   age, sex, BMI, current smoking, hypertension, DM, CKD, prior PCI, prior CABG,
#   prior PVD, prior CVA, prior MI, AF, stable angina, proximal-to-mid CTO,
#   ISR CTO, Rentrop grade 3, J-CTO score, concurrent non-CTO PCI
#
#  본 연구의 조정 (사유는 보고서 §3-1 에 기술):
#   (+) index period      — 등록기간 2014-2023 의 시대 효과(기법·약제 변화) 보정
#   (+) center volume     — 시설 경험도. 동시에 center_anon 로 clustered robust SE
#   (+) 다혈관 질환       — CTO 외 질환 부담. parent 는 CTO PCI '시도' 예측 모형이라
#                            미포함이었으나 본 연구는 '예후' 모형이므로 필수
#   (-) height / weight / BSA — 성별과 사실상 결정적으로 연관(SMD -1.2 ~ -2.4).
#       positivity 를 위반하므로 교란변수가 아니라 노출의 구성요소로 간주하고
#       주 모델에서 제외. 체격은 BMI(SMD -0.19)로만 대표. (보고서 §3-1, §5-2)
#       ⚠️ 단 시술전략 모형에서는 체격이 임상적으로 정당한 결정요인일 수 있으므로
#          05_qa_sensitivity.R(A3) 에서 공통지지 제한 재분석으로 별도 검증한다.
#   (-) 시술 변수(성공, 완전재관류, 항혈소판 요법) — 노출(성별) 이후에 결정되는
#       **잠재적 매개변수**. 총효과 추정 시 보정하면 편향되므로 주 모델에서 제외하고
#       Model 5(탐색적)에서만 추가.
#   (~) LVEF — 예후에 중요하나 결측 20.3%. 주 모델 제외, Model 4L 및 다중대치에서 평가.
###############################################################################

COV_AGE   <- c("Age")
COV_CLIN  <- c("BMI_re", "Current_sm", "HTN", "DM", "CKD_2",
               "prePCI", "preCABG", "preMI", "preCVA", "prePVD", "AF",
               "stable_angina")
COV_LES   <- c("JCTO_re", "CTO_location_pm", "ISR_CTO", "Rentrop3",
               "non_CTO.PCI", "Vessel.Disease_f")
COV_ADMIN <- c("index_period_f", "Center_volume_f")

COV_M2 <- COV_AGE
COV_M3 <- c(COV_AGE, COV_CLIN)
COV_M4 <- c(COV_AGE, COV_CLIN, COV_LES, COV_ADMIN)              # PRIMARY
COV_M4L<- c(COV_M4, "preEF")                                     # LVEF 추가(민감도)
COV_M5 <- c(COV_M4, "cto_success", "CR_yes", "dc_antiplt2_f", "dc_statin")  # 탐색적

# 희소 사건(events per parameter < 10) 전용 절약형 보정 세트.
COV_PARS <- c("Age", "DM", "CKD_2", "AF")

# IPTW / PSM 용 balancing 변수 (Model 4 와 동일; 체격은 BMI 로만 대표)
COV_PS <- COV_M4

# 센터를 평균 구조에 직접 넣는 모형(고정효과 / 랜덤절편)에서는 Center_volume_f 가
# 센터에 완전히 흡수되므로 제외한다 (05_qa_sensitivity.R A1).
COV_M4_NOCENTER <- setdiff(COV_M4, "Center_volume_f")

TAU5 <- 1825   # 5년 (일)

###############################################################################
# outcome 목록 (endpoint 는 전부 derive_endpoints() 파생본)
###############################################################################

OUT <- list(
  DMS       = list(ev = "DMS",            tm = "DMS_fudur",            lab = "DMS (death, MI, or stroke) - PRIMARY"),
  MACCE_pub = list(ev = "MACCE_pub",      tm = "MACCE_pub_fudur",      lab = "MACCE (published definition)"),
  MACE_pub  = list(ev = "MACE_pub",       tm = "MACE_pub_fudur",       lab = "MACE (published definition)"),
  NACE_pub  = list(ev = "NACE_pub",       tm = "NACE_pub_fudur",       lab = "NACE (published definition)"),
  Death     = list(ev = "Death",          tm = "death_fudur",          lab = "All-cause death"),
  CVdeath   = list(ev = "CVdeath",        tm = "death_fudur",          lab = "Cardiovascular death"),
  NCVdeath  = list(ev = "NCVdeath",       tm = "death_fudur",          lab = "Non-cardiovascular death"),
  AnyRevasc = list(ev = "AnyRevasc",      tm = "AnyRevasc_fudur",      lab = "Any coronary revascularization"),
  MI        = list(ev = "MI",             tm = "MI_fudur",             lab = "Myocardial infarction"),
  Stroke    = list(ev = "CVA_re",         tm = "CVA_fudur",            lab = "Stroke"),
  MajBld    = list(ev = "Bleeding_major", tm = "Bleeding_major_fudur", lab = "TIMI major bleeding")
)

# 관리 전략(management strategy) 지표 — 단일 검정 family (05 A2 의 다중비교 단위)
STRAT_VARS <- list(
  c("IVUS",               "Intravascular ultrasound use"),
  c("Bilati_inj",         "Bilateral injection"),
  c("retro_attempt_any",  "Retrograde approach attempted"),
  c("ante_approach_only", "Antegrade approach only"),
  c("dapt_potent",        "Potent P2Y12 inhibitor at discharge (among DAPT patients)"),
  c("dc_statin",          "Statin at discharge"))

###############################################################################
# 공통 모델 helper
###############################################################################

# 생존 데이터 구성 (day-0 처리 포함)
make_surv <- function(dat, key, day0 = "half_day") {
  o <- OUT[[key]]
  a <- apply_day0(dat[[o$tm]], dat[[o$ev]], rule = day0)
  keep <- a$keep & !is.na(dat$sex_f)
  d <- dat[keep, , drop = FALSE]
  d$.time <- a$time[keep]; d$.event <- a$event[keep]
  d
}

# Cox 적합 + female 항 추출 (center_anon clustered robust SE)
fit_cox <- function(d, covs, label, outcome_lab, extra_note = "") {
  rhs <- if (length(covs)) paste(c("female", covs), collapse = " + ") else "female"
  f <- stats::as.formula(paste0("Surv(.time, .event) ~ ", rhs))
  m <- survival::coxph(f, data = d, cluster = d$center_anon,
                       robust = TRUE, x = TRUE, ties = "efron")
  s <- summary(m)
  i <- which(rownames(s$coefficients) == "female")
  data.frame(
    outcome = outcome_lab, model = label,
    n = m$n, n_events = m$nevent,
    n_female = sum(d$female == 1), events_female = sum(d$.event[d$female == 1] == 1),
    n_male = sum(d$female == 0), events_male = sum(d$.event[d$female == 0] == 1),
    HR = unname(s$conf.int[i, 1]),
    lo = unname(s$conf.int[i, 3]), hi = unname(s$conf.int[i, 4]),
    HR_fmt = sprintf("%.2f (%.2f-%.2f)", s$conf.int[i, 1], s$conf.int[i, 3], s$conf.int[i, 4]),
    se_robust = unname(s$coefficients[i, "robust se"]),
    p_value = unname(s$coefficients[i, "Pr(>|z|)"]),
    p_fmt = fmt_p(unname(s$coefficients[i, "Pr(>|z|)"])),
    note = extra_note, stringsAsFactors = FALSE)
}

# 전체 계수표 추출
tidy_cox <- function(m, label, outcome_lab) {
  s <- summary(m)
  data.frame(outcome = outcome_lab, model = label,
             term = rownames(s$coefficients),
             HR = s$conf.int[, 1], lo = s$conf.int[, 3], hi = s$conf.int[, 4],
             HR_fmt = sprintf("%.2f (%.2f-%.2f)", s$conf.int[, 1], s$conf.int[, 3], s$conf.int[, 4]),
             se_robust = s$coefficients[, "robust se"],
             p_value = s$coefficients[, "Pr(>|z|)"],
             p_fmt = fmt_p(s$coefficients[, "Pr(>|z|)"]),
             stringsAsFactors = FALSE, row.names = NULL)
}

cox_model <- function(d, covs) {
  rhs <- if (length(covs)) paste(c("female", covs), collapse = " + ") else "female"
  survival::coxph(stats::as.formula(paste0("Surv(.time, .event) ~ ", rhs)),
                  data = d, cluster = d$center_anon, robust = TRUE, x = TRUE, ties = "efron")
}

# 보정 로지스틱 (Model 4 공변량 기본) + center_anon cluster-robust SE
logit_sex <- function(dat, yvar, lab, covs = COV_M4) {
  d <- dat[!is.na(dat[[yvar]]), , drop = FALSE]
  d$.y <- d[[yvar]]
  f <- stats::as.formula(paste0(".y ~ female + ", paste(covs, collapse = " + ")))
  m <- stats::glm(f, data = d, family = stats::binomial())
  vc <- sandwich::vcovCL(m, cluster = d$center_anon)
  ct <- lmtest::coeftest(m, vcov. = vc)
  b <- ct["female", "Estimate"]; se <- ct["female", "Std. Error"]
  data.frame(outcome = lab, model = "Adjusted logistic (Model 4 covariates, cluster-robust SE)",
             n = nrow(d), n_events = sum(d$.y == 1),
             events_male = sum(d$.y == 1 & d$female == 0),
             events_female = sum(d$.y == 1 & d$female == 1),
             OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
             OR_fmt = sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)),
             p_value = ct["female", "Pr(>|z|)"], p_fmt = fmt_p(ct["female", "Pr(>|z|)"]),
             stringsAsFactors = FALSE)
}

# 미보정 로지스틱 (cluster-robust SE)
crude_or <- function(dat, yvar, lab) {
  d <- dat[!is.na(dat[[yvar]]), , drop = FALSE]; d$.y <- d[[yvar]]
  m <- stats::glm(.y ~ female, data = d, family = stats::binomial())
  vc <- sandwich::vcovCL(m, cluster = d$center_anon)
  ct <- lmtest::coeftest(m, vcov. = vc)
  b <- ct["female", "Estimate"]; se <- ct["female", "Std. Error"]
  data.frame(outcome = lab, model = "Unadjusted logistic (cluster-robust SE)",
             n = nrow(d), n_events = sum(d$.y == 1),
             events_male = sum(d$.y == 1 & d$female == 0),
             events_female = sum(d$.y == 1 & d$female == 1),
             OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
             OR_fmt = sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)),
             p_value = ct["female", "Pr(>|z|)"], p_fmt = fmt_p(ct["female", "Pr(>|z|)"]),
             stringsAsFactors = FALSE)
}
