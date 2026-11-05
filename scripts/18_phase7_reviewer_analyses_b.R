# =============================================================================
# 18_phase7_reviewer_analyses_b.R
#   Phase 7 (peer review) 추가 분석 2차 — JACC: Advances JADV102226-3033-T
#   연구자 승인 2026-11-05 (registry §4 권고안 승인) 후 수행
#
# 실행:  Rscript scripts/18_phase7_reviewer_analyses_b.R
# 선행:  scripts/17_phase7_reviewer_analyses.R (A-3·A-6·A-9·A-8·A-7·A-5·A-2①②)
#
# 본 스크립트의 세 분석은 모두 **심사 의견에 답하기 위해 추가**되었다.
# ⛔ "prespecified / pre-registered" 긍정 주장 금지 (CLAUDE.md §1-2).
#
#   A-1  (R2-M1 / R1-W5) 전체 레지스트리 2014-2023: 성별 x CTO PCI 시행 여부
#                        baseline 표 + 시행 odds(여 vs 남) 단계적 보정
#   A-4  (R2-M4)         bilateral injection 의 "적응 가능 상황" 제한 분석 +
#                        조영제량·투시시간·신장 안전성 성별 비교
#   A-2③ (R2-M2)         person-time 미기여 환자에 대한 delta-adjusted
#                        (tipping-point) 민감도, 0-1년 창
#
# 코호트: 시행군 = flag_cto_pci & index_yr 2014-2023 (N = 6,777, 본 연구 코호트)
#         미시행군 = 같은 기간의 CTO PCI 미시행 환자 (n = 2,008) — **단일군**
#         (구성은 각주로 공개: 30일 내 CABG / 비CTO PCI / 약물치료)
#
# ⛔ PHI: 센터는 익명 정수 center_anon 으로만 보고한다.
# ⚠️ 관찰연구다. 모든 추정치는 연관성(association)이며 인과 해석·용어 금지.
#
# 출력: output/tables/table35_A1_registry_pci_vs_nopci_by_sex.csv
#       output/tables/table36_A4_bilateral_indication_and_resource.csv
#       output/tables/table37_A2c_delta_adjusted_tipping_point.csv
#       interim_report/phase7_rebuttal/reports/phase7b_analyses_log.txt
#       interim_report/phase7_rebuttal/reports/phase7b_sessionInfo.txt
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(survival); library(lme4); library(sandwich); library(lmtest)
})
set.seed(42)
SEED <- 42L
M_IMP <- 20L        # delta-adjusted 민감도의 대치 횟수

PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "00_cohort_prep.R"))
DIR_TAB <- file.path(PROJ, "output", "tables")
DIR_REP <- file.path(PROJ, "interim_report", "phase7_rebuttal", "reports")
dir.create(DIR_REP, recursive = TRUE, showWarnings = FALSE)

LOG <- character(0)
logmsg <- function(...) { s <- sprintf(...); LOG <<- c(LOG, s); message(s) }

CAV <- paste0(
  "Added in response to peer review (JACC: Advances, reviewer comments received 2026-11-04); ",
  "not part of the analysis set of the submitted manuscript. Observational secondary analysis of a ",
  "registry: all estimates are associations.")

rw <- function(analysis, outcome, model, n = NA_integer_, n_female = NA_integer_,
               estimate = NA_real_, lo = NA_real_, hi = NA_real_,
               est_fmt = NA_character_, p_fmt = NA_character_,
               note = "", caveat = CAV, ...) {
  base <- data.frame(analysis = analysis, outcome = outcome, model = model, n = n,
                     n_female = n_female, estimate = estimate, lo = lo, hi = hi,
                     est_fmt = est_fmt, p_fmt = p_fmt, note = note, caveat = caveat,
                     stringsAsFactors = FALSE)
  ex <- list(...); for (nm in names(ex)) base[[nm]] <- ex[[nm]]
  base
}
fmt2 <- function(e, l, h) sprintf("%.2f (%.2f-%.2f)", e, l, h)

# output 인덱스 확인 (관찰 0122: 누적 디렉토리에서 파일 인덱스는 상태다).
# table35-37 은 본 스크립트의 산출물이며, 최초 실행 전 비어 있음을 확인했다.
# 재실행 시에는 자신의 산출물만 덮어쓴다.
OUT_FILES <- c("table35_A1_registry_pci_vs_nopci_by_sex.csv",
               "table36_A4_bilateral_indication_and_resource.csv",
               "table37_A2c_delta_adjusted_tipping_point.csv")
for (f in OUT_FILES) if (file.exists(file.path(DIR_TAB, f)))
  message("[note] overwriting this script's own earlier output: ", f)

###############################################################################
# 코호트
###############################################################################
co   <- build_cohort()                                  # 시행군 6,777
full <- build_cohort(filter_cohort = FALSE, verbose = FALSE)
reg  <- full[!is.na(full$index_yr) & full$index_yr >= 2014 & full$index_yr <= 2023, , drop = FALSE]
stopifnot(nrow(reg) == 8785, sum(reg$flag_cto_pci) == 6777)
reg$pci <- as.integer(reg$flag_cto_pci)
an <- co[stats::complete.cases(co[, COV_M4, drop = FALSE]), , drop = FALSE]
logmsg("[cohort] registry 2014-2023 N = %d | CTO PCI performed %d (women %d) | not performed %d (women %d) | centres %d (PCI %d / no-PCI %d)",
       nrow(reg), sum(reg$pci == 1), sum(reg$pci == 1 & reg$female == 1),
       sum(reg$pci == 0), sum(reg$pci == 0 & reg$female == 1),
       length(unique(reg$center_anon)), length(unique(reg$center_anon[reg$pci == 1])),
       length(unique(reg$center_anon[reg$pci == 0])))

###############################################################################
# A-1 (R2-M1 / R1-W5) — 전체 레지스트리: 성별 x CTO PCI 시행 여부
###############################################################################
message("\n================ A-1: registry, PCI performed vs not performed ================")

# --- 미시행군 구성 (단일군으로 분석하고 구성만 공개) -----------------------
np <- reg[reg$pci == 0, , drop = FALSE]
cabg30 <- np$CABG_re == 1 & np$CABG_fudur <= 30
nctopci <- np$non_CTO.PCI == 1
cabg30[is.na(cabg30)] <- FALSE; nctopci[is.na(nctopci)] <- FALSE
COMP_NOTE <- sprintf(paste0(
  "The group in whom CTO percutaneous coronary intervention was not performed is analysed as a ",
  "single group and is not subdivided. Its composition within 2014-2023: coronary artery bypass ",
  "grafting within 30 days n = %d, percutaneous coronary intervention of a non-CTO lesion only ",
  "n = %d (overlap with early surgery n = %d), medical therapy alone n = %d. The corresponding ",
  "whole-registry (2003-2023) counts are 74 and 675."),
  sum(cabg30), sum(nctopci), sum(cabg30 & nctopci), sum(!cabg30 & !nctopci))
logmsg("[A-1] no-PCI composition (2014-2023): early CABG %d | non-CTO PCI %d | overlap %d | medical %d",
       sum(cabg30), sum(nctopci), sum(cabg30 & nctopci), sum(!cabg30 & !nctopci))

a1 <- list(); addA1 <- function(x) a1[[length(a1) + 1]] <<- x
addA1(rw("A-1(composition)", "Group definition",
         "CTO PCI not performed, 2014-2023 (single group)", nrow(np), sum(np$female == 1),
         NA_real_, NA_real_, NA_real_,
         sprintf("n = %d (women %d, %.1f%%)", nrow(np), sum(np$female == 1),
                 100 * mean(np$female == 1)), NA_character_, COMP_NOTE,
         measure = "group size"))

# --- baseline 표: 변수 x 성별, 시행 vs 미시행 (SMD + P) --------------------
smd_bin <- function(p1, p0) (p1 - p0) / sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
smd_con <- function(x1, x0) (mean(x1) - mean(x0)) / sqrt((stats::var(x1) + stats::var(x0)) / 2)

row_bin <- function(d, v, lab) {
  out <- list()
  for (g in c("Male", "Female")) {
    s <- d[d$sex_f == g & !is.na(d[[v]]), , drop = FALSE]
    x1 <- s[[v]][s$pci == 1]; x0 <- s[[v]][s$pci == 0]
    if (!length(x1) || !length(x0)) next
    tb <- table(factor(s$pci, 0:1), factor(s[[v]], c(0, 1)))
    p <- tryCatch(if (any(suppressWarnings(stats::chisq.test(tb)$expected) < 5))
      stats::fisher.test(tb)$p.value else stats::chisq.test(tb)$p.value,
      error = function(e) NA_real_)
    out[[length(out) + 1]] <- rw("A-1(baseline)", lab,
      sprintf("%s: CTO PCI performed vs not performed", ifelse(g == "Male", "Men", "Women")),
      nrow(s), sum(s$sex_f == "Female"), smd_bin(mean(x1), mean(x0)), NA_real_, NA_real_,
      sprintf("%d/%d (%.1f%%) vs %d/%d (%.1f%%)", sum(x1), length(x1), 100 * mean(x1),
              sum(x0), length(x0), 100 * mean(x0)), fmt_p(p),
      paste0("Counts and percentages with the standardized mean difference between the group in ",
             "whom CTO PCI was performed and the group in whom it was not, computed within sex. ",
             "P value: chi-squared test (Fisher's exact test if any expected count < 5). ", COMP_NOTE),
      measure = "proportion (%), PCI performed vs not performed", sex = g,
      pct_pci = round(100 * mean(x1), 1), pct_nopci = round(100 * mean(x0), 1),
      n_pci = length(x1), n_nopci = length(x0),
      n_missing_pci = sum(d$sex_f == g & d$pci == 1 & is.na(d[[v]])),
      n_missing_nopci = sum(d$sex_f == g & d$pci == 0 & is.na(d[[v]])))
  }
  dplyr::bind_rows(out)
}
row_con <- function(d, v, lab) {
  out <- list()
  for (g in c("Male", "Female")) {
    s <- d[d$sex_f == g & !is.na(d[[v]]), , drop = FALSE]
    x1 <- s[[v]][s$pci == 1]; x0 <- s[[v]][s$pci == 0]
    if (length(x1) < 3 || length(x0) < 3) next
    p <- tryCatch(stats::t.test(x1, x0)$p.value, error = function(e) NA_real_)
    out[[length(out) + 1]] <- rw("A-1(baseline)", lab,
      sprintf("%s: CTO PCI performed vs not performed", ifelse(g == "Male", "Men", "Women")),
      nrow(s), sum(s$sex_f == "Female"), smd_con(x1, x0), NA_real_, NA_real_,
      sprintf("%.1f (%.1f) vs %.1f (%.1f)", mean(x1), stats::sd(x1), mean(x0), stats::sd(x0)),
      fmt_p(p),
      paste0("Mean (SD) with the standardized mean difference, computed within sex. P value: ",
             "Welch t-test. ", COMP_NOTE),
      measure = "mean (SD), PCI performed vs not performed", sex = g,
      mean_pci = round(mean(x1), 2), mean_nopci = round(mean(x0), 2),
      n_pci = length(x1), n_nopci = length(x0),
      n_missing_pci = sum(d$sex_f == g & d$pci == 1 & is.na(d[[v]])),
      n_missing_nopci = sum(d$sex_f == g & d$pci == 0 & is.na(d[[v]])))
  }
  dplyr::bind_rows(out)
}
row_cat <- function(d, v, lab) {
  out <- list()
  for (g in c("Male", "Female")) {
    s <- d[d$sex_f == g & !is.na(d[[v]]), , drop = FALSE]
    tb <- table(factor(s$pci, 0:1), droplevels(factor(s[[v]])))
    p <- tryCatch(if (any(suppressWarnings(stats::chisq.test(tb)$expected) < 5))
      stats::fisher.test(tb, simulate.p.value = TRUE, B = 2e4)$p.value
      else stats::chisq.test(tb)$p.value, error = function(e) NA_real_)
    for (lv in colnames(tb)) {
      x1 <- s[[v]][s$pci == 1] == lv; x0 <- s[[v]][s$pci == 0] == lv
      out[[length(out) + 1]] <- rw("A-1(baseline)", sprintf("%s: %s", lab, lv),
        sprintf("%s: CTO PCI performed vs not performed", ifelse(g == "Male", "Men", "Women")),
        nrow(s), sum(s$sex_f == "Female"), smd_bin(mean(x1), mean(x0)), NA_real_, NA_real_,
        sprintf("%d/%d (%.1f%%) vs %d/%d (%.1f%%)", sum(x1), length(x1), 100 * mean(x1),
                sum(x0), length(x0), 100 * mean(x0)), fmt_p(p),
        paste0("Category of a multi-level variable; the P value is the overall test across all ",
               "categories within sex and is repeated on each category row. ", COMP_NOTE),
        measure = "proportion (%), PCI performed vs not performed", sex = g,
        pct_pci = round(100 * mean(x1), 1), pct_nopci = round(100 * mean(x0), 1),
        n_pci = length(x1), n_nopci = length(x0),
        n_missing_pci = sum(d$sex_f == g & d$pci == 1 & is.na(d[[v]])),
        n_missing_nopci = sum(d$sex_f == g & d$pci == 0 & is.na(d[[v]])))
    }
  }
  dplyr::bind_rows(out)
}

reg$jcto2  <- ifelse(is.na(reg$JCTO_re), NA_real_, as.numeric(reg$JCTO_re >= 2))
reg$lvef40 <- ifelse(is.na(reg$preEF), NA_real_, as.numeric(reg$preEF < 40))
addA1(row_con(reg, "Age", "Age, years"))
addA1(row_bin(reg, "DM", "Diabetes mellitus"))
addA1(row_bin(reg, "HTN", "Hypertension"))
addA1(row_bin(reg, "CKD_2", "Chronic kidney disease (eGFR <60 mL/min/1.73 m2)"))
addA1(row_bin(reg, "preCABG", "Previous CABG"))
addA1(row_bin(reg, "preMI", "Previous myocardial infarction"))
addA1(row_bin(reg, "lvef40", "Left ventricular ejection fraction <40%"))
addA1(row_con(reg, "preEF", "Left ventricular ejection fraction, %"))
addA1(row_cat(reg, "Clindx_f", "Clinical presentation"))
addA1(row_con(reg, "JCTO_re", "J-CTO score"))
addA1(row_bin(reg, "jcto2", "J-CTO score >=2"))
addA1(row_bin(reg, "Rentrop3", "Rentrop grade 3 collateral"))
addA1(row_cat(reg, "CTO_Vessel_f", "CTO vessel"))
addA1(row_bin(reg, "multivessel", "Multivessel disease"))
addA1(row_cat(reg, "index_period_f", "Index treatment-era period"))

# 결측 보고 (변수별·군별)
for (v in c("preEF", "JCTO_re", "Rentrop", "CKD_2", "Clindx", "multivessel")) {
  for (g in c("Male", "Female")) for (pc in 1:0) {
    s <- reg[reg$sex_f == g & reg$pci == pc, , drop = FALSE]
    addA1(rw("A-1(missingness)", v,
      sprintf("%s, CTO PCI %s", ifelse(g == "Male", "men", "women"),
              ifelse(pc == 1, "performed", "not performed")),
      nrow(s), ifelse(g == "Female", nrow(s), 0L), 100 * mean(is.na(s[[v]])), NA_real_, NA_real_,
      sprintf("%d/%d (%.1f%%)", sum(is.na(s[[v]])), nrow(s), 100 * mean(is.na(s[[v]]))),
      NA_character_,
      "Proportion not recorded, by sex and group. Adjusted models below are complete-case; the left ventricular ejection fraction enters only as a not-recorded indicator for that reason.",
      measure = "proportion not recorded (%)", sex = g,
      group = ifelse(pc == 1, "CTO PCI performed", "CTO PCI not performed")))
  }
}

# --- 시행률 및 시행 odds (여 vs 남) ---------------------------------------
for (g in c("Male", "Female")) {
  s <- reg[reg$sex_f == g, , drop = FALSE]
  ci <- stats::binom.test(sum(s$pci), nrow(s))$conf.int
  addA1(rw("A-1(rate)", "CTO PCI performed",
    sprintf("Crude rate, %s", ifelse(g == "Male", "men", "women")), nrow(s),
    ifelse(g == "Female", nrow(s), 0L), 100 * mean(s$pci), 100 * ci[1], 100 * ci[2],
    sprintf("%d/%d (%.1f%%, 95%% CI %.1f-%.1f)", sum(s$pci), nrow(s), 100 * mean(s$pci),
            100 * ci[1], 100 * ci[2]), NA_character_,
    "Crude proportion in whom CTO PCI was performed, registry 2014-2023, with an exact binomial 95% CI.",
    measure = "proportion (%)", sex = g))
}
reg$lvef_missing <- as.integer(is.na(reg$preEF))
PCI_SPECS <- list(
  list(cv = character(0), lab = "unadjusted"),
  list(cv = "Age",        lab = "age-adjusted"),
  list(cv = c("Age", "DM", "CKD_2", "preCABG", "lvef_missing", "Clindx_f", "JCTO_re"),
       lab = "age + diabetes, CKD, previous CABG, LVEF-not-recorded indicator, presentation, J-CTO score"),
  list(cv = c("Age", "DM", "CKD_2", "preCABG", "lvef40", "Clindx_f", "JCTO_re"),
       lab = "sensitivity: same set with LVEF <40% in place of the not-recorded indicator (complete case)"))
for (sp in PCI_SPECS) {
  d <- reg
  if (length(sp$cv)) d <- droplevels(d[stats::complete.cases(d[, sp$cv, drop = FALSE]), , drop = FALSE])
  f <- stats::as.formula(paste("pci ~ female", if (length(sp$cv)) paste("+", paste(sp$cv, collapse = " + ")) else ""))
  m <- stats::glm(f, data = d, family = stats::binomial())
  ct <- lmtest::coeftest(m, vcov. = sandwich::vcovCL(m, cluster = d$center_anon))
  b <- ct["female", "Estimate"]; se <- ct["female", "Std. Error"]
  addA1(rw("A-1(odds)", "CTO PCI performed (women vs men)",
    sprintf("Logistic regression, %s, centre-clustered robust SE", sp$lab),
    nrow(d), sum(d$female == 1), exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se),
    fmt2(exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)), fmt_p(ct["female", "Pr(>|z|)"]),
    paste0("Odds of having undergone CTO PCI for women versus men in the registry population of ",
           "2014-2023. Complete case for the covariates listed; the left ventricular ejection ",
           "fraction is represented by a not-recorded indicator in the primary adjusted model ",
           "because it is missing in about one fifth of patients. Centres contributing only ",
           "treated patients remain in the clustering structure. This is a description of who was ",
           "treated, not an estimate of an effect of sex on the treatment decision. ", COMP_NOTE),
    measure = "odds ratio (women vs men)",
    n_events = sum(d$pci == 1), n_parameters = length(stats::coef(m))))
}

# --- 맥락용: 미시행군의 1년 전사망 (기술 목적 한정) -----------------------
np2 <- np; np2$.t <- np2$death_fudur; np2$.e <- np2$Death
np2$.t[np2$.e == 1 & np2$.t == 0] <- 0.5
fit <- survival::survfit(Surv(.t, .e) ~ sex_f, data = np2, conf.type = "log-log")
sm <- summary(fit, times = 365, extend = TRUE)
lr <- survival::survdiff(Surv(.t, .e) ~ sex_f, data = np2)
p_lr <- stats::pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
for (g in c("Male", "Female")) {
  i <- which(as.character(sm$strata) == paste0("sex_f=", g))
  addA1(rw("A-1(context)", "All-cause death at 1 year, CTO PCI not performed",
    sprintf("Kaplan-Meier, %s", ifelse(g == "Male", "men", "women")),
    sum(np2$sex_f == g), sum(np2$sex_f == "Female"),
    100 * (1 - sm$surv[i]), 100 * (1 - sm$upper[i]), 100 * (1 - sm$lower[i]),
    sprintf("%.1f%% (%.1f-%.1f)", 100 * (1 - sm$surv[i]), 100 * (1 - sm$upper[i]),
            100 * (1 - sm$lower[i])), fmt_p(p_lr),
    paste0("DESCRIPTIVE CONTEXT ONLY, not an adjusted comparison and not a treatment comparison: ",
           "1-year Kaplan-Meier mortality within the group in whom CTO PCI was not performed, by ",
           "sex (log-rank P over the whole follow-up). Indication, frailty and the reason for not ",
           "proceeding are not recorded, so this must not be read as an effect of sex or of ",
           "treatment. ", COMP_NOTE),
    measure = "cumulative incidence (%)", sex = g,
    n_events = sum(np2$.e[np2$sex_f == g] == 1)))
}

A1 <- dplyr::bind_rows(a1)
write.csv(A1, file.path(DIR_TAB, "table35_A1_registry_pci_vs_nopci_by_sex.csv"),
          row.names = FALSE, na = "")
logmsg("[A-1] wrote %d rows to table35_A1_registry_pci_vs_nopci_by_sex.csv", nrow(A1))

###############################################################################
# A-4 (R2-M4) — bilateral injection 의 "적응 가능 상황" 제한 + 자원·안전성
#
# 대측 공여혈관(contralateral donor) 정의: CTO 가 RCA 면 좌측(LAD 또는 LCX) 공여,
#   CTO 가 LAD·LCX·LM 이면 RCA 공여. 결측 0%.
#   검증: 대측 존재군의 bilateral 시행률 남 39.0% / 여 29.9% — 사전 감사치와 일치.
###############################################################################
message("\n================ A-4: bilateral injection, indication-restricted ================")

a4 <- list(); addA4 <- function(x) if (!is.null(x)) a4[[length(a4) + 1]] <<- x

co$contra_donor <- with(co, ifelse(CTO_Vessel == 3,
                                   as.integer(lad_donor == 1 | lcx_donor == 1),
                                   as.integer(rca_donor == 1)))
stopifnot(sum(is.na(co$contra_donor)) == 0)
# ⚠️ 코호트는 tibble 이라 rownames 가 위치 문자열이다 — rownames 기반 조인 금지.
#    같은 식을 an 에 직접 적용한다.
an$contra_donor <- with(an, ifelse(CTO_Vessel == 3,
                                   as.integer(lad_donor == 1 | lcx_donor == 1),
                                   as.integer(rca_donor == 1)))
stopifnot(sum(is.na(an$contra_donor)) == 0)
chk <- tapply(co$Bilati_inj[co$contra_donor == 1], co$sex_f[co$contra_donor == 1],
              function(x) 100 * mean(x, na.rm = TRUE))
logmsg("[A-4] contralateral donor present: n = %d (%.1f%%); bilateral injection men %.1f%% / women %.1f%% (audit: 39.0 / 29.9)",
       sum(co$contra_donor == 1), 100 * mean(co$contra_donor == 1), chk[["Male"]], chk[["Female"]])
stopifnot(abs(chk[["Male"]] - 39.0) < 0.1, abs(chk[["Female"]] - 29.9) < 0.1)

# bilateral 전혀 기록되지 않은 센터
mz <- tapply(is.na(co$Bilati_inj), co$center_anon, mean)
CEN_NOREC <- as.integer(names(mz)[mz == 1])
cen_sizes <- table(co$center_anon)
CEN_DROP  <- CEN_NOREC[which.max(cen_sizes[as.character(CEN_NOREC)])]
logmsg("[A-4] centres with bilateral injection never recorded: %s (n = %s); the sensitivity analysis removes centre %d (n = %d)",
       paste(CEN_NOREC, collapse = ", "),
       paste(as.integer(cen_sizes[as.character(CEN_NOREC)]), collapse = ", "),
       CEN_DROP, as.integer(cen_sizes[as.character(CEN_DROP)]))

A4_NOTE <- paste0(
  "Restriction defined on pre-procedural anatomy: a contralateral collateral donor vessel is ",
  "present (left-system donor for a right coronary artery CTO, right coronary artery donor for a ",
  "left-system CTO; recorded in all patients) and the Rentrop collateral grade is at least 1. ",
  "The restriction is a judgement about when bilateral injection is anatomically available, not a ",
  "formal indication; the indication itself, collateral suitability, operator rationale, access ",
  "strategy and contrast or renal concerns are not recorded.")

bil_models <- function(d, tag, note = "") {
  out <- list()
  dd <- droplevels(d[!is.na(d$Bilati_inj), , drop = FALSE])
  pm <- 100 * mean(dd$Bilati_inj[dd$female == 0]); pw <- 100 * mean(dd$Bilati_inj[dd$female == 1])
  out[[1]] <- rw("A-4", "Bilateral contrast injection", sprintf("%s | crude rate by sex", tag),
    nrow(dd), sum(dd$female == 1), pw - pm, NA_real_, NA_real_,
    sprintf("men %d/%d (%.1f%%) vs women %d/%d (%.1f%%)",
            sum(dd$Bilati_inj[dd$female == 0]), sum(dd$female == 0), pm,
            sum(dd$Bilati_inj[dd$female == 1]), sum(dd$female == 1), pw),
    NA_character_, paste0(A4_NOTE, " ", note),
    measure = "crude proportion (%), men vs women",
    pct_men = round(pm, 1), pct_women = round(pw, 1))
  # (a) 센터 랜덤절편 (table18 의 주 사양)
  m <- tryCatch(lme4::glmer(stats::as.formula(paste0(
        "Bilati_inj ~ female + ", paste(COV_M4_NOCENTER, collapse = " + "), " + (1 | center_anon)")),
        data = dd, family = stats::binomial(), nAGQ = 0L,
        control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5))),
        error = function(e) NULL)
  if (!is.null(m)) {
    s <- summary(m)$coefficients
    b <- s["female", "Estimate"]; se <- s["female", "Std. Error"]
    sig2 <- as.numeric(lme4::VarCorr(m)$center_anon[1]); icc <- sig2 / (sig2 + pi^2 / 3)
    cvg <- m@optinfo$conv$lme4$messages
    out[[length(out) + 1]] <- rw("A-4", "Bilateral contrast injection",
      sprintf("%s | Model 4 covariates + centre random intercept", tag),
      nrow(dd), sum(dd$female == 1), exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se),
      fmt2(exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)), fmt_p(s["female", "Pr(>|z|)"]),
      paste0(A4_NOTE, " ", note, " Covariates: Model 4 minus centre volume (absorbed by the centre ",
             "term), the specification of table18; fast approximation (nAGQ = 0).",
             if (length(cvg)) paste0(" CONVERGENCE NOTE: ", paste(cvg, collapse = "; ")) else ""),
      measure = "adjusted odds ratio (women vs men)",
      re_variance = round(sig2, 4), ICC = round(icc, 4),
      n_centres = length(unique(dd$center_anon)))
  }
  # (b) Model 4 공변량 + 센터 clustered robust SE
  m2 <- stats::glm(stats::as.formula(paste0("Bilati_inj ~ female + ", paste(COV_M4, collapse = " + "))),
                   data = dd, family = stats::binomial())
  ct <- lmtest::coeftest(m2, vcov. = sandwich::vcovCL(m2, cluster = dd$center_anon))
  b <- ct["female", "Estimate"]; se <- ct["female", "Std. Error"]
  out[[length(out) + 1]] <- rw("A-4", "Bilateral contrast injection",
    sprintf("%s | Model 4 covariates, centre-clustered robust SE", tag),
    nrow(dd), sum(dd$female == 1), exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se),
    fmt2(exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)), fmt_p(ct["female", "Pr(>|z|)"]),
    paste0(A4_NOTE, " ", note), measure = "adjusted odds ratio (women vs men)")
  # (c) 센터 고정효과 + clustered robust SE
  m3 <- stats::glm(stats::as.formula(paste0("Bilati_inj ~ female + ",
        paste(c(COV_M4_NOCENTER, "factor(center_anon)"), collapse = " + "))),
        data = dd, family = stats::binomial())
  ct3 <- lmtest::coeftest(m3, vcov. = sandwich::vcovCL(m3, cluster = dd$center_anon))
  b <- ct3["female", "Estimate"]; se <- ct3["female", "Std. Error"]
  out[[length(out) + 1]] <- rw("A-4", "Bilateral contrast injection",
    sprintf("%s | centre fixed effects, centre-clustered robust SE", tag),
    nrow(dd), sum(dd$female == 1), exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se),
    fmt2(exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)), fmt_p(ct3["female", "Pr(>|z|)"]),
    paste0(A4_NOTE, " ", note, " Within-centre contrast with centre indicators in the mean structure."),
    measure = "adjusted odds ratio (women vs men)")
  dplyr::bind_rows(out)
}

d_i   <- an[an$contra_donor == 1 & !is.na(an$Rentrop) & an$Rentrop >= 1, , drop = FALSE]
d_ii  <- an[an$contra_donor == 1, , drop = FALSE]
d_iii <- d_i[d_i$center_anon != CEN_DROP, , drop = FALSE]
logmsg("[A-4] restriction (i) n = %d (women %d) | (ii) n = %d (women %d) | (iii) n = %d (women %d)",
       nrow(d_i), sum(d_i$female), nrow(d_ii), sum(d_ii$female), nrow(d_iii), sum(d_iii$female))
addA4(bil_models(d_i,  "(i) PRIMARY restriction: contralateral donor present and Rentrop grade >=1"))
addA4(bil_models(d_ii, "(ii) sensitivity: contralateral donor present (any Rentrop grade)"))
addA4(bil_models(d_iii, sprintf("(iii) sensitivity: restriction (i) excluding centre %d", CEN_DROP),
  sprintf(paste0("Centre %d (n = %d in the cohort) never recorded bilateral injection, so all of ",
                 "its patients are already excluded from every bilateral model by complete-case ",
                 "estimation; this row is expected to be numerically identical to (i) and is ",
                 "reported to document that."), CEN_DROP, as.integer(cen_sizes[as.character(CEN_DROP)]))))

# --- (iv) 자원 사용 · 신장 안전성 ------------------------------------------
NOT_REC <- paste0("Procedure duration, radiation dose, vascular access site and dialysis are not ",
                  "recorded in the registry and cannot be compared.")
cen_ok <- function(v) { m <- tapply(is.na(co[[v]]), co$center_anon, mean)
                        as.integer(names(m)[m < 0.5]) }

res_models <- function(v, lab, unit, dig = 0) {
  keep <- cen_ok(v)
  d <- droplevels(an[an$center_anon %in% keep & !is.na(an[[v]]) & an[[v]] > 0, , drop = FALSE])
  out <- list()
  for (g in c("Male", "Female")) {
    x <- d[[v]][d$sex_f == g]
    out[[length(out) + 1]] <- rw("A-4(iv)", lab,
      sprintf("Median [IQR], %s", ifelse(g == "Male", "men", "women")),
      length(x), sum(d$female == 1), stats::median(x), unname(stats::quantile(x, .25)),
      unname(stats::quantile(x, .75)),
      sprintf(paste0("%.", dig, "f [%.", dig, "f-%.", dig, "f] %s"), stats::median(x),
              stats::quantile(x, .25), stats::quantile(x, .75), unit), NA_character_,
      sprintf(paste0("Restricted to the %d of 33 centres that record this variable in more than ",
                     "half of their patients (%d patients of the Model 4 complete-case sample; ",
                     "values as recorded). %s"), length(keep), nrow(d), NOT_REC),
      measure = sprintf("median [IQR], %s", unit), sex = g, n_centres = length(keep))
  }
  p <- stats::wilcox.test(d[[v]] ~ d$female)$p.value
  m <- lme4::lmer(stats::as.formula(paste0("log(", v, ") ~ female + ",
        paste(COV_M4_NOCENTER, collapse = " + "), " + (1 | center_anon)")), data = d, REML = TRUE)
  s <- summary(m)$coefficients
  b <- s["female", "Estimate"]; se <- s["female", "Std. Error"]
  pw <- 2 * stats::pnorm(-abs(b / se))
  sig2 <- as.numeric(lme4::VarCorr(m)$center_anon[1])
  out[[length(out) + 1]] <- rw("A-4(iv)", lab,
    "Adjusted ratio of geometric means (women / men), log-linear mixed model with centre random intercept",
    nrow(d), sum(d$female == 1), exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se),
    sprintf("%.3f (%.3f-%.3f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)), fmt_p(pw),
    sprintf(paste0("The outcome is log-transformed because both measures are right-skewed, so the ",
                   "exponentiated sex coefficient is a multiplicative difference: a value below 1 ",
                   "means lower use in women. Model 4 covariates minus centre volume plus a centre ",
                   "random intercept (variance %.3f). Crude Mann-Whitney P = %s. Restricted to the ",
                   "%d centres recording the variable in more than half of their patients. %s"),
            sig2, fmt_p(p), length(keep), NOT_REC),
    measure = "ratio of adjusted geometric means (women / men)",
    re_variance = round(sig2, 4), n_centres = length(keep), crude_p_fmt = fmt_p(p))
  dplyr::bind_rows(out)
}
addA4(res_models("dye", "Contrast media volume", "mL"))
addA4(res_models("fluorotime", "Fluoroscopy time", "min"))
# 브리프 외 추가(기술 목적): 체중당 조영제량 — "여성이 작아서" 라는 당연한 반론에
# 미리 답하기 위한 행. 주 모형은 BMI 를 공변량으로 포함하나 체중 자체는 포함하지
# 않는다(positivity 위반으로 주 모형에서 제외, supplement §1).
an$dye_per_kg <- ifelse(is.na(an$dye) | is.na(an$Weight) | an$Weight <= 0,
                        NA_real_, an$dye / an$Weight)
co$dye_per_kg <- ifelse(is.na(co$dye) | is.na(co$Weight) | co$Weight <= 0,
                        NA_real_, co$dye / co$Weight)
addA4(res_models("dye_per_kg", "Contrast media volume per kilogram of body weight", "mL/kg", dig = 2))

# 신장 안전성: CIN (출판값 재현) + 기록 센터 제한 + 랜덤절편 + Cr 상승 정의
mk_cr_rise <- function(d) ifelse(is.na(d$cr_increase) & is.na(d$cr_change), NA_real_,
                                 as.numeric((!is.na(d$cr_increase) & d$cr_increase >= 0.5) |
                                            (!is.na(d$cr_change) & d$cr_change >= 1.25)))
an$cr_rise <- mk_cr_rise(an); co$cr_rise <- mk_cr_rise(co)
disc <- sum(an$CIN_2 != an$cr_rise, na.rm = TRUE)
logmsg("[A-4(iv)] CIN_2 vs (creatinine rise >=0.5 mg/dL or >=25%%): %d discordant of %d jointly recorded",
       disc, sum(!is.na(an$CIN_2) & !is.na(an$cr_rise)))

renal_models <- function(v, lab, extra = "") {
  out <- list()
  specs <- list(
    list(d = an, tag = "all centres, Model 4 covariates, centre-clustered robust SE", re = FALSE),
    list(d = an[an$center_anon %in% cen_ok(v), , drop = FALSE],
         tag = sprintf("centres recording the variable in >50%% of patients, Model 4 covariates, centre-clustered robust SE"), re = FALSE),
    list(d = an[an$center_anon %in% cen_ok(v), , drop = FALSE],
         tag = "centres recording the variable in >50% of patients, Model 4 covariates + centre random intercept", re = TRUE))
  for (sp in specs) {
    d <- droplevels(sp$d[!is.na(sp$d[[v]]), , drop = FALSE]); d$.y <- d[[v]]
    pm <- 100 * mean(d$.y[d$female == 0]); pw <- 100 * mean(d$.y[d$female == 1])
    if (!sp$re) {
      m <- stats::glm(stats::as.formula(paste0(".y ~ female + ", paste(COV_M4, collapse = " + "))),
                      data = d, family = stats::binomial())
      ct <- lmtest::coeftest(m, vcov. = sandwich::vcovCL(m, cluster = d$center_anon))
      b <- ct["female", "Estimate"]; se <- ct["female", "Std. Error"]; pv <- ct["female", "Pr(>|z|)"]
      icc <- NA_real_
    } else {
      m <- lme4::glmer(stats::as.formula(paste0(".y ~ female + ",
            paste(COV_M4_NOCENTER, collapse = " + "), " + (1 | center_anon)")),
            data = d, family = stats::binomial(), nAGQ = 0L,
            control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
      s <- summary(m)$coefficients
      b <- s["female", "Estimate"]; se <- s["female", "Std. Error"]; pv <- s["female", "Pr(>|z|)"]
      sig2 <- as.numeric(lme4::VarCorr(m)$center_anon[1]); icc <- sig2 / (sig2 + pi^2 / 3)
    }
    out[[length(out) + 1]] <- rw("A-4(iv)", lab, sp$tag, nrow(d), sum(d$female == 1),
      exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se),
      fmt2(exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)), fmt_p(pv),
      paste0("Adjusted odds ratio for women versus men. Crude rates: men ", sprintf("%.1f%%", pm),
             " versus women ", sprintf("%.1f%%", pw), ". ", extra, " ", NOT_REC),
      measure = "adjusted odds ratio (women vs men)",
      n_events = sum(d$.y == 1), pct_men = round(pm, 1), pct_women = round(pw, 1),
      ICC = if (is.na(icc)) NA_real_ else round(icc, 4))
  }
  dplyr::bind_rows(out)
}
addA4(renal_models("CIN_2", "Contrast-induced nephropathy",
  "The first row reproduces the estimate reported in the submitted manuscript (1.27, 0.90-1.77; n = 5,122, 285 events)."))
addA4(renal_models("cr_rise",
  "Creatinine rise >=0.5 mg/dL or >=25% at day 2 (registry nephropathy criteria, recomputed)",
  sprintf(paste0("Recomputed from the day-2 creatinine fields rather than read from the recorded ",
                 "nephropathy flag; the two disagree in %d of the jointly recorded patients, so ",
                 "both are reported."), disc)))

A4 <- dplyr::bind_rows(a4)
write.csv(A4, file.path(DIR_TAB, "table36_A4_bilateral_indication_and_resource.csv"),
          row.names = FALSE, na = "")
logmsg("[A-4] wrote %d rows to table36_A4_bilateral_indication_and_resource.csv", nrow(A4))

###############################################################################
# A-2③ (R2-M2) — delta-adjusted (tipping-point) 민감도, 0-1년 창
#
# 질문: person-time 을 기여하지 않는 환자(추적 미제출)의 1년 사건 위험이 성별로
#       달랐다면 0-1년 창 추정치가 어디까지 움직이는가?
# 절차: ① 관측자료(창 내 person-time 기여자)로 Cox 모형을 적합
#       ② 미기여 환자 각자의 모형 예측 1년 위험 r_i 에 delta(성별)를 곱해
#          p_i = min(0.999, delta_sex * r_i) 로 사건 여부를 추출하고, 사건이면
#          적합된 기저위험의 (0,365일] 조건부 분포에서 시각을 역변환 추출
#       ③ 20회 대치마다 같은 모형을 재적합하고 Rubin 규칙으로 결합
# 격자: delta_women ∈ {0.5,0.75,1,1.25,1.5,2} (delta_men=1) 및 그 거울 격자
###############################################################################
message("\n================ A-2(3): delta-adjusted tipping-point sensitivity ================")

a2c <- list(); addA2c <- function(x) a2c[[length(a2c) + 1]] <<- x

TIP_NOTE <- paste0(
  "Delta-adjusted (tipping-point) sensitivity analysis for patients who contribute no person-time ",
  "to the first-year model because no follow-up was submitted for them. Their first-year outcome ",
  "is imputed ", M_IMP, " times from the risk predicted for them by the model fitted to the ",
  "patients with follow-up, multiplied by a sex-specific factor delta (delta = 1 reproduces the ",
  "model-predicted risk, delta = 2 doubles it); estimates are combined by Rubin's rules, and ",
  "all delta values reuse the same random draws so that differences across the grid come only ",
  "from delta. The ",
  "analysis asks how far the first-year estimate can be moved by an assumed sex-specific ",
  "difference in the unobserved risk; it is not an estimate of that difference.")

rubin <- function(q, se) {
  m <- length(q); qbar <- mean(q); Ubar <- mean(se^2); B <- stats::var(q)
  if (!is.finite(B) || B <= 0) B <- 0
  T <- Ubar + (1 + 1 / m) * B
  df <- if (B > 0) (m - 1) * (1 + Ubar / ((1 + 1 / m) * B))^2 else Inf
  tq <- if (is.finite(df)) stats::qt(.975, df) else 1.96
  p  <- if (is.finite(df)) 2 * stats::pt(-abs(qbar) / sqrt(T), df)
        else 2 * stats::pnorm(-abs(qbar) / sqrt(T))
  c(est = qbar, lo = qbar - tq * sqrt(T), hi = qbar + tq * sqrt(T), p = p, df = df)
}
interp_cross <- function(x, y) {           # y 가 0 을 지나는 x (선형보간)
  o <- order(x); x <- x[o]; y <- y[o]
  for (i in seq_len(length(x) - 1)) {
    if (is.finite(y[i]) && is.finite(y[i + 1]) && y[i] * y[i + 1] <= 0 && y[i] != y[i + 1])
      return(x[i] + (0 - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]))
  }
  NA_real_
}

tip_one <- function(key, covs, cov_lab) {
  d_all <- make_surv(an, key)
  vars <- intersect(unique(c("female", "center_anon", covs)), names(d_all))
  d_all <- droplevels(d_all[stats::complete.cases(d_all[, vars, drop = FALSE]), , drop = FALSE])
  d_all$.wt <- pmin(d_all$.time, 365)
  d_all$.we <- ifelse(d_all$.time <= 365 & d_all$.event == 1, 1L, 0L)
  obs <- d_all[d_all$.time > 0, , drop = FALSE]
  mis <- d_all[d_all$.time == 0, , drop = FALSE]
  obs$.time <- obs$.wt; obs$.event <- obs$.we
  rhs <- if (length(covs)) paste(c("female", covs), collapse = " + ") else "female"
  f <- stats::as.formula(paste0("Surv(.time, .event) ~ ", rhs))
  # 관측자료 적합 (추정치용: table22 와 동일한 cox_model) + 기저위험용 Breslow 적합
  m_rob <- cox_model(obs, covs)
  m_bre <- survival::coxph(f, data = obs, ties = "breslow")
  b <- stats::coef(m_bre); b[is.na(b)] <- 0
  bh <- survival::basehaz(m_bre, centered = FALSE)
  bh <- bh[bh$time <= 365, , drop = FALSE]
  H <- bh$hazard; TT <- bh$time; H365 <- if (length(H)) max(H) else 0
  # 설계행렬은 전체 창 표본에서 한 번 만든다 (대치 대상만으로 만들면 특정 factor
  # 수준이 빠져 열이 사라질 수 있다)
  Xall <- stats::model.matrix(f, data = d_all)[, -1, drop = FALSE]
  Xm <- Xall[which(d_all$.time == 0), names(b), drop = FALSE]
  lp <- as.vector(Xm %*% b)
  r  <- 1 - exp(-H365 * exp(lp))                         # 모형 예측 1년 위험
  logmsg("[A-2(3)] %s / %s: observed n = %d (events %d), imputed n = %d (women %d); predicted 1-year risk mean %.3f (men %.3f / women %.3f)",
         key, cov_lab, nrow(obs), sum(obs$.event), nrow(mis), sum(mis$female == 1),
         mean(r), mean(r[mis$female == 0]), mean(r[mis$female == 1]))
  GRID <- rbind(data.frame(dw = c(0.5, 0.75, 1, 1.25, 1.5, 2), dm = 1),
                data.frame(dw = 1, dm = c(0.5, 0.75, 1.25, 1.5, 2)))
  # 공통 난수: 모든 delta 칸이 같은 균등난수를 쓰므로 칸 사이 차이는 delta 에서만
  # 오고(몬테카를로 잡음 제거) tipping-point 보간이 단조로워진다.
  set.seed(SEED)
  U <- matrix(stats::runif(M_IMP * nrow(mis)), M_IMP, nrow(mis))   # 사건 여부용
  V <- matrix(stats::runif(M_IMP * nrow(mis)), M_IMP, nrow(mis))   # 사건 시각용
  out <- list()
  for (gi in seq_len(nrow(GRID))) {
    dw <- GRID$dw[gi]; dm <- GRID$dm[gi]
    delta <- ifelse(mis$female == 1, dw, dm)
    p_i <- pmin(0.999, delta * r)
    qs <- numeric(M_IMP); ses <- numeric(M_IMP)
    for (k in seq_len(M_IMP)) {
      ev <- as.integer(U[k, ] < p_i)
      tm <- rep(365, length(p_i))
      if (any(ev == 1) && length(TT)) {
        for (i in which(ev == 1)) {
          Fi <- 1 - exp(-H * exp(lp[i]))
          tgt <- V[k, i] * Fi[length(Fi)]
          j <- which(Fi >= tgt)[1]
          tm[i] <- if (is.na(j)) 365 else max(TT[j], 0.5)
        }
      }
      di <- mis; di$.time <- tm; di$.event <- ev
      dc <- rbind(obs[, names(di), drop = FALSE], di)
      mk <- tryCatch(cox_model(dc, covs), error = function(e) NULL)
      if (is.null(mk)) { qs[k] <- NA; ses[k] <- NA; next }
      sk <- summary(mk)$coefficients
      qs[k] <- sk["female", "coef"]; ses[k] <- sk["female", "robust se"]
    }
    ok <- !is.na(qs); rb <- rubin(qs[ok], ses[ok])
    out[[length(out) + 1]] <- rw("A-2(3)", OUT[[key]]$lab,
      sprintf("0-1 year, %s, delta(women) = %.2f / delta(men) = %.2f", cov_lab, dw, dm),
      nrow(obs) + nrow(mis), sum(d_all$female == 1),
      exp(rb[["est"]]), exp(rb[["lo"]]), exp(rb[["hi"]]),
      fmt2(exp(rb[["est"]]), exp(rb[["lo"]]), exp(rb[["hi"]])), fmt_p(rb[["p"]]),
      sprintf("%s Imputations combined: %d of %d. Observed n = %d with %d events; imputed n = %d (women %d).",
              TIP_NOTE, sum(ok), M_IMP, nrow(obs), sum(obs$.event), nrow(mis), sum(mis$female == 1)),
      measure = "hazard ratio (women vs men), first year",
      delta_women = dw, delta_men = dm, HR = exp(rb[["est"]]), p_value = rb[["p"]],
      n_observed = nrow(obs), n_events_observed = sum(obs$.event), n_imputed = nrow(mis),
      n_imputed_women = sum(mis$female == 1), rubin_df = round(rb[["df"]], 1))
  }
  tab <- dplyr::bind_rows(out)
  # tipping point (delta_women 격자에서 보간)
  w <- tab[tab$delta_men == 1, , drop = FALSE]
  x_hr <- interp_cross(w$delta_women, log(w$HR))
  x_lo <- interp_cross(w$delta_women, log(w$lo))
  tab <- dplyr::bind_rows(tab, rw("A-2(3)", OUT[[key]]$lab,
    sprintf("0-1 year, %s, tipping point in delta(women) (linear interpolation on the log scale)", cov_lab),
    nrow(obs) + nrow(mis), sum(d_all$female == 1), x_hr, NA_real_, NA_real_,
    if (is.na(x_hr)) sprintf("point estimate does not reach 1.0 for delta(women) up to 2.0 (highest grid value: %s)",
                             w$est_fmt[which.max(w$delta_women)])
    else sprintf("delta(women) = %.2f for the point estimate to reach 1.0%s", x_hr,
                 if (is.na(x_lo)) "; the lower confidence limit does not reach 1.0 within the grid"
                 else sprintf("; delta(women) = %.2f for the lower confidence limit to exceed 1.0", x_lo)),
    NA_character_,
    paste0(TIP_NOTE, " The tipping point is the factor by which the unobserved first-year risk of ",
           "the women without follow-up would have to be multiplied, relative to the risk the ",
           "model predicts for them, before the first-year estimate reaches the null."),
    measure = "tipping-point delta(women)",
    delta_women = x_hr, delta_men = 1,
    HR = 1, p_value = NA_real_, n_imputed = nrow(mis), n_imputed_women = sum(mis$female == 1),
    tipping_delta_point_estimate = x_hr, tipping_delta_lower_CI = x_lo))
  tab
}

for (key in c("MACCE_pub", "DMS")) {
  for (sp in list(list(cv = COV_PARS, lb = "reduced (age, DM, CKD, AF)"),
                  list(cv = COV_M4,   lb = "Model 4"))) {
    message(sprintf("  [A-2(3)] %s / %s ...", key, sp$lb))
    addA2c(tip_one(key, sp$cv, sp$lb))
  }
}
A2C <- dplyr::bind_rows(a2c)
# 요청 규격의 tidy 컬럼: analysis, outcome, model, delta_women, delta_men, HR, lo, hi, p
# (⚠️ rw() 의 ... 는 formal 뒤에 오지 않으므로 extra 이름 "p" 는 formal "p_fmt" 로
#  부분일치해 버린다 — 그래서 p_value 로 받고 여기서 p 컬럼을 명시 생성한다.)
stopifnot(is.character(A2C$p_fmt))
A2C <- data.frame(A2C, p = A2C$p_value, check.names = FALSE, stringsAsFactors = FALSE)
KEY <- c("analysis", "outcome", "model", "delta_women", "delta_men", "HR", "lo", "hi", "p")
A2C <- A2C[, c(KEY, setdiff(names(A2C), KEY))]
write.csv(A2C, file.path(DIR_TAB, "table37_A2c_delta_adjusted_tipping_point.csv"),
          row.names = FALSE, na = "")
logmsg("[A-2(3)] wrote %d rows to table37_A2c_delta_adjusted_tipping_point.csv", nrow(A2C))

###############################################################################
# 로그 · sessionInfo
###############################################################################
writeLines(c(sprintf("18_phase7_reviewer_analyses_b.R  run %s", format(Sys.time())),
             sprintf("seed = %d | imputations per delta cell = %d", SEED, M_IMP), "", LOG),
           file.path(DIR_REP, "phase7b_analyses_log.txt"))
utils::capture.output(utils::sessionInfo(), file = file.path(DIR_REP, "phase7b_sessionInfo.txt"))
message("\n== done ==")
