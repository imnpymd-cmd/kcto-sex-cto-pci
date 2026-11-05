###############################################################################
# 05_qa_sensitivity.R — Phase 2.5 QA 요구 추가 분석 (A1–A10)
#
# 실행:  Rscript scripts/05_qa_sensitivity.R
# 입력:  data/processed/kcto_sex_analysis.rds  (via scripts/00_cohort_prep.R)
#        output/tables/table6f_*.csv, table7_*.csv  (04_primary_analysis.R 산출물)
# 출력:  output/tables/table14* .. table20*, tableS12* .. tableS15*
#        interim_report/phase2_stats/reports/qa_sensitivity_05_*.txt
#
# 코호트: flag_cto_pci & index_yr 2014-2023 -> N = 6,777 (CLAUDE.md §1)
#         ⛔ flag_analysis_primary 사용 금지
# Primary outcome: DMS (death + MI + stroke)
#
# ⚠️ 본 분석은 관찰연구다. 모든 추정치는 연관성(association)이며 인과 해석 금지.
#    "explained by" 류의 매개 함의 표현도 사용하지 않는다.
# ⛔ 사전 계획을 시사하는 표현 금지 (CLAUDE.md §1-2 금지어 목록) (CLAUDE.md §1-2).
#    본 분석은 사전등록되지 않았다. 아래 검정 family 는 "본 분석에서 정의한 것"이다.
#
# 목적 — Phase 2 의 치료 과정 격차 5개 소견에 대한 네 가지 대안 설명을 검정한다:
#   (1) 센터 구성 artifact 인가?                   -> A1
#   (2) 다중비교로 인한 우연인가?                   -> A2
#   (3) 여성의 작은 체격 때문에 정당한 선택인가?     -> A3
#   (4) prasugrel 라벨 제한이 만든 인위적 격차인가?  -> A4
#   (+) 에스컬레이션 적응증 자체가 적었던 것인가?    -> A5
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(lme4)
})

set.seed(42)

PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "00_cohort_prep.R"))

DIR_TAB  <- file.path(PROJ, "output", "tables")
DIR_REPX <- file.path(PROJ, "interim_report", "phase2_stats", "reports")
dir.create(DIR_TAB,  recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_REPX, recursive = TRUE, showWarnings = FALSE)

co <- build_cohort()

# 주 분석 표본 (Model 4 complete case) — 04_primary_analysis.R 과 동일
an <- co[stats::complete.cases(co[, COV_M4, drop = FALSE]), , drop = FALSE]
message(sprintf("[analysis set] Model 4 complete case: N = %d (F = %d)", nrow(an), sum(an$female)))
stopifnot(nrow(co) == 6777)

LOG <- list()
logmsg <- function(...) {
  s <- sprintf(...)
  LOG[[length(LOG) + 1]] <<- s
  message(s)
}

###############################################################################
# 0. 공통 helper — 표준화된 결과행 · 모델 fitter
###############################################################################

.counts <- function(d, yvar) {
  y <- d[[yvar]]
  list(n = nrow(d), n_events = sum(y == 1, na.rm = TRUE),
       n_male = sum(d$female == 0), n_female = sum(d$female == 1),
       events_male   = sum(y == 1 & d$female == 0, na.rm = TRUE),
       events_female = sum(y == 1 & d$female == 1, na.rm = TRUE))
}

std_row <- function(analysis, outcome, model, cnt, measure, est, lo, hi, p, note = "") {
  data.frame(
    analysis = analysis, outcome = outcome, model = model,
    n = cnt$n, n_events = cnt$n_events,
    n_male = cnt$n_male, n_female = cnt$n_female,
    events_male = cnt$events_male, events_female = cnt$events_female,
    pct_male   = round(100 * cnt$events_male / cnt$n_male, 1),
    pct_female = round(100 * cnt$events_female / cnt$n_female, 1),
    measure = measure,
    estimate = est, lo = lo, hi = hi,
    est_fmt = ifelse(is.na(est), NA_character_, sprintf("%.2f (%.2f-%.2f)", est, lo, hi)),
    p_value = p, p_fmt = fmt_p(p), note = note,
    stringsAsFactors = FALSE, row.names = NULL)
}

# 모형에 실제 투입되는 complete-case 부분집합
.cc_subset <- function(dat, yvar, covs) {
  vars <- unique(c(yvar, "female", "center_anon", covs))
  vars <- intersect(vars, names(dat))
  dat[stats::complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE]
}

# (i) 로지스틱 + cluster-robust SE (센터 평균구조 보정 없음 = 기존 Model 4 방식)
fit_logit_cr <- function(dat, yvar, covs, analysis, outcome_lab, model_lab, note = "") {
  d <- .cc_subset(dat, yvar, covs); d$.y <- d[[yvar]]
  if (nrow(d) < 30 || length(unique(d$.y)) < 2) return(NULL)
  rhs <- if (length(covs)) paste(c("female", covs), collapse = " + ") else "female"
  m <- stats::glm(stats::as.formula(paste0(".y ~ ", rhs)), data = d, family = stats::binomial())
  vc <- tryCatch(sandwich::vcovCL(m, cluster = d$center_anon), error = function(e) stats::vcov(m))
  ct <- lmtest::coeftest(m, vcov. = vc)
  b <- ct["female", "Estimate"]; se <- ct["female", "Std. Error"]
  std_row(analysis, outcome_lab, model_lab, .counts(d, yvar), "OR",
          exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se), ct["female", "Pr(>|z|)"], note)
}

# (ii) 센터 랜덤절편 로지스틱 (within-center; ICC 병기)
#  nagq = 1 : 완전 Laplace 근사 (고정효과 26개를 함께 최적화 -> 지표당 1-5분)
#  nagq = 0 : 고정효과를 PIRLS 로 profile out 한 빠른 근사 (지표당 <1초).
#  A1 에서 두 방식을 같은 표에 병기해 등가성을 명시적으로 보여주고,
#  적합 횟수가 많은 A3/A4/A5 에서는 nagq = 0 을 사용한다(각 행 note 에 표기).
fit_glmer_sex <- function(dat, yvar, covs, analysis, outcome_lab, model_lab, note = "",
                          nagq = 0L) {
  d <- .cc_subset(dat, yvar, covs); d$.y <- d[[yvar]]
  if (nrow(d) < 30 || length(unique(d$.y)) < 2) return(NULL)
  if (length(unique(d$center_anon)) < 3) return(NULL)
  rhs <- if (length(covs)) paste(c("female", covs), collapse = " + ") else "female"
  f <- stats::as.formula(paste0(".y ~ ", rhs, " + (1 | center_anon)"))
  m <- tryCatch(lme4::glmer(f, data = d, family = stats::binomial(), nAGQ = nagq,
                            control = lme4::glmerControl(optimizer = "bobyqa",
                                                         optCtrl = list(maxfun = 3e5))),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  s <- summary(m)
  if (!"female" %in% rownames(s$coefficients)) return(NULL)
  b <- s$coefficients["female", "Estimate"]; se <- s$coefficients["female", "Std. Error"]
  sig2 <- as.numeric(lme4::VarCorr(m)$center_anon[1])
  icc  <- sig2 / (sig2 + pi^2 / 3)
  cvg <- m@optinfo$conv$lme4$messages
  nte <- paste0(note,
                sprintf(" | centre random-intercept variance = %.3f, ICC = %.3f, centres = %d, nAGQ = %d",
                        sig2, icc, length(unique(d$center_anon)), nagq),
                if (length(cvg)) paste0(" | CONVERGENCE NOTE: ", paste(cvg, collapse = "; ")) else "")
  r <- std_row(analysis, outcome_lab, model_lab, .counts(d, yvar), "OR",
               exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se), s$coefficients["female", "Pr(>|z|)"], nte)
  r$re_variance <- round(sig2, 4); r$ICC <- round(icc, 4)
  r
}

# (iii) 센터 고정효과 로지스틱 (within-center), cluster-robust SE
#      결과 변이가 없는 센터(전원 0 / 전원 1)와 한 성별만 있는 센터는 within-center
#      정보를 제공하지 않으므로 제외한다(조건부 로지스틱과 동일한 취급).
fit_fe_logit <- function(dat, yvar, covs, analysis, outcome_lab, model_lab, note = "") {
  d <- .cc_subset(dat, yvar, covs); d$.y <- d[[yvar]]
  ny  <- tapply(d$.y, d$center_anon, function(z) length(unique(z)))
  nsx <- tapply(d$female, d$center_anon, function(z) length(unique(z)))
  keep_ctr <- names(ny)[ny > 1 & nsx > 1]
  n_drop <- length(unique(d$center_anon)) - length(keep_ctr)
  d <- d[as.character(d$center_anon) %in% keep_ctr, , drop = FALSE]
  if (nrow(d) < 30 || length(unique(d$.y)) < 2) return(NULL)
  d$center_f <- factor(d$center_anon)
  rhs <- if (length(covs)) paste(c("female", covs), collapse = " + ") else "female"
  m <- stats::glm(stats::as.formula(paste0(".y ~ ", rhs, " + center_f")),
                  data = d, family = stats::binomial())
  vc <- tryCatch(sandwich::vcovCL(m, cluster = d$center_anon), error = function(e) stats::vcov(m))
  ct <- tryCatch(lmtest::coeftest(m, vcov. = vc), error = function(e) NULL)
  if (is.null(ct) || !"female" %in% rownames(ct)) return(NULL)
  b <- ct["female", "Estimate"]; se <- ct["female", "Std. Error"]
  std_row(analysis, outcome_lab, model_lab, .counts(d, yvar), "OR",
          exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se), ct["female", "Pr(>|z|)"],
          paste0(note, sprintf(" | centers retained = %d (dropped %d with no within-centre variation in outcome or sex)",
                               length(keep_ctr), n_drop)))
}

MGMT <- STRAT_VARS   # 6개 관리전략 지표 (00_cohort_prep.R 정의) = 단일 검정 family

###############################################################################
# A1. 센터 내(within-center) 재추정 — 치료 과정 격차 ★최우선
#
# 왜: 기존 logit_sex() 는 센터를 SE 에만 반영하고 평균 구조에는 Center_volume_f
#     (4범주)만 넣었다. frailty variance 0.661 / 33개 센터로 이질성이 크다.
#     여성이 저-IVUS 센터에 몰려 있으면 소견 전체가 센터 구성 artifact 일 수 있다.
###############################################################################

message("== A1: within-centre re-estimation of management-strategy indicators ==")

a1_rows <- list()
for (z in MGMT) {
  v <- z[1]; lab <- z[2]
  t0 <- Sys.time()
  a1_rows[[length(a1_rows) + 1]] <- fit_logit_cr(
    an, v, character(0), "A1", lab, "(0) Unadjusted, cluster-robust SE",
    "centre used for SE only")
  a1_rows[[length(a1_rows) + 1]] <- fit_logit_cr(
    an, v, COV_M4, "A1", lab, "(1) Model 4 reference: centre volume only, cluster-robust SE",
    "reference = Phase 2 table13; centre identity NOT in the mean structure")
  a1_rows[[length(a1_rows) + 1]] <- fit_glmer_sex(
    an, v, COV_M4_NOCENTER, "A1", lab, "(2) Within-centre: random intercept for centre (full Laplace)",
    "Model 4 covariates minus centre volume (absorbed by centre)", nagq = 1L)
  a1_rows[[length(a1_rows) + 1]] <- fit_glmer_sex(
    an, v, COV_M4_NOCENTER, "A1", lab, "(2b) Within-centre: random intercept, fast approximation (nAGQ=0)",
    "equivalence check for the fast approximation used in A3-A5", nagq = 0L)
  a1_rows[[length(a1_rows) + 1]] <- fit_fe_logit(
    an, v, COV_M4_NOCENTER, "A1", lab, "(3) Within-centre: centre fixed effects, cluster-robust SE",
    "Model 4 covariates minus centre volume (absorbed by centre)")
  logmsg("[A1] %s done (%.0f sec)", lab, as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
a1_tab <- dplyr::bind_rows(a1_rows)

# --- 주 예후 모형(DMS Model 4)의 센터 랜덤절편 버전 병기 --------------------
d_dms <- make_surv(an, "DMS")
m_cox_ref <- cox_model(d_dms, COV_M4)
s_ref <- summary(m_cox_ref)
cnt_dms <- list(n = m_cox_ref$n, n_events = m_cox_ref$nevent,
                n_male = sum(d_dms$female == 0), n_female = sum(d_dms$female == 1),
                events_male = sum(d_dms$.event[d_dms$female == 0] == 1),
                events_female = sum(d_dms$.event[d_dms$female == 1] == 1))
a1_cox <- std_row("A1", OUT$DMS$lab, "(1) Model 4 reference: Cox, cluster-robust SE", cnt_dms, "HR",
                  s_ref$conf.int["female", 1], s_ref$conf.int["female", 3],
                  s_ref$conf.int["female", 4], s_ref$coefficients["female", "Pr(>|z|)"],
                  "primary prognostic model")

m_fr <- survival::coxph(stats::as.formula(paste0(
  "Surv(.time, .event) ~ female + ", paste(COV_M4_NOCENTER, collapse = " + "),
  " + frailty(center_anon)")), data = d_dms, ties = "efron")
b_fr <- coef(m_fr)["female"]; se_fr <- sqrt(m_fr$var[1, 1])
theta <- tryCatch(m_fr$history[[1]]$theta, error = function(e) NA_real_)
a1_cox_fr <- std_row("A1", OUT$DMS$lab, "(2) Within-centre: Cox with shared gamma frailty by centre",
                     cnt_dms, "HR", exp(b_fr), exp(b_fr - 1.96 * se_fr), exp(b_fr + 1.96 * se_fr),
                     summary(m_fr)$coefficients["female", "p"],
                     sprintf("frailty variance = %.4f", theta))

m_str <- survival::coxph(stats::as.formula(paste0(
  "Surv(.time, .event) ~ female + ", paste(COV_M4_NOCENTER, collapse = " + "),
  " + strata(center_anon)")), data = d_dms, cluster = d_dms$center_anon,
  robust = TRUE, ties = "efron")
s_str <- summary(m_str)
a1_cox_st <- std_row("A1", OUT$DMS$lab, "(3) Within-centre: Cox stratified by centre (centre-specific baseline hazard)",
                     cnt_dms, "HR", s_str$conf.int["female", 1], s_str$conf.int["female", 3],
                     s_str$conf.int["female", 4], s_str$coefficients["female", "Pr(>|z|)"],
                     "exact within-centre comparison of the prognostic model")

a1_tab <- dplyr::bind_rows(a1_tab, a1_cox, a1_cox_fr, a1_cox_st)
write.csv(a1_tab, file.path(DIR_TAB, "table14_within_center_management.csv"),
          row.names = FALSE, na = "")

# --- 성별 × 센터 분포 진단 ---------------------------------------------------
# ⛔ PHI: 센터 식별자는 출력하지 않는다. 규모 순위로 재부여한 익명 인덱스만 사용.
ctr <- co %>%
  dplyr::group_by(center_anon) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_female = sum(female == 1),
    pct_female = 100 * mean(female == 1),
    centre_volume_category = names(sort(table(as.character(Center_volume_f)), decreasing = TRUE))[1],
    pct_IVUS = 100 * mean(IVUS, na.rm = TRUE),
    pct_bilateral = 100 * mean(Bilati_inj, na.rm = TRUE),
    pct_retrograde = 100 * mean(retro_attempt_any, na.rm = TRUE),
    pct_antegrade_only = 100 * mean(ante_approach_only, na.rm = TRUE),
    pct_potent_P2Y12 = 100 * mean(dapt_potent, na.rm = TRUE),
    pct_statin = 100 * mean(dc_statin, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n)) %>%
  dplyr::mutate(centre_index = dplyr::row_number()) %>%
  dplyr::select(centre_index, n, n_female, pct_female, centre_volume_category,
                dplyr::starts_with("pct_")) %>%
  dplyr::mutate(dplyr::across(dplyr::starts_with("pct_"), ~ round(.x, 1)))
ctr$note <- "centre_index = anonymised rank by centre size; NOT the clustering key used in models"
write.csv(ctr, file.path(DIR_TAB, "tableS12_center_distribution_by_sex.csv"),
          row.names = FALSE, na = "")

qs <- stats::quantile(ctr$pct_female, c(0, .25, .5, .75, 1), na.rm = TRUE)
cor_rows <- do.call(rbind, lapply(
  c("pct_IVUS", "pct_bilateral", "pct_retrograde", "pct_antegrade_only",
    "pct_potent_P2Y12", "pct_statin"), function(v) {
      ct2 <- suppressWarnings(stats::cor.test(ctr$pct_female, ctr[[v]], method = "spearman"))
      data.frame(item = paste0("Spearman rho: centre % female vs centre ", v),
                 value = sprintf("%.3f (p = %s)", unname(ct2$estimate), fmt_p(ct2$p.value)),
                 stringsAsFactors = FALSE)
    }))
icc_rows <- a1_tab %>%
  dplyr::filter(grepl("^\\(2\\) Within-centre: random intercept", model)) %>%
  dplyr::transmute(item = paste0("ICC (centre) for ", outcome),
                   value = sprintf("%.3f (random-intercept variance %.3f)", ICC, re_variance))
s12b <- rbind(
  data.frame(item = "Number of contributing centres", value = as.character(nrow(ctr))),
  data.frame(item = "Centre size: min / median / max",
             value = sprintf("%d / %d / %d", min(ctr$n), as.integer(stats::median(ctr$n)), max(ctr$n))),
  data.frame(item = "Centre-level % female: min / Q1 / median / Q3 / max",
             value = sprintf("%.1f / %.1f / %.1f / %.1f / %.1f", qs[1], qs[2], qs[3], qs[4], qs[5])),
  data.frame(item = "Centres with % female below 10% or above 30%",
             value = sprintf("%d of %d", sum(ctr$pct_female < 10 | ctr$pct_female > 30), nrow(ctr))),
  cor_rows, as.data.frame(icc_rows))
write.csv(s12b, file.path(DIR_TAB, "tableS12b_center_heterogeneity_summary.csv"),
          row.names = FALSE, na = "")

###############################################################################
# A2. 다중비교 정식 처리 + 복합 procedural intensity score
#
# 왜: 5개 양성 소견 중 3개(IVUS / potent P2Y12 / retrograde)가 검정 family 정의에
#     따라 생사가 갈린다. 또한 ante_approach_only 와 retro_attempt_any 는 사실상
#     같은 대비의 중복 계수다.
###############################################################################

message("== A2: multiplicity correction + procedural intensity score ==")

# --- (a) 검정 family 선언 및 보정 -------------------------------------------
# family = 본 분석에서 관리전략 지표로 정의한 6개 (00_cohort_prep.R::STRAT_VARS).
# 참고 family(m=14) = 위 6개 + 원내/시술 결과 8개 (table13 의 보정 비교 전체).
INHOS_VARS <- list(
  c("cto_success", "Technical success of CTO PCI"),
  c("inhos_any",   "Any in-hospital adverse event"),
  c("inhos_death", "In-hospital death"),
  c("PeriMI",      "Periprocedural myocardial infarction"),
  c("inhos_CVA",   "In-hospital stroke"),
  c("CIN_2",       "Contrast-induced nephropathy"),
  c("Final_TIMI3", "Final TIMI 3 flow"),
  c("CR_yes",      "Complete revascularization"))

m4_all <- dplyr::bind_rows(lapply(c(MGMT, INHOS_VARS), function(z)
  fit_logit_cr(an, z[1], COV_M4, "A2", z[2], "Model 4 adjusted logistic, cluster-robust SE")))
m4_all$family_m14 <- TRUE
m4_all$family_m6  <- m4_all$outcome %in% sapply(MGMT, `[`, 2)

mult <- m4_all %>% dplyr::filter(family_m6) %>%
  dplyr::select(outcome, n, n_events, pct_male, pct_female, estimate, lo, hi, est_fmt, p_value)
mult$p_BH_m6   <- stats::p.adjust(mult$p_value, method = "BH")
mult$p_holm_m6 <- stats::p.adjust(mult$p_value, method = "holm")
bh14 <- stats::p.adjust(m4_all$p_value, method = "BH")
mult$p_BH_m14 <- bh14[match(mult$outcome, m4_all$outcome)]

# within-centre(랜덤절편) 추정치에도 동일한 보정 적용
ri <- a1_tab %>% dplyr::filter(grepl("^\\(2\\) Within-centre: random intercept", model)) %>%
  dplyr::select(outcome, ri_OR = estimate, ri_lo = lo, ri_hi = hi, ri_p = p_value)
mult <- dplyr::left_join(mult, ri, by = "outcome")
mult$ri_est_fmt   <- sprintf("%.2f (%.2f-%.2f)", mult$ri_OR, mult$ri_lo, mult$ri_hi)
mult$ri_p_BH_m6   <- stats::p.adjust(mult$ri_p, method = "BH")
mult$ri_p_holm_m6 <- stats::p.adjust(mult$ri_p, method = "holm")

for (v in c("p_value", "p_BH_m6", "p_holm_m6", "p_BH_m14", "ri_p", "ri_p_BH_m6", "ri_p_holm_m6"))
  mult[[paste0(v, "_fmt")]] <- fmt_p(mult[[v]])
mult$survives_BH_m6    <- mult$p_BH_m6   < 0.05
mult$survives_holm_m6  <- mult$p_holm_m6 < 0.05
mult$survives_BH_m14   <- mult$p_BH_m14  < 0.05
mult$survives_BH_m6_within_centre   <- mult$ri_p_BH_m6   < 0.05
mult$survives_holm_m6_within_centre <- mult$ri_p_holm_m6 < 0.05

# 중복 계수(동일 대비 재표기) 정량화
ov <- table(an$ante_approach_only, an$retro_attempt_any, useNA = "no")
n_disc <- sum(an$ante_approach_only == 1 & an$retro_attempt_any == 1, na.rm = TRUE)
n_tot  <- sum(!is.na(an$ante_approach_only) & !is.na(an$retro_attempt_any))
mult$redundancy_note <- ifelse(
  mult$outcome %in% c("Antegrade approach only", "Retrograde approach attempted"),
  sprintf(paste0("Same contrast reported in opposite directions: antegrade-only equals the complement of ",
                 "retrograde-attempted in %d of %d patients (%.1f%%); they differ only for the %d patients ",
                 "in whom retrograde crossing was attempted but failed and the case was completed antegradely. ",
                 "These two rows are NOT independent findings."),
          n_tot - n_disc, n_tot, 100 * (n_tot - n_disc) / n_tot, n_disc), "")
mult$family_definition <- paste0(
  "Family of six management-strategy indicators defined for this analysis ",
  "(coded as STRAT_VARS in scripts/00_cohort_prep.R). BH and Holm reported side by side; ",
  "BH over the wider family of 14 adjusted comparisons in table13 also shown. ",
  "This analysis was not prospectively registered.")
write.csv(mult, file.path(DIR_TAB, "table13b_management_multiplicity.csv"),
          row.names = FALSE, na = "")

# --- (b) Procedural intensity score (0-4) ------------------------------------
# 구성: bilateral injection + retrograde attempted + IVUS + potent P2Y12
# 결측 처리: 주 분석 = 4개 성분 모두 관측된 환자 (complete case)
#            민감도  = 결측을 "시행하지 않음(0)" 으로 간주
IC <- c("Bilati_inj", "retro_attempt_any", "IVUS", "dapt_potent")
an$pis_cc <- ifelse(stats::complete.cases(an[, IC]), rowSums(an[, IC]), NA_real_)
an$pis_mz <- rowSums(sapply(IC, function(v) ifelse(is.na(an[[v]]), 0, an[[v]])))

pis_desc <- do.call(rbind, lapply(c("pis_cc", "pis_mz"), function(v) {
  d <- an[!is.na(an[[v]]), ]
  do.call(rbind, lapply(c("Male", "Female"), function(g) {
    x <- d[[v]][d$sex_f == g]
    data.frame(score_definition = ifelse(v == "pis_cc", "Complete case (all 4 components observed)",
                                         "Missing component treated as not performed"),
               sex = g, n = length(x),
               mean_sd = sprintf("%.2f (%.2f)", mean(x), stats::sd(x)),
               median_iqr = sprintf("%.0f (%.0f-%.0f)", stats::median(x),
                                    stats::quantile(x, .25), stats::quantile(x, .75)),
               pct_0 = round(100 * mean(x == 0), 1), pct_1 = round(100 * mean(x == 1), 1),
               pct_2 = round(100 * mean(x == 2), 1), pct_3 = round(100 * mean(x == 3), 1),
               pct_4 = round(100 * mean(x == 4), 1),
               stringsAsFactors = FALSE)
  }))
}))
write.csv(pis_desc, file.path(DIR_TAB, "table15b_procedural_intensity_distribution.csv"),
          row.names = FALSE, na = "")

pis_models <- function(scorevar, score_lab) {
  d <- .cc_subset(an, scorevar, COV_M4_NOCENTER)
  d$.s <- d[[scorevar]]
  cnt <- list(n = nrow(d), n_events = NA_integer_,
              n_male = sum(d$female == 0), n_female = sum(d$female == 1),
              events_male = NA_integer_, events_female = NA_integer_)
  rhs <- paste(c("female", COV_M4_NOCENTER), collapse = " + ")
  out <- list()

  # (1) 선형혼합모형 — 보정 평균 점수차 (센터 랜덤절편)
  m1 <- tryCatch(lme4::lmer(stats::as.formula(paste0(".s ~ ", rhs, " + (1 | center_anon)")),
                            data = d, REML = TRUE), error = function(e) NULL)
  if (!is.null(m1)) {
    s <- summary(m1); b <- s$coefficients["female", "Estimate"]; se <- s$coefficients["female", "Std. Error"]
    p <- 2 * stats::pnorm(-abs(b / se))
    r <- std_row("A2", score_lab, "Linear mixed model, centre random intercept", cnt,
                 "adjusted mean difference in score (female - male)",
                 b, b - 1.96 * se, b + 1.96 * se, p,
                 "single primary test for the management-strategy family")
    r$est_fmt <- sprintf("%.3f (%.3f to %.3f)", b, b - 1.96 * se, b + 1.96 * se)
    out[[length(out) + 1]] <- r
  }
  # (2) Poisson GLM (count), cluster-robust SE -> RR
  m2 <- stats::glm(stats::as.formula(paste0(".s ~ ", rhs)), data = d, family = stats::poisson())
  vc <- sandwich::vcovCL(m2, cluster = d$center_anon)
  ct <- lmtest::coeftest(m2, vcov. = vc)
  b <- ct["female", "Estimate"]; se <- ct["female", "Std. Error"]
  out[[length(out) + 1]] <- std_row("A2", score_lab, "Poisson regression, cluster-robust SE", cnt,
                                    "RR (score as count)", exp(b), exp(b - 1.96 * se),
                                    exp(b + 1.96 * se), ct["female", "Pr(>|z|)"], "")
  # (3) Poisson GLMM, 센터 랜덤절편
  m3 <- tryCatch(lme4::glmer(stats::as.formula(paste0(".s ~ ", rhs, " + (1 | center_anon)")),
                             data = d, family = stats::poisson(), nAGQ = 1L,
                             control = lme4::glmerControl(optimizer = "bobyqa",
                                                          optCtrl = list(maxfun = 3e5))),
                 error = function(e) NULL)
  if (!is.null(m3)) {
    s <- summary(m3); b <- s$coefficients["female", "Estimate"]; se <- s$coefficients["female", "Std. Error"]
    out[[length(out) + 1]] <- std_row("A2", score_lab, "Poisson GLMM, centre random intercept", cnt,
                                      "RR (score as count)", exp(b), exp(b - 1.96 * se),
                                      exp(b + 1.96 * se), s$coefficients["female", "Pr(>|z|)"],
                                      sprintf("centre random-intercept variance = %.4f",
                                              as.numeric(lme4::VarCorr(m3)$center_anon[1])))
  }
  # (4) 순서형 로지스틱 (비례오즈), cluster-robust SE
  d$.so <- factor(d$.s, ordered = TRUE)
  m4 <- tryCatch(MASS::polr(stats::as.formula(paste0(".so ~ ", rhs)), data = d, Hess = TRUE,
                            method = "logistic"), error = function(e) NULL)
  if (!is.null(m4)) {
    vcp <- tryCatch(sandwich::vcovCL(m4, cluster = d$center_anon), error = function(e) stats::vcov(m4))
    b <- coef(m4)["female"]; se <- sqrt(vcp["female", "female"])
    p <- 2 * stats::pnorm(-abs(b / se))
    out[[length(out) + 1]] <- std_row("A2", score_lab, "Ordinal logistic (proportional odds), cluster-robust SE",
                                      cnt, "OR (one-category-higher intensity)", exp(b),
                                      exp(b - 1.96 * se), exp(b + 1.96 * se), p, "")
  }
  dplyr::bind_rows(out)
}
pis_tab <- dplyr::bind_rows(
  pis_models("pis_cc", "Procedural intensity score 0-4 (complete case)"),
  pis_models("pis_mz", "Procedural intensity score 0-4 (missing = not performed)"))
pis_tab$note <- paste0(pis_tab$note,
  " | Score = bilateral injection + retrograde attempted + IVUS + potent P2Y12 (0-4). ",
  "Covariates: Model 4 minus centre volume (absorbed by centre).")
write.csv(pis_tab, file.path(DIR_TAB, "table15_procedural_intensity_score.csv"),
          row.names = FALSE, na = "")

###############################################################################
# A3. 체격 공통지지(common support) 제한 분석 ★핵심 방어
#
# 왜: 신장 SMD -2.37 로 positivity 위반 -> 모든 모델에서 배제됨. 그러나 여성의
#     stent 직경이 유의하게 작아(2.88 vs 2.98 mm) 혈관경 차이는 실재하며,
#     IVUS 통과성·양측 access·septal channel 크기의 임상적으로 정당한 결정요인이다.
###############################################################################

message("== A3: common-support (height) restricted re-analysis ==")

h <- co$Height
p5m  <- unname(stats::quantile(h[co$female == 0], .05, na.rm = TRUE))
p10m <- unname(stats::quantile(h[co$female == 0], .10, na.rm = TRUE))
p90f <- unname(stats::quantile(h[co$female == 1], .90, na.rm = TRUE))
p95f <- unname(stats::quantile(h[co$female == 1], .95, na.rm = TRUE))
CS_STRICT <- c(p5m, p95f)                      # 남성 5백분위수 ~ 여성 95백분위수
CS_WIDE   <- c(unname(stats::quantile(h[co$female == 0], .025, na.rm = TRUE)),
               unname(stats::quantile(h[co$female == 1], .975, na.rm = TRUE)))
logmsg("[A3] height common support: strict %.0f-%.0f cm, wide %.0f-%.0f cm (male P5=%.0f, P10=%.0f; female P90=%.0f, P95=%.0f)",
       CS_STRICT[1], CS_STRICT[2], CS_WIDE[1], CS_WIDE[2], p5m, p10m, p90f, p95f)

in_cs <- function(d, rng) !is.na(d$Height) & d$Height >= rng[1] & d$Height <= rng[2]

# --- balance 표: 제한 전/후 체격 변수 |SMD| ---------------------------------
smd_cont <- function(x, g) {
  x1 <- x[g == 1]; x0 <- x[g == 0]
  (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) /
    sqrt((stats::var(x1, na.rm = TRUE) + stats::var(x0, na.rm = TRUE)) / 2)
}
bal_rows <- list()
for (cfg in list(list(lab = "Full cohort (no height restriction)", d = an),
                 list(lab = sprintf("Common support, strict (%.0f-%.0f cm)", CS_STRICT[1], CS_STRICT[2]),
                      d = an[in_cs(an, CS_STRICT), ]),
                 list(lab = sprintf("Common support, wide (%.0f-%.0f cm)", CS_WIDE[1], CS_WIDE[2]),
                      d = an[in_cs(an, CS_WIDE), ]))) {
  dd <- cfg$d
  for (v in c("Height", "Weight", "bsa", "BMI_re", "Age")) {
    bal_rows[[length(bal_rows) + 1]] <- data.frame(
      restriction = cfg$lab, variable = v,
      n_male = sum(dd$female == 0 & !is.na(dd[[v]])),
      n_female = sum(dd$female == 1 & !is.na(dd[[v]])),
      mean_male = round(mean(dd[[v]][dd$female == 0], na.rm = TRUE), 2),
      mean_female = round(mean(dd[[v]][dd$female == 1], na.rm = TRUE), 2),
      SMD = round(smd_cont(dd[[v]], dd$female), 3),
      abs_SMD = round(abs(smd_cont(dd[[v]], dd$female)), 3),
      stringsAsFactors = FALSE)
  }
}
bal_cs <- dplyr::bind_rows(bal_rows)
bal_cs$balance_achieved <- ifelse(bal_cs$abs_SMD < 0.10, "yes (|SMD|<0.10)",
                                  ifelse(bal_cs$abs_SMD < 0.25, "acceptable (|SMD|<0.25)", "no"))
write.csv(bal_cs, file.path(DIR_TAB, "tableS13_balance_common_support.csv"),
          row.names = FALSE, na = "")

# --- 제한 코호트에서 6개 지표 재적합 -----------------------------------------
# 표본이 작아지므로 (i) 랜덤절편 + Model 4 공변량, (ii) 절약형 공변량 + cluster-robust
COV_PARSIM_PROC <- c("Age", "DM", "CKD_2", "JCTO_re", "Rentrop3", "CTO_location_pm",
                     "Vessel.Disease_f", "index_half_f")
a3_rows <- list()
for (cfg in list(list(tag = sprintf("Common support, strict (%.0f-%.0f cm)", CS_STRICT[1], CS_STRICT[2]),
                      d = an[in_cs(an, CS_STRICT), ]),
                 list(tag = sprintf("Common support, wide (%.0f-%.0f cm)", CS_WIDE[1], CS_WIDE[2]),
                      d = an[in_cs(an, CS_WIDE), ]))) {
  for (z in MGMT) {
    a3_rows[[length(a3_rows) + 1]] <- fit_glmer_sex(
      cfg$d, z[1], COV_M4_NOCENTER, "A3", z[2],
      paste0(cfg$tag, " | Model 4 covariates + centre random intercept"))
    a3_rows[[length(a3_rows) + 1]] <- fit_logit_cr(
      cfg$d, z[1], COV_PARSIM_PROC, "A3", z[2],
      paste0(cfg$tag, " | parsimonious covariates, cluster-robust SE"),
      "age, DM, CKD, J-CTO, Rentrop 3, prox/mid CTO, vessel disease, index half")
  }
}
# --- 보조: 전 코호트에서 신장을 공변량으로 추가 (외삽 경고 병기) -------------
for (z in MGMT) {
  a3_rows[[length(a3_rows) + 1]] <- fit_logit_cr(
    an, z[1], c(COV_M4, "Height"), "A3", z[2],
    "Full cohort with height as a covariate, cluster-robust SE",
    paste0("EXTRAPOLATION WARNING: height distributions overlap only marginally (|SMD| ",
           sprintf("%.2f", abs(smd_cont(an$Height, an$female))),
           "); this model relies on the parametric form outside the region of common support"))
}
a3_tab <- dplyr::bind_rows(a3_rows)
write.csv(a3_tab, file.path(DIR_TAB, "table16_common_support_height_management.csv"),
          row.names = FALSE, na = "")

###############################################################################
# A4. Potent P2Y12 격차의 분해 — 약제별 + 라벨 제한 하위군
#
# 왜: 격차가 prasugrel 주도(2.8% vs 7.1%)이고 prasugrel 은 >=75세 / <60 kg 라벨
#     제한이 있다. 여성 median 체중 57.0 kg, >=75세 43.6% -> 라벨 준수 처방과
#     구분되지 않는다.
###############################################################################

message("== A4: antiplatelet decomposition and label-unrestricted subgroup ==")

# 약제별 이진 결과 (DAPT 분모 내)
an$tica_vs_clop  <- dplyr::case_when(an$dc_antiplt_2 == "at" ~ 1, an$dc_antiplt_2 == "ac" ~ 0, TRUE ~ NA_real_)
an$pras_vs_clop  <- dplyr::case_when(an$dc_antiplt_2 == "ap" ~ 1, an$dc_antiplt_2 == "ac" ~ 0, TRUE ~ NA_real_)
# 분모 민감도: 퇴원 항혈소판 요법이 기록된 전체
an$potent_any <- ifelse(is.na(an$dc_antiplt_2) | (is.na(an$dc_prasugrel) & is.na(an$dc_ticagrelor)),
                        NA_real_,
                        as.numeric((!is.na(an$dc_prasugrel) & an$dc_prasugrel == 1) |
                                     (!is.na(an$dc_ticagrelor) & an$dc_ticagrelor == 1)))
# 라벨 비제한 하위군 (prasugrel 라벨: >=75세 또는 <60 kg 에서 제한/감량)
an$label_free <- !is.na(an$Age) & !is.na(an$Weight) & an$Age < 75 & an$Weight >= 60

a4 <- list()
add4 <- function(x) if (!is.null(x)) a4[[length(a4) + 1]] <<- x

# 1) 약제별 분리 모형 (DAPT 분모)
add4(fit_logit_cr(an, "dapt_potent", COV_M4, "A4", "Potent P2Y12 vs clopidogrel (DAPT denominator)",
                  "Model 4 reference, cluster-robust SE"))
add4(fit_glmer_sex(an, "dapt_potent", COV_M4_NOCENTER, "A4", "Potent P2Y12 vs clopidogrel (DAPT denominator)",
                   "Model 4 covariates + centre random intercept"))
add4(fit_logit_cr(an, "tica_vs_clop", COV_M4, "A4", "Ticagrelor vs clopidogrel (within DAPT)",
                  "Model 4 reference, cluster-robust SE"))
add4(fit_glmer_sex(an, "tica_vs_clop", COV_M4_NOCENTER, "A4", "Ticagrelor vs clopidogrel (within DAPT)",
                   "Model 4 covariates + centre random intercept"))
add4(fit_logit_cr(an, "pras_vs_clop", COV_M4, "A4", "Prasugrel vs clopidogrel (within DAPT)",
                  "Model 4 reference, cluster-robust SE"))
add4(fit_glmer_sex(an, "pras_vs_clop", COV_M4_NOCENTER, "A4", "Prasugrel vs clopidogrel (within DAPT)",
                   "Model 4 covariates + centre random intercept"))

# 2) 라벨 비제한 하위군 (Age < 75 & Weight >= 60)
an_lf <- an[an$label_free, , drop = FALSE]
logmsg("[A4] label-unrestricted subgroup (age<75 & weight>=60 kg): N = %d (M %d / F %d)",
       nrow(an_lf), sum(an_lf$female == 0), sum(an_lf$female == 1))
add4(fit_logit_cr(an_lf, "dapt_potent", COV_M4, "A4", "Potent P2Y12 (DAPT denominator)",
                  "Label-unrestricted subgroup: age <75 y AND weight >=60 kg, cluster-robust SE",
                  "prasugrel label restricts use at age >=75 y or weight <60 kg"))
add4(fit_glmer_sex(an_lf, "dapt_potent", COV_M4_NOCENTER, "A4", "Potent P2Y12 (DAPT denominator)",
                   "Label-unrestricted subgroup + centre random intercept"))
add4(fit_logit_cr(an_lf, "tica_vs_clop", COV_M4, "A4", "Ticagrelor vs clopidogrel (within DAPT)",
                  "Label-unrestricted subgroup: age <75 y AND weight >=60 kg, cluster-robust SE",
                  "ticagrelor has no weight-based label restriction"))
add4(fit_glmer_sex(an_lf, "tica_vs_clop", COV_M4_NOCENTER, "A4", "Ticagrelor vs clopidogrel (within DAPT)",
                   "Label-unrestricted subgroup + centre random intercept"))
add4(fit_logit_cr(an_lf, "pras_vs_clop", COV_M4, "A4", "Prasugrel vs clopidogrel (within DAPT)",
                  "Label-unrestricted subgroup: age <75 y AND weight >=60 kg, cluster-robust SE"))

# 3) 분모 민감도 — 퇴원 항혈소판 요법이 기록된 전체
add4(fit_logit_cr(an, "potent_any", COV_M4, "A4", "Any potent P2Y12 (all recorded discharge regimens)",
                  "Denominator sensitivity: all patients with a recorded discharge regimen, cluster-robust SE",
                  "avoids post-treatment conditioning on DAPT"))
add4(fit_glmer_sex(an, "potent_any", COV_M4_NOCENTER, "A4", "Any potent P2Y12 (all recorded discharge regimens)",
                   "Denominator sensitivity + centre random intercept"))
add4(fit_logit_cr(an_lf, "potent_any", COV_M4, "A4", "Any potent P2Y12 (all recorded discharge regimens)",
                  "Denominator sensitivity, label-unrestricted subgroup, cluster-robust SE"))

# 4) 체중을 공변량으로 추가 (외삽 경고 병기)
add4(fit_logit_cr(an, "dapt_potent", c(COV_M4, "Weight"), "A4", "Potent P2Y12 vs clopidogrel (DAPT denominator)",
                  "Model 4 + weight as covariate, cluster-robust SE",
                  sprintf("EXTRAPOLATION WARNING: weight |SMD| = %.2f between sexes; positivity is marginal",
                          abs(smd_cont(an$Weight, an$female)))))
add4(fit_logit_cr(an, "pras_vs_clop", c(COV_M4, "Weight"), "A4", "Prasugrel vs clopidogrel (within DAPT)",
                  "Model 4 + weight as covariate, cluster-robust SE",
                  "weight is the label criterion itself; this model removes the label pathway"))

a4_tab <- dplyr::bind_rows(a4)

# 기술통계 보강: 라벨 기준 충족 분포
lab_desc <- data.frame(
  item = c("Median weight, male / female (kg)",
           "% weight <60 kg, male / female",
           "% age >=75 y, male / female",
           "% meeting either prasugrel label restriction (age >=75 or weight <60), male / female",
           "Label-unrestricted subgroup n (male / female)"),
  value = c(sprintf("%.1f / %.1f", stats::median(an$Weight[an$female == 0], na.rm = TRUE),
                    stats::median(an$Weight[an$female == 1], na.rm = TRUE)),
            sprintf("%.1f%% / %.1f%%", 100 * mean(an$Weight[an$female == 0] < 60, na.rm = TRUE),
                    100 * mean(an$Weight[an$female == 1] < 60, na.rm = TRUE)),
            sprintf("%.1f%% / %.1f%%", 100 * mean(an$age75[an$female == 0], na.rm = TRUE),
                    100 * mean(an$age75[an$female == 1], na.rm = TRUE)),
            sprintf("%.1f%% / %.1f%%",
                    100 * mean((an$Age >= 75 | an$Weight < 60)[an$female == 0], na.rm = TRUE),
                    100 * mean((an$Age >= 75 | an$Weight < 60)[an$female == 1], na.rm = TRUE)),
            sprintf("%d / %d", sum(an_lf$female == 0), sum(an_lf$female == 1))),
  stringsAsFactors = FALSE)
write.csv(lab_desc, file.path(DIR_TAB, "table17b_prasugrel_label_criteria_by_sex.csv"),
          row.names = FALSE, na = "")
write.csv(a4_tab, file.path(DIR_TAB, "table17_antiplatelet_decomposition.csv"),
          row.names = FALSE, na = "")

###############################################################################
# A5. 에스컬레이션 "적응증" 조건부 분석
#
# 왜: 역행성·양측조영은 antegrade 실패 시의 bail-out 이다. 여성은 초기 와이어링
#     실패가 적었고(6.0% vs 9.0%) 병변도 덜 복잡했다. "에스컬레이션이 적었다"가
#     "필요가 적었다"의 결과라면 disparity 해석은 성립하지 않는다.
#
# ⚠️ Fail_tech / J-CTO 층화는 노출(성별) 이후 또는 시술 중 결정되는 변수에
#    조건을 거는 것이므로 collider 위험이 있다. 결과는 "적응증이 유사한 환자에서의
#    기술적(descriptive) 대비"로만 기술하며 인과 해석하지 않는다.
###############################################################################

message("== A5: escalation conditional on indication ==")

ESC <- list(c("retro_attempt_any", "Retrograde approach attempted"),
            c("Bilati_inj",        "Bilateral injection"),
            c("IVUS",              "Intravascular ultrasound use"))
COV_PARSIM_ESC <- c("Age", "JCTO_re", "Rentrop3", "CTO_Vessel_f", "index_half_f")

a5 <- list()
add5 <- function(x) if (!is.null(x)) a5[[length(a5) + 1]] <<- x

STRATA <- list(
  list(tag = "Indication stratum: failure of initial wiring technique (Fail_tech = 1)",
       d = an[an$Fail_tech %in% 1, , drop = FALSE]),
  list(tag = "Indication stratum: J-CTO score >=2",
       d = an[an$JCTO_re >= 2 & !is.na(an$JCTO_re), , drop = FALSE]),
  list(tag = "Indication stratum: occlusion length >=20 mm",
       d = an[an$Length20 %in% 1, , drop = FALSE]))
for (st in STRATA) {
  logmsg("[A5] %s: N = %d (M %d / F %d)", st$tag, nrow(st$d),
         sum(st$d$female == 0), sum(st$d$female == 1))
  for (z in ESC) {
    add5(fit_logit_cr(st$d, z[1], COV_PARSIM_ESC, "A5", z[2],
                      paste0(st$tag, " | parsimonious covariates, cluster-robust SE"),
                      "post-exposure stratification; descriptive contrast only, not a causal estimate"))
    add5(fit_logit_cr(st$d, z[1], COV_M4, "A5", z[2],
                      paste0(st$tag, " | Model 4 covariates, cluster-robust SE"),
                      "post-exposure stratification; may be unstable when the stratum is small"))
    add5(fit_glmer_sex(st$d, z[1], COV_M4_NOCENTER, "A5", z[2],
                       paste0(st$tag, " | Model 4 covariates + centre random intercept"),
                       "post-exposure stratification"))
  }
}

# 4) 역행성 성공률의 보정 분석 — 해부학 설명의 직접 검증
an_retro <- an[an$retro_attempt_any %in% 1, , drop = FALSE]
logmsg("[A5] retrograde attempted subgroup: N = %d (M %d / F %d); Retro_success observed in %d",
       nrow(an_retro), sum(an_retro$female == 0), sum(an_retro$female == 1),
       sum(!is.na(an_retro$Retro_success)))
add5(fit_logit_cr(an_retro, "Retro_success", character(0), "A5",
                  "Retrograde success (among patients in whom it was attempted)",
                  "Unadjusted, cluster-robust SE"))
add5(fit_logit_cr(an_retro, "Retro_success", COV_PARSIM_ESC, "A5",
                  "Retrograde success (among patients in whom it was attempted)",
                  "Parsimonious covariates (age, J-CTO, Rentrop 3, CTO vessel, index half), cluster-robust SE",
                  "direct test of the anatomical explanation for less frequent retrograde use in women"))
a5_tab <- dplyr::bind_rows(a5)
write.csv(a5_tab, file.path(DIR_TAB, "table18_escalation_conditional_on_indication.csv"),
          row.names = FALSE, na = "")

###############################################################################
# A6. 후기(>3년) 신호의 정직한 제시 — 구간별 노출량 + 3년 landmark
###############################################################################

message("== A6: interval exposure counts and 3-year landmark analysis ==")

TV_KEYS <- c("DMS", "Death", "NCVdeath", "MACCE_pub", "MACE_pub", "NACE_pub")
CUTS <- c(0, 365, 1095, Inf)
PER_LAB <- c("0-1 year", "1-3 years", ">3 years")

# --- 구간별 event 수 · 구간 시작 at-risk · person-time -----------------------
interval_exposure <- function(dat, key) {
  d <- make_surv(dat, key); d <- d[d$.time > 0, , drop = FALSE]
  do.call(rbind, lapply(seq_along(PER_LAB), function(i) {
    lo <- CUTS[i]; hi <- CUTS[i + 1]
    do.call(rbind, lapply(c("Male", "Female"), function(g) {
      ds <- d[d$sex_f == g, , drop = FALSE]
      at_risk <- sum(ds$.time > lo)
      ev <- sum(ds$.event == 1 & ds$.time > lo & ds$.time <= hi)
      pt <- sum(pmax(pmin(ds$.time, hi) - lo, 0)) / 365.25
      data.frame(outcome = OUT[[key]]$lab, outcome_key = key, period = PER_LAB[i], sex = g,
                 at_risk_at_interval_start = at_risk, events_in_interval = ev,
                 person_years_in_interval = round(pt, 1),
                 rate_per_100py = round(100 * ev / pt, 2),
                 stringsAsFactors = FALSE)
    }))
  }))
}
iv_tab <- dplyr::bind_rows(lapply(TV_KEYS, function(k) interval_exposure(an, k)))
iv_wide <- iv_tab %>%
  tidyr::pivot_wider(id_cols = c(outcome, outcome_key, period), names_from = sex,
                     values_from = c(at_risk_at_interval_start, events_in_interval,
                                     person_years_in_interval, rate_per_100py))

# 기존 table8b (구간별 보정 HR) 와 병합 + 18개 구간검정 FDR
f8b <- file.path(DIR_TAB, "table8b_time_varying_female_HR.csv")
if (file.exists(f8b)) {
  tv <- utils::read.csv(f8b, stringsAsFactors = FALSE)
  tv$period <- trimws(tv$period)
  tv$p_BH_18 <- stats::p.adjust(tv$p_value, method = "BH")
  tv$p_BH_18_fmt <- fmt_p(tv$p_BH_18)
  tv$multiplicity_note <- sprintf(
    "%d time-interval tests (6 outcomes x 3 intervals); Benjamini-Hochberg q reported. These interval analyses were not prospectively registered.",
    nrow(tv))
  iv_wide <- dplyr::left_join(tv, iv_wide, by = c("outcome", "outcome_key", "period"))
}
write.csv(iv_wide, file.path(DIR_TAB, "table8b2_time_varying_with_exposure.csv"),
          row.names = FALSE, na = "")

# --- 3년 landmark 분석 --------------------------------------------------------
LM <- 1095
landmark_fit <- function(key) {
  d <- make_surv(an, key)
  d <- d[d$.time > LM, , drop = FALSE]           # 3년 시점 event-free 생존자
  d$.time <- d$.time - LM
  r1 <- fit_cox(d, character(0), "Landmark at 3 years: unadjusted", OUT[[key]]$lab,
                "restricted to patients event-free and still under observation at 3 years")
  r2 <- fit_cox(d, COV_M4, "Landmark at 3 years: Model 4 (re-fitted from the 3-year landmark)",
                OUT[[key]]$lab, "depletion-of-susceptibles check for the >3-year signal")
  rbind(r1, r2)
}
lm_tab <- dplyr::bind_rows(lapply(c("DMS", "Death", "MACCE_pub", "MACE_pub", "NACE_pub"),
                                  function(k) tryCatch(landmark_fit(k), error = function(e) NULL)))
lm_tab$caveat <- "Landmark analyses condition on 3-year survival and are subject to depletion of susceptibles; they were not prospectively registered and are hypothesis-generating."
write.csv(lm_tab, file.path(DIR_TAB, "table19_landmark_3yr.csv"), row.names = FALSE, na = "")

###############################################################################
# A7. 중복 레코드 제거 민감도
#
# 왜: CLAUDE.md §2-1 은 완전동일행 141건 "제거", 재시술 103건 "index 채택"을 방침으로
#     명시했으나 00_cohort_prep.R 구현에는 어떤 dedup 도 없다.
# ⛔ PHI: patient_key / id_code / id 의 값은 계산에만 쓰고 출력하지 않는다.
###############################################################################

message("== A7: duplicate-record sensitivity ==")

dup_flags <- local({
  raw <- readRDS(file.path(PROJ, "data", "processed", "kcto_sex_analysis.rds"))
  stopifnot(nrow(raw) == 9316)
  raw$index_yr <- suppressWarnings(as.numeric(as.character(raw$index_yr)))
  sel <- raw$flag_cto_pci %in% TRUE & !is.na(raw$index_yr) &
    raw$index_yr >= 2014 & raw$index_yr <= 2023
  stopifnot(sum(sel) == nrow(co))
  sub <- raw[sel, , drop = FALSE]
  key <- paste(sub$id_code, sub$id, sep = "|")   # PHI — 계산 전용, 출력 금지

  # (a) 기록된 모든 변수가 동일한 행.
  #     idx / row_id 는 cleaning 단계에서 부여된 행 고유 일련번호이므로 제외한다
  #     (포함하면 어떤 두 행도 같아질 수 없어 검정 자체가 무의미해진다).
  dup_allcol <- duplicated(sub[, setdiff(names(sub), c("idx", "row_id")), drop = FALSE])

  # (b) 중복 입력 (CLAUDE.md §2-1 "완전동일행" 정의):
  #     동일 식별자 키 내에서 (시술일, 나이, 성별)이 모두 같은 행 -> 같은 시술의 중복 입력
  dup_entry <- duplicated(paste(key, as.character(sub$CTO_intervention_date),
                                sub$Age, sub$Sex, sep = "||"))

  # (c) 진짜 재시술: (b) 제거 후 남은 행 중 같은 키 + 같은 성별 + 나이궤적이
  #     등록연도 차이와 정합 -> 같은 환자의 재시술. index(최초) 시술만 유지.
  #     성별 불일치 / 나이궤적 비타당 그룹은 id 재사용이므로 병합하지 않는다.
  drop_repeat <- rep(FALSE, nrow(sub)); n_distinct_excess <- 0L
  keep1 <- which(!dup_entry)
  for (g in split(keep1, key[keep1])) {
    if (length(g) < 2) next
    sx <- sub$Sex[g]; ag <- suppressWarnings(as.numeric(sub$Age[g])); yr <- sub$index_yr[g]
    same_sex <- length(unique(sx[!is.na(sx)])) <= 1
    plausible <- !any(is.na(ag)) && !any(is.na(yr)) &&
      abs((max(ag) - min(ag)) - (max(yr) - min(yr))) <= 1
    if (same_sex && plausible) {
      drop_repeat[g[order(sub$index_yr[g])][-1]] <- TRUE
    } else n_distinct_excess <- n_distinct_excess + length(g) - 1L
  }
  list(dup_allcol = dup_allcol, dup_entry = dup_entry, drop_repeat = drop_repeat,
       n_distinct_excess = n_distinct_excess)
})
logmsg(paste0("[A7] analysis cohort: rows identical on every recorded variable = %d; ",
              "duplicate data-entry rows (same key + date/age/sex) = %d; ",
              "repeat-procedure rows = %d; excess rows retained as distinct patients = %d"),
       sum(dup_flags$dup_allcol), sum(dup_flags$dup_entry),
       sum(dup_flags$drop_repeat), dup_flags$n_distinct_excess)

CO_I  <- co[!dup_flags$dup_entry, , drop = FALSE]
CO_II <- co[!(dup_flags$dup_entry | dup_flags$drop_repeat), , drop = FALSE]

a7 <- list()
for (cfg in list(list(tag = sprintf("Primary cohort (no de-duplication), N = %d", nrow(co)), d = co),
                 list(tag = sprintf("(i) Duplicate data-entry rows removed, N = %d", nrow(CO_I)), d = CO_I),
                 list(tag = sprintf("(ii) Duplicate entries removed + index procedure retained for repeat procedures, N = %d",
                                    nrow(CO_II)), d = CO_II))) {
  dd <- cfg$d[stats::complete.cases(cfg$d[, COV_M4, drop = FALSE]), , drop = FALSE]
  r <- fit_cox(make_surv(dd, "DMS"), COV_M4, cfg$tag, OUT$DMS$lab, "primary outcome, Model 4")
  a7[[length(a7) + 1]] <- data.frame(analysis = "A7", cohort = cfg$tag, outcome = r$outcome,
                                     measure = "HR", n = r$n, n_events = r$n_events,
                                     estimate = r$HR, lo = r$lo, hi = r$hi, est_fmt = r$HR_fmt,
                                     p_value = r$p_value, p_fmt = r$p_fmt, stringsAsFactors = FALSE)
  for (z in MGMT) {
    rr <- fit_logit_cr(dd, z[1], COV_M4, "A7", z[2], cfg$tag)
    if (is.null(rr)) next
    a7[[length(a7) + 1]] <- data.frame(analysis = "A7", cohort = cfg$tag, outcome = rr$outcome,
                                       measure = "OR", n = rr$n, n_events = rr$n_events,
                                       estimate = rr$estimate, lo = rr$lo, hi = rr$hi,
                                       est_fmt = rr$est_fmt, p_value = rr$p_value,
                                       p_fmt = rr$p_fmt, stringsAsFactors = FALSE)
  }
}
a7_tab <- dplyr::bind_rows(a7)
a7_tab$strobe_note <- sprintf(
  paste0("STROBE flow footnote: within the analysis cohort of %d, no two records were identical on every recorded variable (%d rows). ",
         "%d rows shared an identifier with another record and had the same procedure date, age and sex, and were treated as duplicate data entry of the same procedure. ",
         "A further %d rows represented a repeat procedure in a patient already contributing an index procedure (same identifier, same sex, age trajectory consistent with the interval between index years); ",
         "for these, only the index procedure was retained. %d further excess rows shared an identifier but had an inconsistent sex or age trajectory and were retained as distinct patients, ",
         "because the registry identifier is re-used across patients (CLAUDE.md 2-1)."),
  nrow(co), sum(dup_flags$dup_allcol), sum(dup_flags$dup_entry),
  sum(dup_flags$drop_repeat), dup_flags$n_distinct_excess)
write.csv(a7_tab, file.path(DIR_TAB, "tableS14_duplicate_record_sensitivity.csv"),
          row.names = FALSE, na = "")

###############################################################################
# A8. 추적 0일 환자 관련 보완 (새 모형 없음 — 전부 기술통계)
#
# 왜: 추적 0일 환자는 person-time 을 기여하지 않아 이미 모든 time-to-event 모형에서
#     사실상 배제되어 있다. 보고 n=6,259 는 정보 기여 표본이 아니다.
###############################################################################

message("== A8: person-time contributing sample and zero-follow-up baseline ==")

pt_rows <- do.call(rbind, lapply(names(OUT), function(k) {
  d <- make_surv(an, k)
  dpos <- d[d$.time > 0, , drop = FALSE]
  data.frame(outcome = OUT[[k]]$lab, outcome_key = k,
             n_model_reported = nrow(d),
             n_contributing_person_time = nrow(dpos),
             n_zero_follow_up_excluded = nrow(d) - nrow(dpos),
             n_contributing_male = sum(dpos$female == 0),
             n_contributing_female = sum(dpos$female == 1),
             events = sum(d$.event == 1),
             person_years = round(sum(d$.time) / 365.25, 1),
             stringsAsFactors = FALSE)
}))
pt_rows$note <- paste0("Rows with zero follow-up and no event contribute to no risk set and are ",
                       "therefore uninformative in any time-to-event model; the reported model n ",
                       "overstates the information-contributing sample.")
write.csv(pt_rows, file.path(DIR_TAB, "tableS15b_person_time_contributing_n.csv"),
          row.names = FALSE, na = "")

# --- 추적 0일군 vs 나머지 baseline 비교 (전체 + 성별 층화) -------------------
Z_VARS <- list(
  c("Age", "cont", "Age, years"), c("age75", "bin", "Age >=75 years"),
  c("female", "bin", "Female sex"), c("BMI_re", "cont", "Body mass index"),
  c("HTN", "bin", "Hypertension"), c("DM", "bin", "Diabetes"),
  c("CKD_2", "bin", "CKD (eGFR<60)"), c("AF", "bin", "Atrial fibrillation"),
  c("preMI", "bin", "Prior myocardial infarction"), c("prePCI", "bin", "Prior PCI"),
  c("preCABG", "bin", "Prior CABG"), c("stable_angina", "bin", "Stable angina at presentation"),
  c("acute_MI", "bin", "Myocardial infarction at presentation"),
  c("JCTO_re", "cont", "J-CTO score"), c("JCTO_bi", "bin", "J-CTO score >=2"),
  c("Rentrop3", "bin", "Rentrop grade 3"),
  c("cto_success", "bin", "Technical success"),
  c("inhos_any", "bin", "Any in-hospital adverse event"),
  c("Center_volume", "cont", "Centre volume category (1-4)"))

# 추적 0일군의 정의는 tableS7 / CLAUDE.md §1-2 와 동일하게 맞춘다:
#   death_fudur == 0 & Death == 0  (person-time 을 전혀 기여하지 않는 환자)
# ⚠️ recode_common() 의 fu_zero(= fudur2 == 0)는 약간 다른 집합이며 원고에서
#    두 정의를 섞지 말 것. 두 수치를 아래 각주에 함께 명시한다.
co$zero_fu <- as.numeric(co$death_fudur == 0 & co$Death == 0)

zf_compare <- function(dat, tag) {
  g <- dat$zero_fu
  do.call(rbind, lapply(Z_VARS, function(z) {
    v <- z[1]; ty <- z[2]; lab <- z[3]
    x <- suppressWarnings(as.numeric(dat[[v]]))
    x1 <- x[g == 1]; x0 <- x[g == 0]
    if (ty == "cont") {
      val1 <- sprintf("%.1f (%.1f)", mean(x1, na.rm = TRUE), stats::sd(x1, na.rm = TRUE))
      val0 <- sprintf("%.1f (%.1f)", mean(x0, na.rm = TRUE), stats::sd(x0, na.rm = TRUE))
      p <- tryCatch(stats::t.test(x1, x0)$p.value, error = function(e) NA_real_)
    } else {
      val1 <- sprintf("%d (%.1f%%)", sum(x1 == 1, na.rm = TRUE), 100 * mean(x1, na.rm = TRUE))
      val0 <- sprintf("%d (%.1f%%)", sum(x0 == 1, na.rm = TRUE), 100 * mean(x0, na.rm = TRUE))
      p <- tryCatch(suppressWarnings(stats::chisq.test(table(x, g))$p.value), error = function(e) NA_real_)
    }
    data.frame(stratum = tag, variable = lab,
               n_zero_follow_up = sum(g == 1), n_with_follow_up = sum(g == 0),
               zero_follow_up = val1, with_follow_up = val0,
               SMD = round(smd_cont(x, g), 3), abs_SMD = round(abs(smd_cont(x, g)), 3),
               p_value = p, p_fmt = fmt_p(p), stringsAsFactors = FALSE)
  }))
}
zf_tab <- dplyr::bind_rows(
  zf_compare(co, "All patients"),
  zf_compare(co[co$sex_f == "Male", ], "Men"),
  zf_compare(co[co$sex_f == "Female", ], "Women"))
zf_tab$flow_note <- sprintf(
  paste0("STROBE flow box: no follow-up beyond the index procedure (death_fudur = 0 and no event): ",
         "%d (men %d, %.1f%%; women %d, %.1f%%; chi-square P = %s). ",
         "Under the alternative definition fudur2 = 0 the corresponding count is %d ",
         "(men %d, %.1f%%; women %d, %.1f%%) - report one definition only and state it explicitly."),
  sum(co$zero_fu == 1), sum(co$zero_fu == 1 & co$female == 0),
  100 * mean(co$zero_fu[co$female == 0] == 1), sum(co$zero_fu == 1 & co$female == 1),
  100 * mean(co$zero_fu[co$female == 1] == 1),
  fmt_p(suppressWarnings(stats::chisq.test(table(co$zero_fu, co$female))$p.value)),
  sum(co$fu_zero == 1, na.rm = TRUE), sum(co$fu_zero == 1 & co$female == 0, na.rm = TRUE),
  100 * mean(co$fu_zero[co$female == 0] == 1, na.rm = TRUE),
  sum(co$fu_zero == 1 & co$female == 1, na.rm = TRUE),
  100 * mean(co$fu_zero[co$female == 1] == 1, na.rm = TRUE))
write.csv(zf_tab, file.path(DIR_TAB, "tableS15_zero_followup_baseline.csv"),
          row.names = FALSE, na = "")

###############################################################################
# A9. E-value 벤치마크 재계산 (범주형 연령) + endpoint 별 연령 감쇠율
#
# 왜: 연속형 연령 벤치마크(1.30)는 E-value 1.95 와 자기모순. 범주형(>=75세 vs <65세)
#     으로 바꾸면 bounding factor 1.69 > 1.48 로 정합한다.
###############################################################################

message("== A9: categorical E-value benchmark and endpoint-specific age attenuation ==")

f6f <- file.path(DIR_TAB, "table6f_age_stratified_DMS.csv")
f7  <- file.path(DIR_TAB, "table7_adjusted_HR_all_outcomes.csv")
stopifnot(file.exists(f6f), file.exists(f7))
t6f <- utils::read.csv(f6f, stringsAsFactors = FALSE)
t7  <- utils::read.csv(f7,  stringsAsFactors = FALSE)

# RR(sex -> age>=75): 여성에서 >=75세 비율 / 남성에서 >=75세 비율
p75_f <- 100 * mean(co$age75[co$female == 1], na.rm = TRUE)
p75_m <- 100 * mean(co$age75[co$female == 0], na.rm = TRUE)
RR_EU <- p75_f / p75_m
# RR(age>=75 vs <65 -> 5년 DMS), 남성 내에서 (table6f)
inc75 <- t6f$cum_inc_5yr_pct[t6f$age_stratum == ">=75 years" & t6f$sex == "Male"]
inc65 <- t6f$cum_inc_5yr_pct[t6f$age_stratum == "<65 years"  & t6f$sex == "Male"]
RR_UD <- inc75 / inc65
BF <- (RR_EU * RR_UD) / (RR_EU + RR_UD - 1)

hr_crude <- t7$HR[t7$outcome_key == "DMS" & t7$model == "Unadjusted"]
lo_crude <- t7$lo[t7$outcome_key == "DMS" & t7$model == "Unadjusted"]
ev_crude <- EValue::evalues.HR(est = hr_crude, lo = lo_crude,
                               hi = t7$hi[t7$outcome_key == "DMS" & t7$model == "Unadjusted"],
                               rare = FALSE)

ev_bench <- data.frame(
  item = c("RR (female sex -> age >=75 years)",
           "RR (age >=75 vs <65 years -> 5-year DMS, men only)",
           "Bounding factor implied by categorical age",
           "Crude (unadjusted) female-vs-male HR for DMS",
           "E-value for the crude HR (point estimate)",
           "E-value for the crude HR (lower confidence limit)",
           "Does the measured age contrast exceed the E-value on both arms?",
           "Does the bounding factor exceed the crude HR?"),
  value = c(sprintf("%.2f (%.1f%% vs %.1f%%)", RR_EU, p75_f, p75_m),
            sprintf("%.2f (%.1f%% vs %.1f%%)", RR_UD, inc75, inc65),
            sprintf("%.2f", BF),
            sprintf("%.2f", hr_crude),
            sprintf("%.2f", ev_crude["E-values", "point"]),
            sprintf("%.2f", ev_crude["E-values", "lower"]),
            ifelse(RR_EU >= ev_crude["E-values", "point"] & RR_UD >= ev_crude["E-values", "point"],
                   "yes", "no"),
            ifelse(BF >= hr_crude, "yes", "no")),
  stringsAsFactors = FALSE)
ev_bench$interpretation <- c(rep("", 6),
  sprintf(paste0("Both measured age associations (%.2f and %.2f) exceed the E-value of %.2f required ",
                 "of an unmeasured confounder."), RR_EU, RR_UD, ev_crude["E-values", "point"]),
  sprintf(paste0("A confounder as strongly related to sex and to the outcome as categorical age ",
                 "could account for an association up to HR %.2f, which exceeds the observed crude ",
                 "HR of %.2f. The attenuation of the crude association after adjustment for age is ",
                 "therefore quantitatively consistent with the strength of the measured age contrast. ",
                 "Note the earlier continuous-age benchmark (adjusted HR 1.34 per 10 years, 9-year ",
                 "median difference -> 1.30) understates this and should not be used."), BF, hr_crude))
ev_bench$supersedes <- "Replaces the continuous-age benchmark in table11b_age_association_strength.csv"
write.csv(ev_bench, file.path(DIR_TAB, "table11c_evalue_categorical_benchmark.csv"),
          row.names = FALSE, na = "")

# --- endpoint 별 연령 감쇠율 ---------------------------------------------------
att <- do.call(rbind, lapply(unique(t7$outcome_key), function(k) {
  g <- function(mod, col) t7[[col]][t7$outcome_key == k & t7$model == mod]
  h1 <- g("Unadjusted", "HR"); h2 <- g("Age-adjusted", "HR"); h4 <- g("Fully adjusted", "HR")
  data.frame(outcome = t7$outcome[t7$outcome_key == k][1], outcome_key = k,
             n_events = t7$n_events[t7$outcome_key == k][1],
             crude_HR = sprintf("%s, p %s", g("Unadjusted", "HR_fmt"), g("Unadjusted", "p_fmt")),
             age_adjusted_HR = sprintf("%s, p %s", g("Age-adjusted", "HR_fmt"), g("Age-adjusted", "p_fmt")),
             fully_adjusted_HR = sprintf("%s, p %s", g("Fully adjusted", "HR_fmt"), g("Fully adjusted", "p_fmt")),
             pct_excess_removed_by_age = ifelse(h1 > 1, sprintf("%.1f%%", 100 * (h1 - h2) / (h1 - 1)),
                                                "not applicable (crude HR <=1)"),
             pct_excess_removed_by_full = ifelse(h1 > 1, sprintf("%.1f%%", 100 * (h1 - h4) / (h1 - 1)),
                                                 "not applicable (crude HR <=1)"),
             stringsAsFactors = FALSE)
}))
att$note <- paste0("Attenuation is descriptive of the change in the estimate after adjustment for ",
                   "differences in the age distribution; it is not a mediation analysis and does not ",
                   "support a causal decomposition.")
write.csv(att, file.path(DIR_TAB, "table20_age_attenuation_by_endpoint.csv"),
          row.names = FALSE, na = "")

###############################################################################
# A10. 희소 endpoint 의 경쟁위험 처리 (MI / stroke / TIMI major bleeding)
###############################################################################

message("== A10: competing-risk analysis for sparse endpoints ==")

cr_setup <- function(dat, ev, tm) {
  d <- dat
  e <- suppressWarnings(as.numeric(d[[ev]]))
  t <- suppressWarnings(as.numeric(d[[tm]]))
  dth <- suppressWarnings(as.numeric(d$Death))
  dtm <- suppressWarnings(as.numeric(d$death_fudur))
  status <- ifelse(!is.na(e) & e == 1, 1, ifelse(!is.na(dth) & dth == 1, 2, 0))
  time <- ifelse(status == 1, t, ifelse(status == 2, dtm, t))
  time[status > 0 & !is.na(time) & time == 0] <- 0.5
  ok <- !is.na(time) & !is.na(status) & !is.na(d$female)
  list(d = d[ok, , drop = FALSE], time = time[ok], status = status[ok])
}

a10 <- list()
for (k in c("MI", "Stroke", "MajBld")) {
  o <- OUT[[k]]
  cs <- cr_setup(an, o$ev, o$tm)
  d <- cs$d; tt <- cs$time; st <- cs$status
  # Aalen-Johansen CIF + Gray 검정
  ci <- cmprsk::cuminc(ftime = tt, fstatus = st, group = as.character(d$sex_f), cencode = 0)
  est <- cmprsk::timepoints(ci, times = c(365, 1095, 1825))
  gp <- ci$Tests[, "pv"]
  for (rn in rownames(est$est)) {
    if (!grepl(" 1$", rn)) next          # 관심 사건만 (경쟁사망 CIF 는 별도 행)
    for (j in 1:3) {
      a10[[length(a10) + 1]] <- data.frame(
        endpoint = o$lab, method = "Aalen-Johansen cumulative incidence (death as competing risk)",
        sex = ifelse(grepl("^Male", rn), "Male", "Female"),
        timepoint = c("1-year", "3-year", "5-year")[j],
        n = NA_integer_, n_events = NA_integer_, n_competing_deaths = NA_integer_,
        estimate = round(100 * est$est[rn, j], 2), lo = NA_real_, hi = NA_real_,
        est_fmt = sprintf("%.2f%% (SE %.2f)", 100 * est$est[rn, j], 100 * sqrt(est$var[rn, j])),
        p_value = unname(gp["1"]), p_fmt = fmt_p(unname(gp["1"])),
        note = "Gray test for equality of cumulative incidence functions",
        stringsAsFactors = FALSE)
    }
  }
  # Fine-Gray sHR
  fg_one <- function(covs, lab) {
    X <- stats::model.matrix(stats::as.formula(
      paste0("~ female", if (length(covs)) paste0(" + ", paste(covs, collapse = " + ")) else "")),
      data = d)[, -1, drop = FALSE]
    # model.matrix 는 결측행을 버리므로 원 데이터 내 **위치**로 되돌려 매핑한다
    keep <- match(rownames(X), rownames(d))
    stopifnot(!anyNA(keep))
    fg <- cmprsk::crr(ftime = tt[keep], fstatus = st[keep], cov1 = X, failcode = 1, cencode = 0)
    s <- summary(fg); i <- which(rownames(s$coef) == "female")
    data.frame(endpoint = o$lab, method = lab, sex = "Female vs male",
               timepoint = "over full follow-up",
               n = length(keep), n_events = sum(st[keep] == 1),
               n_competing_deaths = sum(st[keep] == 2),
               estimate = unname(s$conf.int[i, 1]), lo = unname(s$conf.int[i, 3]),
               hi = unname(s$conf.int[i, 4]),
               est_fmt = sprintf("%.2f (%.2f-%.2f)", s$conf.int[i, 1], s$conf.int[i, 3], s$conf.int[i, 4]),
               p_value = unname(s$coef[i, 5]), p_fmt = fmt_p(unname(s$coef[i, 5])),
               note = "Fine-Gray subdistribution hazard ratio", stringsAsFactors = FALSE)
  }
  a10[[length(a10) + 1]] <- tryCatch(fg_one(character(0), "Fine-Gray sHR (unadjusted)"),
                                     error = function(e) NULL)
  a10[[length(a10) + 1]] <- tryCatch(
    fg_one(COV_PARS, "Fine-Gray sHR (parsimoniously adjusted: age, DM, CKD, AF)"),
    error = function(e) NULL)
  # 대조: cause-specific Cox (동일 표본)
  dcs <- d; dcs$.time <- tt; dcs$.event <- as.integer(st == 1)
  r <- fit_cox(dcs, COV_PARS, "Cause-specific Cox (parsimoniously adjusted)", o$lab)
  a10[[length(a10) + 1]] <- data.frame(
    endpoint = o$lab, method = "Cause-specific hazard ratio (parsimoniously adjusted, for contrast)",
    sex = "Female vs male", timepoint = "over full follow-up",
    n = r$n, n_events = r$n_events, n_competing_deaths = sum(st == 2),
    estimate = r$HR, lo = r$lo, hi = r$hi, est_fmt = r$HR_fmt,
    p_value = r$p_value, p_fmt = r$p_fmt,
    note = "cause-specific vs subdistribution: divergence indicates an influence of the competing risk of death",
    stringsAsFactors = FALSE)
}
a10_tab <- dplyr::bind_rows(a10)
write.csv(a10_tab, file.path(DIR_TAB, "table9c_competing_risk_sparse_endpoints.csv"),
          row.names = FALSE, na = "")

###############################################################################
# 마무리 — 진단 로그 · sessionInfo
###############################################################################

capture.output({
  cat("K-CTO sex difference analysis — Phase 2.5 QA additional analyses (A1-A10)\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Cohort:", COHORT_LABEL, "\n")
  cat("Model 4 complete-case analysis set: N =", nrow(an), "\n\n")
  cat("--- run log ---\n"); for (s in LOG) cat(s, "\n")
  cat("\n--- A1: within-centre estimates (management strategy) ---\n")
  print(a1_tab[, c("outcome", "model", "n", "measure", "est_fmt", "p_fmt")])
  cat("\n--- A1: centre heterogeneity summary ---\n"); print(s12b)
  cat("\n--- A2: multiplicity ---\n")
  print(mult[, c("outcome", "est_fmt", "p_value_fmt", "p_BH_m6_fmt", "p_holm_m6_fmt",
                 "p_BH_m14_fmt", "ri_est_fmt", "ri_p_fmt", "ri_p_BH_m6_fmt", "ri_p_holm_m6_fmt")])
  cat("\n--- A2: procedural intensity score ---\n")
  print(pis_desc); print(pis_tab[, c("outcome", "model", "n", "measure", "est_fmt", "p_fmt")])
  cat("\n--- A3: common support balance ---\n"); print(bal_cs)
  cat("\n--- A3: management indicators within common support ---\n")
  print(a3_tab[, c("outcome", "model", "n", "pct_male", "pct_female", "est_fmt", "p_fmt")])
  cat("\n--- A4: antiplatelet decomposition ---\n")
  print(a4_tab[, c("outcome", "model", "n", "pct_male", "pct_female", "est_fmt", "p_fmt")])
  print(lab_desc)
  cat("\n--- A5: escalation conditional on indication ---\n")
  print(a5_tab[, c("outcome", "model", "n", "pct_male", "pct_female", "est_fmt", "p_fmt")])
  cat("\n--- A6: interval exposure and landmark ---\n")
  print(as.data.frame(iv_tab))
  print(lm_tab[, c("outcome", "model", "n", "n_events", "HR_fmt", "p_fmt")])
  cat("\n--- A7: duplicate-record sensitivity ---\n")
  print(a7_tab[, c("cohort", "outcome", "measure", "n", "est_fmt", "p_fmt")])
  cat("\n--- A8: person-time contributing sample ---\n"); print(pt_rows)
  cat("\n--- A9: categorical E-value benchmark ---\n"); print(ev_bench[, c("item", "value")])
  cat("\n--- A9: age attenuation by endpoint ---\n"); print(att)
  cat("\n--- A10: competing risk, sparse endpoints ---\n")
  print(a10_tab[, c("endpoint", "method", "sex", "timepoint", "est_fmt", "p_fmt")])
}, file = file.path(DIR_REPX, "qa_sensitivity_05_results.txt"))

capture.output({
  cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Script: scripts/05_qa_sensitivity.R\n")
  cat("Cohort:", COHORT_LABEL, "\n")
  cat("set.seed(42)\n\n")
  print(sessionInfo())
}, file = file.path(DIR_REPX, "sessionInfo_05_qa_sensitivity.txt"))

message("== DONE (05_qa_sensitivity.R) ==")
