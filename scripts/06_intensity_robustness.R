###############################################################################
# 06_intensity_robustness.R — Phase 2.5 재검토 요구 분석 (A16 · A17)
#
# 실행:  Rscript scripts/06_intensity_robustness.R
# 입력:  data/processed/kcto_sex_analysis.rds  (via scripts/00_cohort_prep.R)
# 출력:  output/tables/table15c_intensity_score_robustness.csv
#        output/tables/table15d_intensity_score_decomposition.csv
#        output/tables/tableS13b_intensity_score_common_support_balance.csv
#        output/tables/table19b_landmark_parsimonious.csv
#        interim_report/phase2_stats/reports/qa_intensity_06_log.txt
#        interim_report/phase2_stats/reports/qa_intensity_06_sessionInfo.txt
#
# 코호트: flag_cto_pci & index_yr 2014-2023 -> N = 6,777 (CLAUDE.md §1)
# Primary outcome: DMS (death + MI + stroke)
#
# ⚠️ 관찰연구다. 모든 추정치는 연관성(association)이며 인과 해석 금지.
#    "explained by" 류의 매개 함의 표현도 사용하지 않는다.
# ⛔ 사전 계획을 시사하는 표현 금지 (CLAUDE.md §1-2). 본 분석은 사전등록되지 않았다.
#    아래의 판정 규칙은 QA 가 본 분석 실행 이전에 문서에 기재해 둔 것이며,
#    결과를 본 뒤 수정하지 않는다 (03_required_additional_analyses.md v2 §2).
#
# 목적 — 복합 procedural intensity score (0-4) 는 지금까지 A1(센터)·A2(다중비교)
#        두 조건에서만 검정되었다. 그 4개 성분 중 3개(retrograde / IVUS /
#        강력 P2Y12)는 각각 체격·센터·라벨 조건에서 개별적으로 무효화되었는데
#        점수 모형에는 신장·체중이 없고 라벨층·적응증층 제한도 없다.
#        A16 은 개별 성분을 무너뜨린 바로 그 도전들을 복합 점수에 그대로 부친다.
#
#   A16-0  참조 재현 + 순서형 비례승산 가정 검정
#   A16-1  Leave-one-out (성분 견인 규명) + 가법 분해
#   A16-2  체격 (신장/체중 공변량, 신장 공통지지 제한)
#   A16-3  prasugrel 라벨 비제한 하위군
#   A16-4  적응증/복잡도 baseline 층 (J-CTO >=2, 폐색길이 >=20 mm)
#   A16-5  점수 구성 민감도 (z-score, 이분화, 확장 6성분, 2성분 core)
#   A17    3년 landmark 절약형 공변량 모형 + 모형별 EPP
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(lme4)
  library(MASS)
  library(ordinal)
})

set.seed(42)

PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "00_cohort_prep.R"))

DIR_TAB  <- file.path(PROJ, "output", "tables")
DIR_REPX <- file.path(PROJ, "interim_report", "phase2_stats", "reports")
dir.create(DIR_TAB,  recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_REPX, recursive = TRUE, showWarnings = FALSE)

co <- build_cohort()
stopifnot(nrow(co) == 6777)

# 주 분석 표본 (Model 4 complete case) — 04/05 스크립트와 동일
an <- co[stats::complete.cases(co[, COV_M4, drop = FALSE]), , drop = FALSE]
message(sprintf("[analysis set] Model 4 complete case: N = %d (F = %d)", nrow(an), sum(an$female)))

LOG <- character(0)
logmsg <- function(...) { s <- sprintf(...); LOG <<- c(LOG, s); message(s) }

###############################################################################
# 0. 공통 helper
###############################################################################

# 모형에 실제 투입되는 complete-case 부분집합 (+ 사용되지 않는 factor level 정리)
.cc_subset <- function(dat, yvar, covs) {
  vars <- unique(c(yvar, "female", "center_anon", covs))
  vars <- intersect(vars, names(dat))
  d <- dat[stats::complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE]
  droplevels(d)
}

# table15 와 동일한 컬럼 규격 + A16 전용 메타 컬럼
sc_row <- function(challenge, score_lab, model_lab, d, measure, est, lo, hi, p,
                   stratum, cov_lab, missing_handling, note = "", est_fmt = NULL) {
  data.frame(
    analysis = challenge, outcome = score_lab, model = model_lab,
    n = nrow(d), n_events = NA_integer_,
    n_male = sum(d$female == 0), n_female = sum(d$female == 1),
    events_male = NA_integer_, events_female = NA_integer_,
    pct_male = NA_real_, pct_female = NA_real_,
    measure = measure,
    estimate = est, lo = lo, hi = hi,
    est_fmt = if (!is.null(est_fmt)) est_fmt else sprintf("%.2f (%.2f-%.2f)", est, lo, hi),
    p_value = p, p_fmt = fmt_p(p), note = note,
    stratum = stratum, covariate_set = cov_lab, missing_handling = missing_handling,
    mean_score_male   = round(mean(d$.s[d$female == 0]), 3),
    mean_score_female = round(mean(d$.s[d$female == 1]), 3),
    stringsAsFactors = FALSE, row.names = NULL)
}

# 점수(계수형) 3종 모형: 선형혼합(센터 RI) / Poisson GLMM(센터 RI) / 순서형(cluster-robust)
score_models <- function(dat, scorevar, covs, challenge, score_lab, stratum,
                         cov_lab, missing_handling, note = "",
                         ordinal = TRUE, nagq = 0L) {
  d <- .cc_subset(dat, scorevar, covs)
  d$.s <- d[[scorevar]]
  if (nrow(d) < 100 || length(unique(d$center_anon)) < 3) return(NULL)
  rhs <- if (length(covs)) paste(c("female", covs), collapse = " + ") else "female"
  out <- list()

  # (1) 선형혼합모형 — 보정 평균 점수차 (센터 랜덤절편)  ★주 지표
  m1 <- tryCatch(lme4::lmer(stats::as.formula(paste0(".s ~ ", rhs, " + (1 | center_anon)")),
                            data = d, REML = TRUE), error = function(e) NULL)
  if (!is.null(m1)) {
    s <- summary(m1)
    b <- s$coefficients["female", "Estimate"]; se <- s$coefficients["female", "Std. Error"]
    p <- 2 * stats::pnorm(-abs(b / se))
    out[[length(out) + 1]] <- sc_row(
      challenge, score_lab, "Linear mixed model, centre random intercept", d,
      "adjusted mean difference in score (female - male)",
      b, b - 1.96 * se, b + 1.96 * se, p, stratum, cov_lab, missing_handling, note,
      est_fmt = sprintf("%.3f (%.3f to %.3f)", b, b - 1.96 * se, b + 1.96 * se))
  }

  # (2) Poisson GLMM (센터 랜덤절편) -> RR ("보정 평균 점수비"로만 해석)
  m2 <- tryCatch(lme4::glmer(stats::as.formula(paste0(".s ~ ", rhs, " + (1 | center_anon)")),
                             data = d, family = stats::poisson(), nAGQ = nagq,
                             control = lme4::glmerControl(optimizer = "bobyqa",
                                                          optCtrl = list(maxfun = 3e5))),
                 error = function(e) NULL)
  if (!is.null(m2)) {
    s <- summary(m2)
    b <- s$coefficients["female", "Estimate"]; se <- s$coefficients["female", "Std. Error"]
    sig2 <- as.numeric(lme4::VarCorr(m2)$center_anon[1])
    out[[length(out) + 1]] <- sc_row(
      challenge, score_lab, "Poisson GLMM, centre random intercept", d,
      "RR (ratio of adjusted mean scores)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se),
      s$coefficients["female", "Pr(>|z|)"], stratum, cov_lab, missing_handling,
      paste0(note, sprintf(" | centre random-intercept variance = %.4f, nAGQ = %d", sig2, nagq)))
  }

  # (3) 순서형 로지스틱 (비례오즈), cluster-robust SE
  if (ordinal && length(unique(d$.s)) >= 3) {
    d$.so <- factor(d$.s, ordered = TRUE)
    m3 <- tryCatch(MASS::polr(stats::as.formula(paste0(".so ~ ", rhs)), data = d,
                              Hess = TRUE, method = "logistic"), error = function(e) NULL)
    if (!is.null(m3)) {
      vcp <- tryCatch(sandwich::vcovCL(m3, cluster = d$center_anon),
                      error = function(e) stats::vcov(m3))
      b <- coef(m3)["female"]; se <- sqrt(vcp["female", "female"])
      p <- 2 * stats::pnorm(-abs(b / se))
      out[[length(out) + 1]] <- sc_row(
        challenge, score_lab, "Ordinal logistic (proportional odds), cluster-robust SE", d,
        "OR (one-category-higher intensity)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se),
        p, stratum, cov_lab, missing_handling, note)
    }
  }
  dplyr::bind_rows(out)
}

# 연속 점수(z-score)용: 선형혼합 + OLS cluster-robust
cont_models <- function(dat, scorevar, covs, challenge, score_lab, stratum,
                        cov_lab, missing_handling, note = "") {
  d <- .cc_subset(dat, scorevar, covs); d$.s <- d[[scorevar]]
  if (nrow(d) < 100) return(NULL)
  rhs <- if (length(covs)) paste(c("female", covs), collapse = " + ") else "female"
  out <- list()
  m1 <- tryCatch(lme4::lmer(stats::as.formula(paste0(".s ~ ", rhs, " + (1 | center_anon)")),
                            data = d, REML = TRUE), error = function(e) NULL)
  if (!is.null(m1)) {
    s <- summary(m1); b <- s$coefficients["female", "Estimate"]; se <- s$coefficients["female", "Std. Error"]
    out[[length(out) + 1]] <- sc_row(
      challenge, score_lab, "Linear mixed model, centre random intercept", d,
      "adjusted mean difference (standard deviation units)",
      b, b - 1.96 * se, b + 1.96 * se, 2 * stats::pnorm(-abs(b / se)),
      stratum, cov_lab, missing_handling, note,
      est_fmt = sprintf("%.3f (%.3f to %.3f)", b, b - 1.96 * se, b + 1.96 * se))
  }
  m2 <- stats::lm(stats::as.formula(paste0(".s ~ ", rhs)), data = d)
  ct <- lmtest::coeftest(m2, vcov. = sandwich::vcovCL(m2, cluster = d$center_anon))
  b <- ct["female", "Estimate"]; se <- ct["female", "Std. Error"]
  out[[length(out) + 1]] <- sc_row(
    challenge, score_lab, "Linear model, cluster-robust SE", d,
    "adjusted mean difference (standard deviation units)",
    b, b - 1.96 * se, b + 1.96 * se, ct["female", "Pr(>|t|)"],
    stratum, cov_lab, missing_handling, note,
    est_fmt = sprintf("%.3f (%.3f to %.3f)", b, b - 1.96 * se, b + 1.96 * se))
  dplyr::bind_rows(out)
}

# 이분화 점수용: 센터 랜덤절편 로지스틱
bin_models <- function(dat, yvar, covs, challenge, score_lab, stratum,
                       cov_lab, missing_handling, note = "") {
  d <- .cc_subset(dat, yvar, covs); d$.s <- d[[yvar]]
  if (nrow(d) < 100 || length(unique(d$.s)) < 2) return(NULL)
  rhs <- if (length(covs)) paste(c("female", covs), collapse = " + ") else "female"
  m <- tryCatch(lme4::glmer(stats::as.formula(paste0(".s ~ ", rhs, " + (1 | center_anon)")),
                            data = d, family = stats::binomial(), nAGQ = 0L,
                            control = lme4::glmerControl(optimizer = "bobyqa",
                                                         optCtrl = list(maxfun = 3e5))),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  s <- summary(m); b <- s$coefficients["female", "Estimate"]; se <- s$coefficients["female", "Std. Error"]
  sig2 <- as.numeric(lme4::VarCorr(m)$center_anon[1])
  r <- sc_row(challenge, score_lab, "Logistic GLMM, centre random intercept", d,
              "OR", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se),
              s$coefficients["female", "Pr(>|z|)"], stratum, cov_lab, missing_handling,
              paste0(note, sprintf(" | centre random-intercept variance = %.4f", sig2)))
  r$events_male   <- sum(d$.s == 1 & d$female == 0)
  r$events_female <- sum(d$.s == 1 & d$female == 1)
  r$pct_male   <- round(100 * r$events_male / r$n_male, 1)
  r$pct_female <- round(100 * r$events_female / r$n_female, 1)
  r
}

###############################################################################
# 1. 점수 구성 — 성분 정의 및 점수 생성기
###############################################################################

COMP <- list(
  bilateral  = c(var = "Bilati_inj",        lab = "bilateral injection"),
  retrograde = c(var = "retro_attempt_any", lab = "retrograde attempted"),
  ivus       = c(var = "IVUS",              lab = "IVUS"),
  potent     = c(var = "dapt_potent",       lab = "potent P2Y12 at discharge")
)
IC <- unname(vapply(COMP, function(z) z[["var"]], character(1)))

# mode = "cc": 성분이 하나라도 결측이면 점수 결측 / "mz": 결측을 미시행(0)으로 간주
mk_score <- function(dat, vars, mode = c("cc", "mz")) {
  mode <- match.arg(mode)
  M <- as.data.frame(dat[, vars, drop = FALSE])
  if (mode == "cc") ifelse(stats::complete.cases(M), rowSums(M), NA_real_)
  else rowSums(sapply(vars, function(v) ifelse(is.na(dat[[v]]), 0, dat[[v]])))
}

an$pis_cc <- mk_score(an, IC, "cc")
an$pis_mz <- mk_score(an, IC, "mz")

# leave-one-out 점수 (0-3) — cc / mz 양쪽
LOO <- lapply(names(COMP), function(k) {
  keep <- setdiff(names(COMP), k)
  list(drop = k, drop_lab = COMP[[k]][["lab"]],
       vars = unname(vapply(COMP[keep], function(z) z[["var"]], character(1))),
       lab  = paste0("Intensity score 0-3, excluding ", COMP[[k]][["lab"]]))
})
names(LOO) <- names(COMP)
for (k in names(LOO)) {
  an[[paste0("loo_cc_", k)]] <- mk_score(an, LOO[[k]]$vars, "cc")
  an[[paste0("loo_mz_", k)]] <- mk_score(an, LOO[[k]]$vars, "mz")
}

# mz 처리에서 실제로 potent 제제를 받았음에도 0 으로 오분류되는 환자 수 (표 note 용)
n_mz_misclass <- sum(an$dc_antiplt_2 %in% c("single", "triple") &
                       ((!is.na(an$dc_prasugrel)  & an$dc_prasugrel  == 1) |
                          (!is.na(an$dc_ticagrelor) & an$dc_ticagrelor == 1)))
logmsg("[A16] missing=0 handling misclassifies %d patients who did receive a potent P2Y12 outside a two-drug regimen (single/triple therapy) as score 0 on that component",
       n_mz_misclass)

MZ_NOTE <- sprintf(paste0("Missing-as-not-performed handling misclassifies %d patients on ",
                          "single-agent or triple antiplatelet therapy who did receive a potent ",
                          "P2Y12 inhibitor as not having received one."), n_mz_misclass)

COV_FULL <- COV_M4_NOCENTER                      # Model 4 minus centre volume
COV_PARSIM_PROC <- c("Age", "DM", "CKD_2", "JCTO_re", "Rentrop3", "CTO_location_pm",
                     "Vessel.Disease_f", "index_half_f")   # A3 과 동일한 절약형 세트
LAB_FULL   <- "Model 4 covariates minus centre volume (absorbed by centre)"
LAB_PARSIM <- "Parsimonious: age, DM, CKD, J-CTO, Rentrop 3, prox/mid CTO, vessel disease, index half"

rows <- list()
add <- function(x) if (!is.null(x) && nrow(x)) rows[[length(rows) + 1]] <<- x

###############################################################################
# A16-0. 참조 재현 + 순서형 비례승산 가정 검정
###############################################################################

message("== A16-0: reference reproduction and proportional-odds assumption ==")

add(score_models(an, "pis_cc", COV_FULL, "A16-0",
                 "Intensity score 0-4 (4 components)", "Full analysis cohort",
                 LAB_FULL, "Complete case (all 4 components observed)",
                 "reference model reproduced from table15", nagq = 0L))
add(score_models(an, "pis_cc", COV_FULL, "A16-0",
                 "Intensity score 0-4 (4 components)", "Full analysis cohort",
                 LAB_FULL, "Complete case (all 4 components observed)",
                 "equivalence check: full Laplace approximation", nagq = 1L))
add(score_models(an, "pis_mz", COV_FULL, "A16-0",
                 "Intensity score 0-4 (4 components)", "Full analysis cohort",
                 LAB_FULL, "Missing component treated as not performed",
                 MZ_NOTE, nagq = 0L))
add(score_models(an, "pis_cc", COV_PARSIM_PROC, "A16-0",
                 "Intensity score 0-4 (4 components)", "Full analysis cohort",
                 LAB_PARSIM, "Complete case (all 4 components observed)", "", nagq = 0L))

# --- 비례승산 가정 검정 (Brant 미설치 -> ordinal::clm 부분비례승산 LR 대조) ----
po_test <- do.call(rbind, lapply(c("pis_cc", "pis_mz"), function(v) {
  d <- .cc_subset(an, v, COV_FULL); d$.so <- factor(d[[v]], ordered = TRUE)
  rhs <- paste(c("female", COV_FULL), collapse = " + ")
  m0 <- tryCatch(ordinal::clm(stats::as.formula(paste0(".so ~ ", rhs)), data = d),
                 error = function(e) NULL)
  m1 <- tryCatch(ordinal::clm(stats::as.formula(paste0(".so ~ ", rhs)),
                              nominal = ~ female, data = d), error = function(e) NULL)
  if (is.null(m0) || is.null(m1)) return(NULL)
  a <- stats::anova(m0, m1)
  data.frame(score_definition = ifelse(v == "pis_cc",
                                       "Complete case (all 4 components observed)",
                                       "Missing component treated as not performed"),
             n = nrow(d),
             test = paste0("Likelihood-ratio test of the proportional-odds assumption for sex ",
                           "(cumulative link model with sex-specific thresholds vs common thresholds)"),
             LR_statistic = round(a$LR.stat[2], 3), df = a$df[2],
             p_value = a$`Pr(>Chisq)`[2], p_fmt = fmt_p(a$`Pr(>Chisq)`[2]),
             interpretation = ifelse(a$`Pr(>Chisq)`[2] < 0.05,
                                     "proportional-odds assumption rejected for sex",
                                     "no evidence against the proportional-odds assumption for sex"),
             stringsAsFactors = FALSE)
}))
logmsg("[A16-0] proportional-odds LR test (sex): %s",
       paste(sprintf("%s p=%s", po_test$score_definition, po_test$p_fmt), collapse = " | "))

###############################################################################
# A16-1. Leave-one-out — 어느 성분이 결과를 견인하는가
#
# (a) 표본 고정 (4성분 모두 관측된 n=5,097): 점수 구성만 바뀌므로 성분 기여를
#     직접 비교할 수 있다.
# (b) 각 점수의 자체 complete case: 실제로 그 점수를 보고할 때의 표본.
###############################################################################

message("== A16-1: leave-one-out ==")

an_fix <- an[!is.na(an$pis_cc), , drop = FALSE]     # 표본 고정용
logmsg("[A16-1] fixed comparison sample (all 4 components observed): N = %d (M %d / F %d)",
       nrow(an_fix), sum(an_fix$female == 0), sum(an_fix$female == 1))

for (k in names(LOO)) {
  # (a) 표본 고정
  add(score_models(an_fix, paste0("loo_cc_", k), COV_FULL, "A16-1",
                   LOO[[k]]$lab, "Full analysis cohort",
                   LAB_FULL, "Complete case, sample fixed to the 4-component complete case",
                   sprintf("leave-one-out: %s removed; sample held fixed at n=%d so that only the score composition changes",
                           LOO[[k]]$drop_lab, nrow(an_fix)), nagq = 0L))
  # (b) 자체 complete case
  add(score_models(an, paste0("loo_cc_", k), COV_FULL, "A16-1",
                   LOO[[k]]$lab, "Full analysis cohort",
                   LAB_FULL, "Complete case for the three retained components",
                   sprintf("leave-one-out: %s removed; own complete case", LOO[[k]]$drop_lab),
                   nagq = 0L))
  # (c) 결측=0
  add(score_models(an, paste0("loo_mz_", k), COV_FULL, "A16-1",
                   LOO[[k]]$lab, "Full analysis cohort",
                   LAB_FULL, "Missing component treated as not performed",
                   sprintf("leave-one-out: %s removed. %s", LOO[[k]]$drop_lab, MZ_NOTE),
                   nagq = 0L))
}

# --- 가법 분해 (crude + 보정) ------------------------------------------------
# 점수는 성분의 합이므로, 동일 표본·동일 설계행렬의 선형모형에서 성분별 female
# 계수의 합은 총 점수의 female 계수와 정확히 일치한다(선형성). 센터는 고정효과로
# 넣어 within-centre 대비로 만든다. 이 항등식은 assertion 으로 확인한다.
message("== A16-1b: additive decomposition of the score gap ==")

decomp_one <- function(dat, comps, tag, missing_handling, labs = NULL) {
  d <- dat
  stopifnot(!any(is.na(d[, comps])))     # 가법 항등식은 동일 표본에서만 성립
  d$.s <- rowSums(sapply(comps, function(v) d[[v]]))
  rhs_adj <- paste(c("female", COV_FULL, "factor(center_anon)"), collapse = " + ")
  fit1 <- function(y) {
    d$.y <- y
    m <- stats::lm(stats::as.formula(paste0(".y ~ ", rhs_adj)), data = d)
    ct <- lmtest::coeftest(m, vcov. = sandwich::vcovCL(m, cluster = d$center_anon))
    c(b = unname(ct["female", "Estimate"]), se = unname(ct["female", "Std. Error"]),
      p = unname(ct["female", "Pr(>|t|)"]))
  }
  tot <- fit1(d$.s)
  parts <- t(sapply(comps, function(v) fit1(d[[v]])))
  stopifnot(abs(sum(parts[, "b"]) - tot["b"]) < 1e-8)   # 가법 항등식 확인
  crude <- sapply(comps, function(v) mean(d[[v]][d$female == 1]) - mean(d[[v]][d$female == 0]))
  crude_tot <- mean(d$.s[d$female == 1]) - mean(d$.s[d$female == 0])
  out <- data.frame(
    decomposition = tag, missing_handling = missing_handling, n = nrow(d),
    n_male = sum(d$female == 0), n_female = sum(d$female == 1),
    component = c(if (!is.null(labs)) labs else vapply(comps, function(v) {
      i <- which(IC == v); if (length(i)) COMP[[i]][["lab"]] else v }, character(1)),
      "TOTAL (sum of components)"),
    pct_male   = round(100 * c(sapply(comps, function(v) mean(d[[v]][d$female == 0])), NA), 1),
    pct_female = round(100 * c(sapply(comps, function(v) mean(d[[v]][d$female == 1])), NA), 1),
    crude_diff_pp = round(100 * c(crude, crude_tot), 2),
    crude_share_pct = round(100 * c(crude, crude_tot) / crude_tot * 100 / 100, 1),
    adj_diff_pp = round(100 * c(parts[, "b"], tot["b"]), 2),
    adj_lo_pp = round(100 * c(parts[, "b"] - 1.96 * parts[, "se"], tot["b"] - 1.96 * tot["se"]), 2),
    adj_hi_pp = round(100 * c(parts[, "b"] + 1.96 * parts[, "se"], tot["b"] + 1.96 * tot["se"]), 2),
    adj_p = c(parts[, "p"], tot["p"]),
    adj_share_pct = round(100 * c(parts[, "b"], tot["b"]) / tot["b"], 1),
    stringsAsFactors = FALSE, row.names = NULL)
  out$adj_p_fmt <- fmt_p(out$adj_p)
  out$adj_fmt <- sprintf("%.2f (%.2f to %.2f)", out$adj_diff_pp, out$adj_lo_pp, out$adj_hi_pp)
  out
}
IC_MZ <- paste0("mz_", IC)
for (i in seq_along(IC)) an[[IC_MZ[i]]] <- ifelse(is.na(an[[IC[i]]]), 0, an[[IC[i]]])
COMP_LABS <- unname(vapply(COMP, function(z) z[["lab"]], character(1)))

dec_tab <- dplyr::bind_rows(
  decomp_one(.cc_subset(an_fix, "pis_cc", COV_FULL), IC,
             "4-component score, adjusted (Model 4 covariates + centre fixed effects, cluster-robust SE)",
             "Complete case (all 4 components observed)", labs = COMP_LABS),
  decomp_one(.cc_subset(an, "pis_mz", COV_FULL), IC_MZ,
             "4-component score, adjusted (Model 4 covariates + centre fixed effects, cluster-robust SE)",
             "Missing component treated as not performed", labs = COMP_LABS))
dec_tab$note <- paste0(
  "Percentage-point difference in the probability of each component for women versus men. ",
  "Because the score is the sum of its components and the models are linear on a common design ",
  "matrix, the component coefficients sum exactly to the coefficient for the total score; this ",
  "identity is verified in code. Shares are relative to the total gap and are descriptive only. ",
  "Observational data: all estimates are associations. ")
dec_tab$note[dec_tab$missing_handling != "Complete case (all 4 components observed)"] <-
  paste0(dec_tab$note[dec_tab$missing_handling != "Complete case (all 4 components observed)"], MZ_NOTE)
write.csv(dec_tab, file.path(DIR_TAB, "table15d_intensity_score_decomposition.csv"),
          row.names = FALSE, na = "")
logmsg("[A16-1b] adjusted decomposition (complete case): %s",
       paste(sprintf("%s %.2f pp (%.0f%%)",
                     dec_tab$component[dec_tab$missing_handling == "Complete case (all 4 components observed)"],
                     dec_tab$adj_diff_pp[dec_tab$missing_handling == "Complete case (all 4 components observed)"],
                     dec_tab$adj_share_pct[dec_tab$missing_handling == "Complete case (all 4 components observed)"]),
             collapse = " | "))

###############################################################################
# A16-2. 체격 (A3 대응)
#   ① 전 코호트 + 신장 공변량  ② + 체중  ③ + 신장·체중
#   ④ strict 공통지지  ⑤ wide 공통지지  (A3 와 동일한 창)
###############################################################################

message("== A16-2: body size ==")

h <- co$Height
CS_STRICT <- c(unname(stats::quantile(h[co$female == 0], .05,  na.rm = TRUE)),
               unname(stats::quantile(h[co$female == 1], .95,  na.rm = TRUE)))
CS_WIDE   <- c(unname(stats::quantile(h[co$female == 0], .025, na.rm = TRUE)),
               unname(stats::quantile(h[co$female == 1], .975, na.rm = TRUE)))
logmsg("[A16-2] height common support (identical to A3): strict %.0f-%.0f cm, wide %.0f-%.0f cm",
       CS_STRICT[1], CS_STRICT[2], CS_WIDE[1], CS_WIDE[2])
in_cs <- function(d, rng) !is.na(d$Height) & d$Height >= rng[1] & d$Height <= rng[2]

smd_cont <- function(x, g) {
  x1 <- x[g == 1]; x0 <- x[g == 0]
  (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) /
    sqrt((stats::var(x1, na.rm = TRUE) + stats::var(x0, na.rm = TRUE)) / 2)
}

# --- balance 표 (점수 분석 표본 기준) ---------------------------------------
bal_rows <- list()
for (cfg in list(
  list(lab = "Score analysis sample, no height restriction", d = an_fix),
  list(lab = sprintf("Common support, strict (%.0f-%.0f cm)", CS_STRICT[1], CS_STRICT[2]),
       d = an_fix[in_cs(an_fix, CS_STRICT), ]),
  list(lab = sprintf("Common support, wide (%.0f-%.0f cm)", CS_WIDE[1], CS_WIDE[2]),
       d = an_fix[in_cs(an_fix, CS_WIDE), ]))) {
  dd <- cfg$d
  for (v in c("Height", "Weight", "bsa", "BMI_re", "Age", "JCTO_re", "pis_cc")) {
    if (!v %in% names(dd)) next
    bal_rows[[length(bal_rows) + 1]] <- data.frame(
      restriction = cfg$lab, variable = v,
      n_male = sum(dd$female == 0 & !is.na(dd[[v]])),
      n_female = sum(dd$female == 1 & !is.na(dd[[v]])),
      mean_male = round(mean(dd[[v]][dd$female == 0], na.rm = TRUE), 3),
      mean_female = round(mean(dd[[v]][dd$female == 1], na.rm = TRUE), 3),
      SMD = round(smd_cont(dd[[v]], dd$female), 3),
      abs_SMD = round(abs(smd_cont(dd[[v]], dd$female)), 3),
      stringsAsFactors = FALSE)
  }
}
bal_cs <- dplyr::bind_rows(bal_rows)
bal_cs$balance_achieved <- ifelse(bal_cs$abs_SMD < 0.10, "yes (|SMD|<0.10)",
                                  ifelse(bal_cs$abs_SMD < 0.25, "acceptable (|SMD|<0.25)", "no"))
bal_cs$note <- paste0("Common-support windows are identical to those used for the individual ",
                      "management-strategy indicators (table16 / tableS13). Restriction on height ",
                      "does not remove the residual height imbalance.")
write.csv(bal_cs, file.path(DIR_TAB, "tableS13b_intensity_score_common_support_balance.csv"),
          row.names = FALSE, na = "")

EXTRAP <- sprintf(paste0("EXTRAPOLATION WARNING: height |SMD| = %.2f and weight |SMD| = %.2f ",
                         "between sexes; models that adjust for body size rely on the parametric ",
                         "form outside the region of common support."),
                  abs(smd_cont(an_fix$Height, an_fix$female)),
                  abs(smd_cont(an_fix$Weight, an_fix$female)))

BODY <- list(
  list(tag = "Full cohort + height as covariate", covs_add = "Height",  d = NULL),
  list(tag = "Full cohort + weight as covariate", covs_add = "Weight",  d = NULL),
  list(tag = "Full cohort + height and weight as covariates",
       covs_add = c("Height", "Weight"), d = NULL)
)
for (b in BODY) {
  for (cs in list(list(cv = COV_FULL, lab = LAB_FULL),
                  list(cv = COV_PARSIM_PROC, lab = LAB_PARSIM))) {
    add(score_models(an, "pis_cc", c(cs$cv, b$covs_add), "A16-2",
                     "Intensity score 0-4 (4 components)", b$tag,
                     paste0(cs$lab, " + ", paste(b$covs_add, collapse = " and ")),
                     "Complete case (all 4 components observed)", EXTRAP, nagq = 0L))
  }
  # 3성분(potent 제외) 점수도 동일하게
  add(score_models(an, "loo_cc_potent", c(COV_FULL, b$covs_add), "A16-2",
                   LOO$potent$lab, b$tag,
                   paste0(LAB_FULL, " + ", paste(b$covs_add, collapse = " and ")),
                   "Complete case for the three retained components", EXTRAP, nagq = 0L))
}
for (cfg in list(list(tag = sprintf("Common support, strict (%.0f-%.0f cm)", CS_STRICT[1], CS_STRICT[2]),
                      d = an[in_cs(an, CS_STRICT), ]),
                 list(tag = sprintf("Common support, wide (%.0f-%.0f cm)", CS_WIDE[1], CS_WIDE[2]),
                      d = an[in_cs(an, CS_WIDE), ]))) {
  logmsg("[A16-2] %s: N = %d (M %d / F %d)", cfg$tag, nrow(cfg$d),
         sum(cfg$d$female == 0), sum(cfg$d$female == 1))
  for (cs in list(list(cv = COV_FULL, lab = LAB_FULL),
                  list(cv = COV_PARSIM_PROC, lab = LAB_PARSIM))) {
    add(score_models(cfg$d, "pis_cc", cs$cv, "A16-2",
                     "Intensity score 0-4 (4 components)", cfg$tag, cs$lab,
                     "Complete case (all 4 components observed)",
                     "height common support restriction; residual height imbalance persists (tableS13b)",
                     nagq = 0L))
    add(score_models(cfg$d, "pis_mz", cs$cv, "A16-2",
                     "Intensity score 0-4 (4 components)", cfg$tag, cs$lab,
                     "Missing component treated as not performed", MZ_NOTE, nagq = 0L))
  }
  add(score_models(cfg$d, "loo_cc_potent", COV_FULL, "A16-2",
                   LOO$potent$lab, cfg$tag, LAB_FULL,
                   "Complete case for the three retained components",
                   "height common support restriction", nagq = 0L))
  add(score_models(cfg$d, "loo_cc_bilateral", COV_FULL, "A16-2",
                   LOO$bilateral$lab, cfg$tag, LAB_FULL,
                   "Complete case for the three retained components",
                   "height common support restriction", nagq = 0L))
}

###############################################################################
# A16-3. prasugrel 라벨 (A4 대응)
#   라벨 제한(>=75세 또는 <60 kg) 비해당 하위군에서 4성분·3성분 점수 재적합
###############################################################################

message("== A16-3: prasugrel label-unrestricted subgroup ==")

an$label_free <- !is.na(an$Age) & !is.na(an$Weight) & an$Age < 75 & an$Weight >= 60
an_lf <- an[an$label_free, , drop = FALSE]
logmsg("[A16-3] label-unrestricted subgroup (age <75 y AND weight >=60 kg): N = %d (M %d / F %d)",
       nrow(an_lf), sum(an_lf$female == 0), sum(an_lf$female == 1))
LF_TAG <- sprintf("Prasugrel label-unrestricted subgroup (age <75 y and weight >=60 kg), n=%d (M %d / F %d)",
                  nrow(an_lf), sum(an_lf$female == 0), sum(an_lf$female == 1))
LF_NOTE <- paste0("Women in this subgroup are a minority of all women in the cohort; ",
                  "an absence of a detectable difference here is a limit of resolution, ",
                  "not evidence of no difference.")

for (cs in list(list(cv = COV_FULL, lab = LAB_FULL),
                list(cv = COV_PARSIM_PROC, lab = LAB_PARSIM))) {
  add(score_models(an_lf, "pis_cc", cs$cv, "A16-3",
                   "Intensity score 0-4 (4 components)", LF_TAG, cs$lab,
                   "Complete case (all 4 components observed)", LF_NOTE, nagq = 0L))
  add(score_models(an_lf, "loo_cc_potent", cs$cv, "A16-3",
                   LOO$potent$lab, LF_TAG, cs$lab,
                   "Complete case for the three retained components", LF_NOTE, nagq = 0L))
}
add(score_models(an_lf, "pis_mz", COV_FULL, "A16-3",
                 "Intensity score 0-4 (4 components)", LF_TAG, LAB_FULL,
                 "Missing component treated as not performed",
                 paste(LF_NOTE, MZ_NOTE), nagq = 0L))
# 라벨 제한 해당군 (대조)
an_lr <- an[!is.na(an$Age) & !is.na(an$Weight) & (an$Age >= 75 | an$Weight < 60), , drop = FALSE]
LR_TAG <- sprintf("Prasugrel label-restricted subgroup (age >=75 y or weight <60 kg), n=%d (M %d / F %d)",
                  nrow(an_lr), sum(an_lr$female == 0), sum(an_lr$female == 1))
logmsg("[A16-3] label-restricted subgroup: N = %d (M %d / F %d)", nrow(an_lr),
       sum(an_lr$female == 0), sum(an_lr$female == 1))
add(score_models(an_lr, "pis_cc", COV_FULL, "A16-3",
                 "Intensity score 0-4 (4 components)", LR_TAG, LAB_FULL,
                 "Complete case (all 4 components observed)",
                 "contrast stratum for the label analysis", nagq = 0L))
add(score_models(an_lr, "loo_cc_potent", COV_FULL, "A16-3",
                 LOO$potent$lab, LR_TAG, LAB_FULL,
                 "Complete case for the three retained components",
                 "contrast stratum for the label analysis", nagq = 0L))

###############################################################################
# A16-4. 적응증 / 복잡도 층 (A5 대응)
#   baseline 층(J-CTO >=2, 폐색길이 >=20 mm)이 주 근거. Fail_tech 층은
#   노출 이후 변수 조건화(collider)이므로 descriptive only 로 별도 표기.
###############################################################################

message("== A16-4: indication and complexity strata ==")

STRATA <- list(
  list(tag = "Baseline complexity stratum: J-CTO score >=2",
       d = an[!is.na(an$JCTO_re) & an$JCTO_re >= 2, , drop = FALSE], desc_only = FALSE),
  list(tag = "Baseline complexity stratum: occlusion length >=20 mm",
       d = an[an$Length20 %in% 1, , drop = FALSE], desc_only = FALSE),
  list(tag = "Baseline complexity stratum: J-CTO score <2 (contrast)",
       d = an[!is.na(an$JCTO_re) & an$JCTO_re < 2, , drop = FALSE], desc_only = FALSE),
  list(tag = "Post-baseline stratum: failure of initial wiring technique (Fail_tech = 1)",
       d = an[an$Fail_tech %in% 1, , drop = FALSE], desc_only = TRUE))
for (st in STRATA) {
  logmsg("[A16-4] %s: N = %d (M %d / F %d)", st$tag, nrow(st$d),
         sum(st$d$female == 0), sum(st$d$female == 1))
  nt <- if (st$desc_only)
    paste0("DESCRIPTIVE CONTRAST ONLY: this stratum is defined by a variable measured after the ",
           "start of the procedure and downstream of the score components; conditioning on it can ",
           "induce or reverse an association (collider). Not used as supporting evidence.")
  else "baseline (pre-procedural) stratum; not conditioned on any post-exposure variable"
  for (cs in list(list(cv = COV_FULL, lab = LAB_FULL),
                  list(cv = COV_PARSIM_PROC, lab = LAB_PARSIM))) {
    add(score_models(st$d, "pis_cc", cs$cv, "A16-4",
                     "Intensity score 0-4 (4 components)", st$tag, cs$lab,
                     "Complete case (all 4 components observed)", nt, nagq = 0L))
  }
  add(score_models(st$d, "loo_cc_potent", COV_FULL, "A16-4",
                   LOO$potent$lab, st$tag, LAB_FULL,
                   "Complete case for the three retained components", nt, nagq = 0L))
  add(score_models(st$d, "loo_cc_bilateral", COV_FULL, "A16-4",
                   LOO$bilateral$lab, st$tag, LAB_FULL,
                   "Complete case for the three retained components", nt, nagq = 0L))
}

###############################################################################
# A16-5. 점수 구성 민감도
#   (a) z-score 평균 (등가중 가정 해제)
#   (b) 이분화 (>=1점, >=2점)
#   (c) 확장 6성분 (미세도관 사용, 와이어 >=3개 추가) — 결측 큼, 민감도 전용
#   (d) 2성분 core (양측조영 + 역행성) / 3성분 순수 시술 (potent 제외)
###############################################################################

message("== A16-5: alternative score constructions ==")

# (a) z-score: 4성분을 각각 표준화(고정 표본 기준) 후 평균
zsc <- sapply(IC, function(v) {
  x <- an_fix[[v]]; (x - mean(x)) / stats::sd(x)
})
an_fix$pis_z <- rowMeans(zsc)
add(cont_models(an_fix, "pis_z", COV_FULL, "A16-5",
                "Intensity score, standardised-component mean (z-score)",
                "Full analysis cohort", LAB_FULL,
                "Complete case (all 4 components observed)",
                paste0("Each component standardised to mean 0 and SD 1 before averaging; ",
                       "removes the assumption that a one-unit change in each component ",
                       "carries the same weight.")))

# (b) 이분화
an_fix$pis_ge1 <- as.numeric(an_fix$pis_cc >= 1)
an_fix$pis_ge2 <- as.numeric(an_fix$pis_cc >= 2)
add(bin_models(an_fix, "pis_ge1", COV_FULL, "A16-5",
               "Intensity score >=1 (any of the four)", "Full analysis cohort",
               LAB_FULL, "Complete case (all 4 components observed)",
               "dichotomised alternative to the ordinal construction"))
add(bin_models(an_fix, "pis_ge2", COV_FULL, "A16-5",
               "Intensity score >=2", "Full analysis cohort",
               LAB_FULL, "Complete case (all 4 components observed)",
               "dichotomised alternative to the ordinal construction"))

# (c) 확장 6성분
an$micro_any  <- ifelse(is.na(an$n.micro.3g), NA_real_, as.numeric(an$n.micro.3g >= 1))
an$wire_ge3   <- ifelse(is.na(an$n.wire.5g),  NA_real_, as.numeric(an$n.wire.5g  >= 3))
IC6 <- c(IC, "micro_any", "wire_ge3")
an$pis6_cc <- mk_score(an, IC6, "cc")
n6 <- sum(!is.na(an$pis6_cc))
logmsg("[A16-5] expanded 6-component score complete case: N = %d (%.1f%% of %d)",
       n6, 100 * n6 / nrow(an), nrow(an))
add(score_models(an, "pis6_cc", COV_FULL, "A16-5",
                 "Expanded intensity score 0-6 (adds microcatheter use and >=3 guidewires)",
                 "Full analysis cohort", LAB_FULL,
                 "Complete case for all six components",
                 sprintf(paste0("Addresses the question of why other procedural-effort markers were ",
                                "not components. Microcatheter count and guidewire count are missing in ",
                                "23.8%% and 30.4%% of the analysis cohort, so this score is available in ",
                                "only %d patients and its missingness may itself differ by sex; ",
                                "sensitivity only."), n6), nagq = 0L))

# (d) 2성분 core / 3성분 순수 시술
an$pis_core2 <- mk_score(an, c("Bilati_inj", "retro_attempt_any"), "cc")
add(score_models(an, "pis_core2", COV_FULL, "A16-5",
                 "Intensity score 0-2 (bilateral injection + retrograde attempt only)",
                 "Full analysis cohort", LAB_FULL,
                 "Complete case for the two retained components",
                 "diagnostic/escalation core; excludes imaging and pharmacotherapy", nagq = 0L))
add(score_models(an, "loo_cc_potent", COV_PARSIM_PROC, "A16-5",
                 LOO$potent$lab, "Full analysis cohort", LAB_PARSIM,
                 "Complete case for the three retained components",
                 "purely procedural score; resolves the naming mismatch of the 4-component score",
                 nagq = 0L))

###############################################################################
# A16 결과 저장
###############################################################################

a16 <- dplyr::bind_rows(rows)
a16$score_components <- "Reference score = bilateral injection + retrograde attempted + IVUS + potent P2Y12 at discharge (0-4), equally weighted."
a16$decision_rule <- paste0(
  "Decision rule recorded by the QA reviewer in interim_report/phase2.5_qa/",
  "03_required_additional_analyses.md (v2, section 2) before these analyses were run: ",
  "(i) if the three-component score excluding potent P2Y12 AND the body-size-adjusted score both ",
  "remain significant, the aggregate claim stands; (ii) if only some conditions hold, the claim is ",
  "downgraded to attenuated but directionally consistent; (iii) if the aggregate is lost under the ",
  "body-size or label condition, the aggregate claim is withdrawn.")
a16$caveat <- paste0("Observational secondary analysis of a registry; estimates are associations. ",
                     "These analyses were not prospectively registered.")
write.csv(a16, file.path(DIR_TAB, "table15c_intensity_score_robustness.csv"),
          row.names = FALSE, na = "")
write.csv(po_test, file.path(DIR_TAB, "table15e_intensity_proportional_odds_test.csv"),
          row.names = FALSE, na = "")
logmsg("[A16] wrote %d model rows to table15c", nrow(a16))

###############################################################################
# A17. 3년 landmark — 절약형 공변량 병기 + 모형별 EPP
#
# 왜: table19 의 보정 landmark 모형은 Model 4(모수 ~26개)를 events 131건(DMS) /
#     108건(사망) 에 적합했다 -> EPP ~ 5. 원고는 다른 곳에서 events/parameter < 10
#     일 때 절약형 모형을 병기하는 방식을 적용했으므로, 이 방어 분석에도 같은
#     방식을 적용해 결과를 대조한다.
###############################################################################

message("== A17: 3-year landmark with parsimonious covariates and EPP ==")

LM <- 1095
landmark_row <- function(key, covs, label, note) {
  d <- make_surv(an, key)
  d <- d[d$.time > LM, , drop = FALSE]        # 3년 시점 event-free 생존자
  d$.time <- d$.time - LM
  d <- .cc_subset(d, ".event", covs)
  m <- tryCatch(cox_model(d, covs), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  s <- summary(m); i <- which(rownames(s$coefficients) == "female")
  npar <- length(stats::coef(m))
  data.frame(
    analysis = "A17", outcome = OUT[[key]]$lab, outcome_key = key, model = label,
    n = m$n, n_events = m$nevent,
    n_female = sum(d$female == 1), events_female = sum(d$.event[d$female == 1] == 1),
    n_male = sum(d$female == 0), events_male = sum(d$.event[d$female == 0] == 1),
    n_parameters = npar,
    events_per_parameter = round(m$nevent / npar, 2),
    HR = unname(s$conf.int[i, 1]), lo = unname(s$conf.int[i, 3]), hi = unname(s$conf.int[i, 4]),
    HR_fmt = sprintf("%.2f (%.2f-%.2f)", s$conf.int[i, 1], s$conf.int[i, 3], s$conf.int[i, 4]),
    se_robust = unname(s$coefficients[i, "robust se"]),
    p_value = unname(s$coefficients[i, "Pr(>|z|)"]),
    p_fmt = fmt_p(unname(s$coefficients[i, "Pr(>|z|)"])),
    note = note, stringsAsFactors = FALSE)
}

LM_SPECS <- list(
  list(covs = character(0), lab = "Landmark at 3 years: unadjusted",
       note = "restricted to patients event-free and still under observation at 3 years"),
  list(covs = COV_AGE, lab = "Landmark at 3 years: age only",
       note = "isolates the contribution of age composition at the landmark"),
  list(covs = COV_PARS, lab = "Landmark at 3 years: parsimonious (age, diabetes, CKD, atrial fibrillation)",
       note = "parsimonious set used elsewhere in this analysis when events per parameter fall below 10"),
  list(covs = COV_M4, lab = "Landmark at 3 years: Model 4 (re-fitted from the 3-year landmark)",
       note = "reference model reported in table19; events per parameter is low, see the EPP column")
)
lm_rows <- list()
for (k in c("DMS", "Death", "MACCE_pub", "MACE_pub", "NACE_pub")) {
  for (sp in LM_SPECS) {
    r <- tryCatch(landmark_row(k, sp$covs, sp$lab, sp$note), error = function(e) NULL)
    if (!is.null(r)) lm_rows[[length(lm_rows) + 1]] <- r
  }
}
lm_tab <- dplyr::bind_rows(lm_rows)
lm_tab$epp_flag <- ifelse(lm_tab$events_per_parameter < 10,
                          "events per parameter <10", "events per parameter >=10")
lm_tab$caveat <- paste0(
  "Landmark analyses condition on being event-free and under observation at 3 years and are ",
  "subject to depletion of susceptibles. Estimates are associations. These analyses were not ",
  "prospectively registered. This table supersedes table19 by adding the parsimonious and ",
  "age-only specifications and the events-per-parameter column.")
write.csv(lm_tab, file.path(DIR_TAB, "table19b_landmark_parsimonious.csv"),
          row.names = FALSE, na = "")

for (k in unique(lm_tab$outcome_key)) {
  ss <- lm_tab[lm_tab$outcome_key == k, ]
  logmsg("[A17] %s: %s", k,
         paste(sprintf("%s HR %s p=%s (EPP %.1f)",
                       sub("Landmark at 3 years: ", "", ss$model), ss$HR_fmt, ss$p_fmt,
                       ss$events_per_parameter), collapse = " | "))
}

###############################################################################
# 로그 · sessionInfo
###############################################################################

writeLines(LOG, file.path(DIR_REPX, "qa_intensity_06_log.txt"))
capture.output(utils::sessionInfo(), file = file.path(DIR_REPX, "qa_intensity_06_sessionInfo.txt"))
message("== 06_intensity_robustness.R done ==")
