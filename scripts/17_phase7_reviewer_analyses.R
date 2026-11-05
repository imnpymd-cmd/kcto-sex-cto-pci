# =============================================================================
# 17_phase7_reviewer_analyses.R
#   Phase 7 (peer review) 추가 분석 — JACC: Advances JADV102226-3033-T
#   편지 수령 2026-11-04 / 본 스크립트 2026-11-05
#
# 실행:  Rscript scripts/17_phase7_reviewer_analyses.R
# 입력:  data/processed/kcto_sex_analysis.rds  (via scripts/00_cohort_prep.R)
#
# 본 스크립트는 **심사 의견에 답하기 위해 추가된** 분석만 담는다 (analysis plan
# 에 없었음). 모든 산출물 행의 caveat 컬럼에 "added in response to peer review"
# 를 기입한다. ⛔ "prespecified / pre-registered" 긍정 주장 금지 (CLAUDE.md §1-2).
#
# 수행 범위 (연구자 결정이 필요 없는 항목만):
#   A-3 (R2-M3)  3성분 procedural intensity score: 분포 · 보정 주변확률 · 성분 분해
#   A-6 (R2-M6)  성분별 성별 누적발생률(사망 경쟁) + 재관류 Fine-Gray
#   A-9 (R2-m4)  g-formula 표준화 보정 절대위험(1·3·5년) + bilateral·성공 주변확률
#   A-8 (R1-W1)  sex x era 상호작용 (3성분 점수 · bilateral)
#   A-7 (R2-M7)  창별 모형 진단(PH·EPP·VIF·연령 선형성) 정리표
#   A-5 (R2-M5)  >3년 창 IPTW 진단(ESS · 가중치 분포 · 사후 SMD · 절단 추정치)
#   A-2 (1)(2)   무추적 환자의 센터 분포 + 해당 센터 제외 landmark 재적합
#   ⛔ 미수행(연구자 정의 대기): A-1, A-4, A-2(3) delta-adjusted tipping point
#
# 코호트: flag_cto_pci & index_yr 2014-2023 -> N = 6,777 (CLAUDE.md §1)
#         Cox 완전관측 Model 4 표본 n = 6,259 / center_anon clustered robust SE
# 창 경계: 365 / 1095일, day-0 사건 0.5일 (기존 규약, 00_cohort_prep.R)
#
# ⛔ PHI: id / id_code / patient_key / Center_name 값은 출력하지 않는다.
#         센터는 익명 정수 center_anon 으로만 보고한다 (00_cohort_prep.R 가 제거).
# ⚠️ 관찰연구다. 모든 추정치는 연관성(association)이며 인과 해석·용어 금지.
#
# 출력: output/tables/table28_* ~ table34_*  (기존 table26/27 이 이미 존재하므로
#       브리프의 "table26 onward" 를 table28 부터로 조정 — 덮어쓰기 방지)
#       interim_report/phase7_rebuttal/reports/phase7_analyses_log.txt
#       interim_report/phase7_rebuttal/reports/phase7_sessionInfo.txt
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(survival); library(lme4); library(MASS)
  library(sandwich); library(lmtest); library(cmprsk)
})

set.seed(42)
SEED <- 42L
B_BOOT <- 500L      # g-formula 부트스트랩 반복
B_MC   <- 1000L     # 고정효과 사후 추출 반복 (GLMM 주변확률 CI)

PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "00_cohort_prep.R"))

DIR_TAB <- file.path(PROJ, "output", "tables")
DIR_REP <- file.path(PROJ, "interim_report", "phase7_rebuttal", "reports")
dir.create(DIR_REP, recursive = TRUE, showWarnings = FALSE)

LOG <- character(0)
logmsg <- function(...) { s <- sprintf(...); LOG <<- c(LOG, s); message(s) }

###############################################################################
# 0. 공통 상수 · helper
###############################################################################

CAV <- paste0(
  "Added in response to peer review (JACC: Advances, reviewer comments received 2026-11-04); ",
  "not part of the analysis set of the submitted manuscript and not prospectively registered. ",
  "Observational secondary analysis of a registry: all estimates are associations.")

# 표준 출력 행 (요청된 공통 컬럼 + 분석별 추가 컬럼)
rw <- function(analysis, outcome, model, n = NA_integer_, n_female = NA_integer_,
               estimate = NA_real_, lo = NA_real_, hi = NA_real_,
               est_fmt = NA_character_, p_fmt = NA_character_,
               note = "", caveat = CAV, ...) {
  base <- data.frame(analysis = analysis, outcome = outcome, model = model,
                     n = n, n_female = n_female,
                     estimate = estimate, lo = lo, hi = hi,
                     est_fmt = est_fmt, p_fmt = p_fmt,
                     note = note, caveat = caveat, stringsAsFactors = FALSE)
  ex <- list(...)
  for (nm in names(ex)) base[[nm]] <- ex[[nm]]
  base
}

fmt3 <- function(e, l, h) sprintf("%.3f (%.3f to %.3f)", e, l, h)
fmt2 <- function(e, l, h) sprintf("%.2f (%.2f-%.2f)", e, l, h)
fmtpp <- function(e, l, h) sprintf("%+.2f pp (%+.2f to %+.2f)", e, l, h)

# 창 분할 (12/13/14 스크립트와 동일)
window_data <- function(dat, key, start, end = Inf) {
  d <- make_surv(dat, key)
  d <- d[d$.time > start, , drop = FALSE]
  d$.time <- d$.time - start
  if (is.finite(end)) {
    w <- end - start
    d$.event <- ifelse(d$.time <= w & d$.event == 1, 1L, 0L)
    d$.time  <- pmin(d$.time, w)
  }
  d
}
cc_subset <- function(dat, covs) {
  vars <- intersect(unique(c(".time", ".event", "female", "center_anon", covs)), names(dat))
  droplevels(dat[stats::complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE])
}
WINDOWS <- list(
  list(start = 0,    end = 365,  lab = "0-1 year (all patients, censored at 1 year)"),
  list(start = 365,  end = 1095, lab = "1-3 years (landmark at 1 year)"),
  list(start = 1095, end = Inf,  lab = ">3 years (landmark at 3 years)")
)
SPECS <- list(
  list(covs = character(0), lab = "unadjusted"),
  list(covs = COV_AGE,      lab = "age only"),
  list(covs = COV_PARS,     lab = "reduced (age, DM, CKD, AF)"),
  list(covs = COV_M4,       lab = "Model 4 (re-fitted within window)")
)

# 센터 층화 재표본 (센터 크기 보존; 창·센터 구성 유지)
boot_idx_center <- function(center) {
  unlist(lapply(split(seq_along(center), center),
                function(ix) if (length(ix) == 1) ix else sample(ix, length(ix), replace = TRUE)),
         use.names = FALSE)
}

###############################################################################
# 코호트
###############################################################################

co <- build_cohort()
stopifnot(nrow(co) == 6777)
an <- co[stats::complete.cases(co[, COV_M4, drop = FALSE]), , drop = FALSE]
logmsg("[cohort] full N = %d (F %d) | Model 4 complete case n = %d (F %d) | centres %d",
       nrow(co), sum(co$female), nrow(an), sum(an$female), length(unique(co$center_anon)))

###############################################################################
# A-3 (R2-M3) — 3성분 procedural intensity score
#   (i)   점수 분포 0/1/2/3 (성별, 빈도·%)
#   (ii)  보정 주변확률(표준화, g-computation): 성분별 · 점수 >=1 · >=2
#   (iii) 성분별 보정 pp 차 — 3성분 점수 차와 정확히 합산 (table15d 의 3성분판)
#   + 주 추정치 -0.088 (-0.140 to -0.036) 재현 확인
###############################################################################

message("\n================ A-3: 3-component procedural intensity score ================")

IC3 <- c("Bilati_inj", "retro_attempt_any", "IVUS")
IC3_LAB <- c("Bilateral contrast injection", "Retrograde approach attempted",
             "Intravascular ultrasound use")
an$pis3 <- ifelse(stats::complete.cases(an[, IC3, drop = FALSE]),
                  rowSums(an[, IC3, drop = FALSE]), NA_real_)
a3 <- droplevels(an[!is.na(an$pis3), , drop = FALSE])
logmsg("[A-3] 3-component score sample: n = %d (men %d / women %d) of the %d Model 4 complete cases",
       nrow(a3), sum(a3$female == 0), sum(a3$female == 1), nrow(an))
stopifnot(nrow(a3) == 5879, sum(a3$female == 1) == 1107)

SCORE_NOTE <- paste0(
  "Three-component procedural intensity score (bilateral contrast injection + retrograde approach ",
  "attempted + intravascular ultrasound use; 0-3, equally weighted); complete case for all three ",
  "components within the Model 4 complete-case sample. Covariates: Model 4 minus centre volume ",
  "(absorbed by the centre term). A higher score denotes a more intensive procedural strategy and ",
  "is not a measure of quality of care.")

a3_rows <- list()
addA3 <- function(x) a3_rows[[length(a3_rows) + 1]] <<- x

## --- (0) 주 추정치 재현 -----------------------------------------------------
m_ref <- lme4::lmer(stats::as.formula(paste0(
  "pis3 ~ female + ", paste(COV_M4_NOCENTER, collapse = " + "), " + (1 | center_anon)")),
  data = a3, REML = TRUE)
s <- summary(m_ref)
b <- s$coefficients["female", "Estimate"]; se <- s$coefficients["female", "Std. Error"]
p <- 2 * stats::pnorm(-abs(b / se))
logmsg("[A-3(0)] reproduction of the primary score estimate: %.4f (%.4f to %.4f), p = %s (published -0.088 [-0.140 to -0.036])",
       b, b - 1.96 * se, b + 1.96 * se, fmt_p(p))
stopifnot(abs(b - (-0.0877)) < 0.002)           # table15c (A16-1) 재현 assertion
addA3(rw("A-3(0)", "3-component procedural intensity score (0-3)",
         "Linear mixed model, centre random intercept (reproduction of the published primary estimate)",
         nrow(a3), sum(a3$female == 1), b, b - 1.96 * se, b + 1.96 * se,
         fmt3(b, b - 1.96 * se, b + 1.96 * se), fmt_p(p),
         paste0("Adjusted mean difference in score, women minus men. Reproduces table15c row A16-1 ",
                "(-0.088, -0.140 to -0.036) exactly; this row is a verification, not a new analysis. ",
                SCORE_NOTE),
         measure = "adjusted mean difference in score (women - men)"))

## --- (i) 분포 ---------------------------------------------------------------
for (g in c("Male", "Female")) {
  x <- a3$pis3[a3$sex_f == g]
  for (k in 0:3) {
    cnt <- sum(x == k); pct <- 100 * cnt / length(x)
    addA3(rw("A-3(i)", sprintf("Score = %d", k),
             sprintf("Crude distribution, %s", ifelse(g == "Male", "men", "women")),
             length(x), ifelse(g == "Female", length(x), 0L),
             pct, NA_real_, NA_real_, sprintf("%d (%.1f%%)", cnt, pct), NA_character_,
             SCORE_NOTE, measure = "percent of patients within sex",
             sex = g, score_level = k, count = cnt, pct = round(pct, 1)))
  }
  addA3(rw("A-3(i)", "Mean score (SD)",
           sprintf("Crude mean, %s", ifelse(g == "Male", "men", "women")),
           length(x), ifelse(g == "Female", length(x), 0L),
           mean(x), NA_real_, NA_real_, sprintf("%.2f (%.2f)", mean(x), stats::sd(x)),
           NA_character_, SCORE_NOTE, measure = "mean score", sex = g))
}
cd <- mean(a3$pis3[a3$female == 1]) - mean(a3$pis3[a3$female == 0])
tt <- stats::t.test(pis3 ~ female, data = a3)
addA3(rw("A-3(i)", "Crude mean difference (women - men)", "Welch t-test (crude, unadjusted)",
         nrow(a3), sum(a3$female == 1), cd, -rev(tt$conf.int)[1], -rev(tt$conf.int)[2],
         fmt3(cd, -tt$conf.int[2], -tt$conf.int[1]), fmt_p(tt$p.value),
         SCORE_NOTE, measure = "crude mean difference in score (women - men)"))

## --- (ii) 보정 주변확률 (표준화 g-computation, 센터 랜덤절편) ---------------
# 각 환자의 관측 공변량과 관측 센터를 유지한 채 성별만 0/1 로 바꿔 예측하고 평균한다
# (Model 4 공변량 분포 + 관측 센터 분포에 대한 표준화). CI 는 고정효과의 표집분포
# 에서 MC 추출(1,000회), 센터 랜덤효과는 조건부 최빈값에 고정.
std_marg_glmm <- function(d, yvar, outcome_lab, note_extra = "", boot_B = 0L) {
  dd <- droplevels(d[!is.na(d[[yvar]]), , drop = FALSE]); dd$.y <- dd[[yvar]]
  f <- stats::as.formula(paste0(".y ~ female + ", paste(COV_M4_NOCENTER, collapse = " + "),
                               " + (1 | center_anon)"))
  m <- lme4::glmer(f, data = dd, family = stats::binomial(), nAGQ = 0L,
                   control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
  X  <- as.matrix(lme4::getME(m, "X"))
  zb <- as.vector(stats::predict(m, re.form = NULL) - stats::predict(m, re.form = NA))
  bb <- lme4::fixef(m); V <- as.matrix(stats::vcov(m)); jf <- which(colnames(X) == "female")
  gfun <- function(bv, Xm, z) {
    X0 <- Xm; X0[, jf] <- 0; X1 <- Xm; X1[, jf] <- 1
    p0 <- mean(stats::plogis(as.vector(X0 %*% bv) + z))
    p1 <- mean(stats::plogis(as.vector(X1 %*% bv) + z))
    c(p0, p1, p1 - p0)
  }
  pt <- gfun(bb, X, zb)
  L <- chol(V)
  dr <- t(vapply(seq_len(B_MC), function(i) gfun(as.vector(bb + t(L) %*% stats::rnorm(length(bb))), X, zb),
                 numeric(3)))
  ci <- apply(dr, 2, stats::quantile, probs = c(.025, .975))
  sm <- summary(m)$coefficients
  orf <- exp(sm["female", "Estimate"]); orl <- exp(sm["female", "Estimate"] - 1.96 * sm["female", "Std. Error"])
  orh <- exp(sm["female", "Estimate"] + 1.96 * sm["female", "Std. Error"])
  pv  <- sm["female", "Pr(>|z|)"]
  nt <- paste0("Model-standardized (g-computation) probability from a logistic model with a centre ",
               "random intercept and Model 4 covariates (minus centre volume): each patient's own ",
               "covariates and own centre are retained and only sex is set to male / female, then ",
               "averaged over the analysis sample. 95% CI by Monte Carlo sampling from the ",
               "asymptotic distribution of the fixed effects (", B_MC, " draws), centre random ",
               "effects held at their conditional modes. The random-intercept model is fitted with ",
               "the fast approximation (nAGQ = 0), as in table14 row 2b; for bilateral contrast ",
               "injection the full Laplace approximation gives an odds ratio of 0.7646 versus ",
               "0.7655 here, a difference confined to the third decimal (0.76 vs 0.77 when rounded ",
               "to two). ", note_extra)
  out <- list(
    rw("A-3(ii)", outcome_lab, "Adjusted marginal probability, men (standardized)",
       nrow(dd), sum(dd$female == 1), 100 * pt[1], 100 * ci[1, 1], 100 * ci[2, 1],
       sprintf("%.1f%% (%.1f-%.1f)", 100 * pt[1], 100 * ci[1, 1], 100 * ci[2, 1]),
       NA_character_, nt, measure = "adjusted marginal probability (%)"),
    rw("A-3(ii)", outcome_lab, "Adjusted marginal probability, women (standardized)",
       nrow(dd), sum(dd$female == 1), 100 * pt[2], 100 * ci[1, 2], 100 * ci[2, 2],
       sprintf("%.1f%% (%.1f-%.1f)", 100 * pt[2], 100 * ci[1, 2], 100 * ci[2, 2]),
       NA_character_, nt, measure = "adjusted marginal probability (%)"),
    rw("A-3(ii)", outcome_lab, "Adjusted difference in marginal probability (women - men)",
       nrow(dd), sum(dd$female == 1), 100 * pt[3], 100 * ci[1, 3], 100 * ci[2, 3],
       fmtpp(100 * pt[3], 100 * ci[1, 3], 100 * ci[2, 3]), fmt_p(pv),
       paste0(nt, " P value is the Wald test for the sex term of the same model; the adjusted odds ",
              "ratio is ", fmt2(orf, orl, orh), "."),
       measure = "adjusted risk difference (percentage points, women - men)",
       adj_OR_fmt = fmt2(orf, orl, orh)))
  # 선택적 센터 층화 부트스트랩 (주요 지표 확인용)
  if (boot_B > 0) {
    set.seed(SEED)
    bs <- matrix(NA_real_, boot_B, 3)
    for (i in seq_len(boot_B)) {
      ix <- boot_idx_center(dd$center_anon)
      db <- droplevels(dd[ix, , drop = FALSE])
      mb <- tryCatch(lme4::glmer(f, data = db, family = stats::binomial(), nAGQ = 0L,
                                 control = lme4::glmerControl(optimizer = "bobyqa",
                                                              optCtrl = list(maxfun = 3e5))),
                     error = function(e) NULL, warning = function(w) NULL)
      if (is.null(mb)) next
      Xb <- as.matrix(lme4::getME(mb, "X"))
      zbb <- as.vector(stats::predict(mb, re.form = NULL) - stats::predict(mb, re.form = NA))
      if (!("female" %in% colnames(Xb))) next
      jfb <- which(colnames(Xb) == "female")
      X0 <- Xb; X0[, jfb] <- 0; X1 <- Xb; X1[, jfb] <- 1
      p0 <- mean(stats::plogis(as.vector(X0 %*% lme4::fixef(mb)) + zbb))
      p1 <- mean(stats::plogis(as.vector(X1 %*% lme4::fixef(mb)) + zbb))
      bs[i, ] <- c(p0, p1, p1 - p0)
    }
    q <- stats::quantile(bs[, 3], c(.025, .975), na.rm = TRUE)
    out[[length(out) + 1]] <- rw(
      "A-3(ii)", outcome_lab,
      "Adjusted difference in marginal probability (women - men), centre-stratified bootstrap CI",
      nrow(dd), sum(dd$female == 1), 100 * pt[3], 100 * q[1], 100 * q[2],
      fmtpp(100 * pt[3], 100 * q[1], 100 * q[2]), NA_character_,
      sprintf(paste0("Same point estimate; 95%% CI from %d centre-stratified bootstrap replications ",
                     "(patients resampled within centre, centre sizes preserved; %d replications ",
                     "converged). Sensitivity check on the Monte Carlo interval."),
              boot_B, sum(!is.na(bs[, 3]))),
      measure = "adjusted risk difference (percentage points, women - men)")
  }
  dplyr::bind_rows(out)
}

a3$pis3_ge1 <- as.numeric(a3$pis3 >= 1)
a3$pis3_ge2 <- as.numeric(a3$pis3 >= 2)
TARG <- list(list(v = "Bilati_inj", lab = IC3_LAB[1], bb = B_BOOT),
             list(v = "retro_attempt_any", lab = IC3_LAB[2], bb = 0L),
             list(v = "IVUS", lab = IC3_LAB[3], bb = 0L),
             list(v = "pis3_ge1", lab = "Score >=1 (any of the three)", bb = B_BOOT),
             list(v = "pis3_ge2", lab = "Score >=2 (at least two of the three)", bb = 0L))
for (z in TARG) {
  message(sprintf("  [A-3(ii)] %s ...", z$lab))
  addA3(std_marg_glmm(a3, z$v, z$lab, SCORE_NOTE, boot_B = z$bb))
}

## --- (iii) 성분별 보정 pp 차 (가법 분해, 항등식 검증) -----------------------
# 점수는 성분의 합이므로, 동일 표본·동일 설계행렬의 선형모형에서 성분별 female
# 계수의 합은 총 점수의 female 계수와 정확히 일치한다. table15d(4성분)의 3성분판.
dec3 <- local({
  d <- a3
  rhs <- paste(c("female", COV_M4_NOCENTER, "factor(center_anon)"), collapse = " + ")
  fit1 <- function(y) {
    d$.y <- y
    m <- stats::lm(stats::as.formula(paste0(".y ~ ", rhs)), data = d)
    ct <- lmtest::coeftest(m, vcov. = sandwich::vcovCL(m, cluster = d$center_anon))
    c(b = unname(ct["female", "Estimate"]), se = unname(ct["female", "Std. Error"]),
      p = unname(ct["female", "Pr(>|t|)"]))
  }
  tot <- fit1(d$pis3)
  parts <- t(vapply(IC3, function(v) fit1(d[[v]]), numeric(3)))
  stopifnot(abs(sum(parts[, "b"]) - tot["b"]) < 1e-8)
  logmsg("[A-3(iii)] additive identity verified: components sum to %.4f, total model %.4f",
         sum(parts[, "b"]), tot["b"])
  M <- rbind(parts, TOTAL = tot)
  labs <- c(IC3_LAB, "TOTAL (3-component score)")
  dplyr::bind_rows(lapply(seq_len(nrow(M)), function(i) {
    e <- 100 * M[i, "b"]; l <- 100 * (M[i, "b"] - 1.96 * M[i, "se"]); h <- 100 * (M[i, "b"] + 1.96 * M[i, "se"])
    rw("A-3(iii)", labs[i],
       "Linear probability model, Model 4 covariates + centre fixed effects, cluster-robust SE",
       nrow(d), sum(d$female == 1), e, l, h, fmtpp(e, l, h), fmt_p(M[i, "p"]),
       paste0("Adjusted percentage-point difference, women minus men, within centre. Because the ",
              "score is the sum of its components and all models share one design matrix, the ",
              "component coefficients sum exactly to the coefficient for the total score; the ",
              "identity is asserted in code. Shares are descriptive. ", SCORE_NOTE),
       measure = "adjusted percentage-point difference (women - men)",
       crude_pct_men = round(100 * mean(if (i <= 3) d[[IC3[i]]][d$female == 0] else d$pis3[d$female == 0]), 1),
       crude_pct_women = round(100 * mean(if (i <= 3) d[[IC3[i]]][d$female == 1] else d$pis3[d$female == 1]), 1),
       share_of_total_pct = round(100 * M[i, "b"] / M[nrow(M), "b"], 1))
  }))
})
addA3(dec3)

A3 <- dplyr::bind_rows(a3_rows)
write.csv(A3, file.path(DIR_TAB, "table28_A3_intensity_score_3component.csv"),
          row.names = FALSE, na = "")
logmsg("[A-3] wrote %d rows to table28_A3_intensity_score_3component.csv", nrow(A3))

###############################################################################
# A-6 (R2-M6) — 성분별 성별 누적발생률 (사망 경쟁) + 재관류 Fine-Gray
#   전사망: Kaplan-Meier (경쟁사건 없음) / MI · 뇌졸중 · 모든 재관류:
#   Aalen-Johansen CIF (사망을 경쟁사건으로) + Gray 검정
#   재관류 Fine-Gray 하위분포 HR: 미보정 / 절약형 / Model 4
#   ⚠️ 등록부는 **첫 사건만** 기록하므로 재발 사건은 분석 불가.
###############################################################################

message("\n================ A-6: component cumulative incidence with death competing ================")

TP <- c(365, 1095, 1825); TPL <- c("1-year", "3-year", "5-year")
FIRST_ONLY <- paste0(
  "The registry records the date of the first occurrence of each event type only; recurrent events ",
  "after a first event are not recorded and therefore cannot be analysed. ")

a6_rows <- list()
addA6 <- function(x) a6_rows[[length(a6_rows) + 1]] <<- x

## --- 전사망: KM -------------------------------------------------------------
dkm <- make_surv(co, "Death")
fkm <- survival::survfit(Surv(.time, .event) ~ sex_f, data = dkm, conf.type = "log-log")
skm <- summary(fkm, times = TP, extend = TRUE)
lr  <- survival::survdiff(Surv(.time, .event) ~ sex_f, data = dkm)
p_lr <- stats::pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
for (g in c("Male", "Female")) for (j in seq_along(TP)) {
  i <- which(as.character(skm$strata) == paste0("sex_f=", g) & skm$time == TP[j])
  addA6(rw("A-6", "All-cause death", sprintf("Kaplan-Meier, %s, %s",
           ifelse(g == "Male", "men", "women"), TPL[j]),
           sum(dkm$sex_f == g), sum(dkm$sex_f == "Female"),
           100 * (1 - skm$surv[i]), 100 * (1 - skm$upper[i]), 100 * (1 - skm$lower[i]),
           sprintf("%.1f%% (%.1f-%.1f)", 100 * (1 - skm$surv[i]), 100 * (1 - skm$upper[i]),
                   100 * (1 - skm$lower[i])), fmt_p(p_lr),
           paste0("Kaplan-Meier failure estimate with log-log 95% CI; no competing event applies to ",
                  "all-cause death. P value is the log-rank test comparing women with men over the ",
                  "whole follow-up. Full cohort; day-0 events placed at 0.5 days. ", FIRST_ONLY),
           measure = "cumulative incidence (%)", sex = g, timepoint = TPL[j],
           n_events = sum(dkm$.event[dkm$sex_f == g] == 1), n_at_risk_start = sum(dkm$sex_f == g)))
}

## --- 비치명 성분: Aalen-Johansen CIF (사망 경쟁) ----------------------------
comp_cif <- function(dat, ev, tm, lab) {
  d <- dat
  t <- d[[tm]]; e <- ifelse(d[[ev]] == 1, 1L, ifelse(d$Death == 1, 2L, 0L))
  bad <- d[[ev]] == 1 & d$Death == 1 & t > d$death_fudur
  nbad <- sum(bad, na.rm = TRUE)
  if (nbad > 0) { t[bad] <- d$death_fudur[bad]; e[bad] <- 2L }   # 사망 후 기록된 사건은 경쟁사망으로 귀속
  t[e > 0 & t == 0] <- 0.5
  ok <- !is.na(t) & !is.na(e)
  d <- d[ok, , drop = FALSE]; t <- t[ok]; e <- e[ok]
  ci <- cmprsk::cuminc(ftime = t, fstatus = e, group = as.character(d$sex_f), cencode = 0)
  est <- cmprsk::timepoints(ci, times = TP)
  gp <- unname(ci$Tests["1", "pv"])
  nt <- paste0("Aalen-Johansen cumulative incidence with all-cause death as the competing event; ",
               "95% CI on the log scale from the variance returned by cmprsk::timepoints. P value is ",
               "Gray's test for the equality of the cause-specific cumulative incidence functions ",
               "between women and men. Full cohort; day-0 events placed at 0.5 days. ",
               if (nbad > 0) sprintf(paste0("%d patient(s) with a recorded component event dated after ",
                                            "the recorded death date were reclassified to the competing ",
                                            "death. "), nbad) else "",
               FIRST_ONLY)
  out <- list()
  for (g in c("Male", "Female")) for (j in seq_along(TP)) {
    rn <- paste(g, 1)
    F <- est$est[rn, j]; v <- est$var[rn, j]
    se <- sqrt(v); lo <- F * exp(-1.96 * se / F); hi <- F * exp(1.96 * se / F)
    out[[length(out) + 1]] <- rw(
      "A-6", lab, sprintf("Aalen-Johansen CIF, %s, %s", ifelse(g == "Male", "men", "women"), TPL[j]),
      sum(d$sex_f == g), sum(d$sex_f == "Female"),
      100 * F, 100 * lo, 100 * hi, sprintf("%.2f%% (%.2f-%.2f)", 100 * F, 100 * lo, 100 * hi),
      fmt_p(gp), nt, measure = "cumulative incidence (%)", sex = g, timepoint = TPL[j],
      n_events = sum(e[d$sex_f == g] == 1), n_competing_deaths = sum(e[d$sex_f == g] == 2))
  }
  dplyr::bind_rows(out)
}
addA6(comp_cif(co, "MI", "MI_fudur", "Myocardial infarction"))
addA6(comp_cif(co, "CVA_re", "CVA_fudur", "Stroke"))
addA6(comp_cif(co, "AnyRevasc", "AnyRevasc_fudur",
               "Any coronary revascularization (TVR, non-target-lesion PCI, or CABG)"))

## --- 재관류 Fine-Gray 하위분포 HR ------------------------------------------
fg_revasc <- function(dat, covs, lab_model) {
  d <- dat
  t <- d$AnyRevasc_fudur
  e <- ifelse(d$AnyRevasc == 1, 1L, ifelse(d$Death == 1, 2L, 0L))
  bad <- d$AnyRevasc == 1 & d$Death == 1 & t > d$death_fudur
  if (any(bad)) { t[bad] <- d$death_fudur[bad]; e[bad] <- 2L }
  t[e > 0 & t == 0] <- 0.5
  keep <- !is.na(t) & !is.na(e)
  if (length(covs)) keep <- keep & stats::complete.cases(d[, covs, drop = FALSE])
  d <- droplevels(d[keep, , drop = FALSE]); t <- t[keep]; e <- e[keep]
  X <- if (length(covs))
    stats::model.matrix(stats::as.formula(paste0("~ female + ", paste(covs, collapse = " + "))),
                        data = d)[, -1, drop = FALSE]
  else matrix(d$female, ncol = 1, dimnames = list(NULL, "female"))
  fg <- cmprsk::crr(ftime = t, fstatus = e, cov1 = X, failcode = 1, cencode = 0)
  s <- summary(fg); i <- which(rownames(s$coef) == "female")
  rw("A-6", "Any coronary revascularization (TVR, non-target-lesion PCI, or CABG)",
     paste0("Fine-Gray subdistribution hazard, ", lab_model),
     nrow(d), sum(d$female == 1), unname(s$conf.int[i, 1]), unname(s$conf.int[i, 3]),
     unname(s$conf.int[i, 4]), fmt2(s$conf.int[i, 1], s$conf.int[i, 3], s$conf.int[i, 4]),
     fmt_p(unname(s$coef[i, 5])),
     paste0("Subdistribution hazard ratio for women versus men, with all-cause death as the ",
            "competing event. ", FIRST_ONLY,
            "Fine-Gray models do not carry the centre-clustered robust variance used for the Cox ",
            "models; standard errors are model based."),
     measure = "subdistribution hazard ratio (women vs men)",
     n_events = sum(e == 1), n_competing_deaths = sum(e == 2))
}
addA6(fg_revasc(co, character(0), "unadjusted (full cohort)"))
addA6(fg_revasc(an, character(0), "unadjusted (Model 4 complete-case sample)"))
addA6(fg_revasc(an, COV_PARS, "reduced (age, DM, CKD, AF)"))
addA6(fg_revasc(an, COV_M4, "Model 4"))

A6 <- dplyr::bind_rows(a6_rows)
write.csv(A6, file.path(DIR_TAB, "table29_A6_component_cuminc_competing.csv"),
          row.names = FALSE, na = "")
logmsg("[A-6] wrote %d rows to table29_A6_component_cuminc_competing.csv", nrow(A6))

###############################################################################
# A-9 (R2-m4) — 모형 표준화(g-formula) 보정 절대위험 + 보정 주변확률
#   (i)  MACCE / hard composite(DMS) 의 1·3·5년 보정 절대위험(성별)과 위험차
#        Model 4 · 절약형, 전기간 모형 + 창별(time-varying) female 항 모형
#        CI: 센터 층화 부트스트랩 500회
#   (ii) bilateral injection · technical success 의 보정 주변확률(성별) + pp 차
###############################################################################

message("\n================ A-9: model-standardized absolute risks ================")

a9_rows <- list()
addA9 <- function(x) a9_rows[[length(a9_rows) + 1]] <<- x

# --- 전기간 Cox 에서의 g-formula 위험 (Breslow 기저위험) ---------------------
risk_all <- function(d, covs, times) {
  f <- stats::as.formula(paste0("Surv(.time, .event) ~ ",
                                paste(c("female", covs), collapse = " + ")))
  m <- survival::coxph(f, data = d, ties = "breslow")
  b <- stats::coef(m); b[is.na(b)] <- 0
  bh <- survival::basehaz(m, centered = FALSE)
  X <- stats::model.matrix(f, data = d)[, -1, drop = FALSE]
  X <- X[, names(b), drop = FALSE]
  jf <- which(colnames(X) == "female")
  X0 <- X; X0[, jf] <- 0; X1 <- X; X1[, jf] <- 1
  lp0 <- as.vector(X0 %*% b); lp1 <- as.vector(X1 %*% b)
  H0 <- function(tt) { i <- findInterval(tt, bh$time); ifelse(i == 0, 0, bh$hazard[pmax(i, 1)]) }
  vapply(times, function(tt) {
    h <- H0(tt)
    c(men = mean(1 - exp(-h * exp(lp0))), women = mean(1 - exp(-h * exp(lp1))))
  }, numeric(2))
}

# --- 창별 female 항(piecewise) Cox 에서의 g-formula 위험 --------------------
# PH 가정을 전기간에 부과하지 않는다: female 효과를 0-1 / 1-3 / >3년 구간별로
# 따로 추정하고(공변량 효과는 공유), 기저위험 1개를 구간별로 적분한다.
risk_pw <- function(d, covs, times, cuts = c(365, 1095)) {
  # survSplit 은 관측시간이 모두 0 보다 커야 한다. person-time 을 기여하지 않는
  # 환자(추적 0일·무사건, 대부분 한 센터의 추적 미제출)는 표준화 모집단에 남기기
  # 위해 0.25일로 절단한다 — 사건이 없으므로 부분우도·기저위험 기여는 없다
  # (전기간 모형과 동일한 표준화 모집단을 유지하기 위한 처리).
  d$.time[d$.event == 0 & d$.time <= 0] <- 0.25
  ds <- survival::survSplit(Surv(.time, .event) ~ ., data = d, cut = cuts,
                            episode = "ep", start = "tstart")
  ds$fem_w1 <- ds$female * (ds$ep == 1)
  ds$fem_w2 <- ds$female * (ds$ep == 2)
  ds$fem_w3 <- ds$female * (ds$ep == 3)
  rhs <- paste(c("fem_w1", "fem_w2", "fem_w3", covs), collapse = " + ")
  m <- survival::coxph(stats::as.formula(paste0("Surv(tstart, .time, .event) ~ ", rhs)),
                       data = ds, ties = "breslow")
  b <- stats::coef(m); b[is.na(b)] <- 0
  bf <- b[c("fem_w1", "fem_w2", "fem_w3")]
  bc <- b[setdiff(names(b), c("fem_w1", "fem_w2", "fem_w3"))]
  bh <- survival::basehaz(m, centered = FALSE)
  Xc <- if (length(bc)) {
    MM <- stats::model.matrix(stats::as.formula(paste0("~ ", paste(covs, collapse = " + "))),
                              data = d)[, -1, drop = FALSE]
    MM[, names(bc), drop = FALSE]
  } else matrix(0, nrow(d), 0)
  lpc <- if (length(bc)) as.vector(Xc %*% bc) else rep(0, nrow(d))
  H0 <- function(tt) { i <- findInterval(tt, bh$time); ifelse(i == 0, 0, bh$hazard[pmax(i, 1)]) }
  bnds <- c(0, cuts, Inf)
  vapply(times, function(tt) {
    dH <- vapply(seq_len(3), function(k)
      max(0, H0(min(tt, bnds[k + 1])) - H0(min(tt, bnds[k]))), numeric(1))
    cum <- function(s) {
      H <- rowSums(vapply(seq_len(3), function(k) dH[k] * exp(lpc + bf[k] * s), numeric(nrow(d))))
      mean(1 - exp(-H))
    }
    c(men = cum(0), women = cum(1))
  }, numeric(2))
}

gcomp <- function(dat, key, covs, cov_lab, times = TP, pw = FALSE, B = B_BOOT) {
  d <- make_surv(dat, key)
  d <- cc_subset(d, covs)
  fun <- if (pw) risk_pw else risk_all
  pt <- fun(d, covs, times)
  set.seed(SEED)
  bs <- array(NA_real_, c(B, 2, length(times)))
  for (i in seq_len(B)) {
    ix <- boot_idx_center(d$center_anon)
    db <- droplevels(d[ix, , drop = FALSE])
    r <- tryCatch(fun(db, covs, times), error = function(e) NULL)
    if (!is.null(r)) bs[i, , ] <- r
  }
  nconv <- sum(!is.na(bs[, 1, 1]))
  mlab <- paste0(cov_lab, if (pw) ", female term estimated separately in each window" else "")
  nt <- paste0("Model-standardized (g-formula) absolute risk: the fitted Cox model and its Breslow ",
               "baseline hazard are applied to every patient in the analysis sample twice, once as ",
               "male and once as female, holding all other covariates at their observed values, and ",
               "the predicted risks are averaged. 95% CI from ", nconv, " of ", B,
               " centre-stratified bootstrap replications (patients resampled within centre). ",
               if (pw) paste0("The sex term is estimated separately within 0-1, 1-3 and >3 years ",
                              "(boundaries 365 and 1095 days) so that no proportional-hazards ",
                              "assumption is imposed on the sex effect across the whole follow-up. ")
               else paste0("This specification imposes proportional hazards for sex over the whole ",
                           "follow-up. That assumption is not satisfied for MACCE: the ",
                           "scaled-Schoenfeld test for the sex term over the whole follow-up is ",
                           "P = 0.014 (table8a), while it is satisfied for the hard composite ",
                           "(P = 0.108). The window-specific specification is reported alongside ",
                           "for this reason. "),
               "Day-0 events placed at 0.5 days.")
  out <- list()
  for (j in seq_along(times)) {
    q <- function(v) stats::quantile(v, c(.025, .975), na.rm = TRUE)
    cm <- q(bs[, 1, j]); cw <- q(bs[, 2, j]); cd <- q(bs[, 2, j] - bs[, 1, j])
    dd <- pt[2, j] - pt[1, j]
    out[[length(out) + 1]] <- rw("A-9(i)", OUT[[key]]$lab,
      sprintf("Standardized risk, men, %s (%s)", TPL[j], mlab),
      nrow(d), sum(d$female == 1), 100 * pt[1, j], 100 * cm[1], 100 * cm[2],
      sprintf("%.1f%% (%.1f-%.1f)", 100 * pt[1, j], 100 * cm[1], 100 * cm[2]), NA_character_,
      nt, measure = "standardized absolute risk (%)", timepoint = TPL[j])
    out[[length(out) + 1]] <- rw("A-9(i)", OUT[[key]]$lab,
      sprintf("Standardized risk, women, %s (%s)", TPL[j], mlab),
      nrow(d), sum(d$female == 1), 100 * pt[2, j], 100 * cw[1], 100 * cw[2],
      sprintf("%.1f%% (%.1f-%.1f)", 100 * pt[2, j], 100 * cw[1], 100 * cw[2]), NA_character_,
      nt, measure = "standardized absolute risk (%)", timepoint = TPL[j])
    out[[length(out) + 1]] <- rw("A-9(i)", OUT[[key]]$lab,
      sprintf("Standardized risk difference (women - men), %s (%s)", TPL[j], mlab),
      nrow(d), sum(d$female == 1), 100 * dd, 100 * cd[1], 100 * cd[2],
      fmtpp(100 * dd, 100 * cd[1], 100 * cd[2]), NA_character_,
      nt, measure = "standardized risk difference (percentage points, women - men)",
      timepoint = TPL[j])
  }
  dplyr::bind_rows(out)
}

for (key in c("MACCE_pub", "DMS")) {
  for (sp in list(list(cv = COV_M4, lb = "Model 4"),
                  list(cv = COV_PARS, lb = "reduced (age, DM, CKD, AF)"))) {
    message(sprintf("  [A-9(i)] %s / %s / all-period ...", key, sp$lb))
    addA9(gcomp(an, key, sp$cv, sp$lb, pw = FALSE))
    message(sprintf("  [A-9(i)] %s / %s / window-specific sex term ...", key, sp$lb))
    addA9(gcomp(an, key, sp$cv, sp$lb, pw = TRUE))
  }
}

## --- (ii) bilateral injection · technical success 의 보정 주변확률 ----------
std_marg_glm <- function(d, yvar, outcome_lab) {
  dd <- droplevels(d[!is.na(d[[yvar]]), , drop = FALSE]); dd$.y <- dd[[yvar]]
  f <- stats::as.formula(paste0(".y ~ female + ", paste(COV_M4, collapse = " + ")))
  m <- stats::glm(f, data = dd, family = stats::binomial())
  V <- sandwich::vcovCL(m, cluster = dd$center_anon)
  b <- stats::coef(m); keep <- !is.na(b); b <- b[keep]
  X <- stats::model.matrix(m)[, keep, drop = FALSE]
  V <- V[keep, keep, drop = FALSE]
  jf <- which(colnames(X) == "female")
  X0 <- X; X0[, jf] <- 0; X1 <- X; X1[, jf] <- 1
  e0 <- as.vector(stats::plogis(X0 %*% b)); e1 <- as.vector(stats::plogis(X1 %*% b))
  p0 <- mean(e0); p1 <- mean(e1); n <- nrow(X)
  g0 <- as.vector(crossprod(X0, e0 * (1 - e0)) / n)
  g1 <- as.vector(crossprod(X1, e1 * (1 - e1)) / n)
  sed <- function(g) sqrt(as.numeric(t(g) %*% V %*% g))
  ct <- lmtest::coeftest(m, vcov. = sandwich::vcovCL(m, cluster = dd$center_anon))
  orf <- exp(ct["female", "Estimate"])
  orl <- exp(ct["female", "Estimate"] - 1.96 * ct["female", "Std. Error"])
  orh <- exp(ct["female", "Estimate"] + 1.96 * ct["female", "Std. Error"])
  nt <- paste0("Model-standardized (g-computation) probability from the adjusted logistic model with ",
               "Model 4 covariates and centre-clustered robust standard errors: each patient's ",
               "covariates are held at their observed values and only sex is set to male / female, ",
               "then the predicted probabilities are averaged. 95% CI by the delta method on the ",
               "centre-clustered robust covariance. Adjusted odds ratio for women versus men ",
               fmt2(orf, orl, orh), ".")
  mk <- function(est, se, lab, meas) rw("A-9(ii)", outcome_lab, lab, nrow(dd), sum(dd$female == 1),
    100 * est, 100 * (est - 1.96 * se), 100 * (est + 1.96 * se),
    if (meas == "d") fmtpp(100 * est, 100 * (est - 1.96 * se), 100 * (est + 1.96 * se))
    else sprintf("%.1f%% (%.1f-%.1f)", 100 * est, 100 * (est - 1.96 * se), 100 * (est + 1.96 * se)),
    if (meas == "d") fmt_p(ct["female", "Pr(>|z|)"]) else NA_character_, nt,
    measure = if (meas == "d") "adjusted risk difference (percentage points, women - men)"
              else "adjusted marginal probability (%)",
    adj_OR_fmt = fmt2(orf, orl, orh))
  dplyr::bind_rows(
    mk(p0, sed(g0), "Adjusted marginal probability, men (standardized, Model 4 logistic)", "p"),
    mk(p1, sed(g1), "Adjusted marginal probability, women (standardized, Model 4 logistic)", "p"),
    mk(p1 - p0, sed(g1 - g0),
       "Adjusted difference in marginal probability (women - men), Model 4 logistic", "d"))
}
for (z in list(list(v = "Bilati_inj", lab = "Bilateral contrast injection"),
               list(v = "cto_success", lab = "Technical success of CTO PCI"))) {
  message(sprintf("  [A-9(ii)] %s ...", z$lab))
  addA9(std_marg_glm(an, z$v, z$lab))
  x <- std_marg_glmm(an, z$v, z$lab,
    paste0("Centre random-intercept specification of the same quantity, reported alongside the ",
           "Model 4 logistic model because the manuscript's headline estimate for bilateral ",
           "injection is the within-centre (random-intercept) one."))
  x$analysis <- "A-9(ii)"
  addA9(x)
}

A9 <- dplyr::bind_rows(a9_rows)
write.csv(A9, file.path(DIR_TAB, "table30_A9_standardized_absolute_risk.csv"),
          row.names = FALSE, na = "")
logmsg("[A-9] wrote %d rows to table30_A9_standardized_absolute_risk.csv", nrow(A9))

###############################################################################
# A-8 (R1-W1) — sex x era 상호작용
#   era(주) = index_half_f (2014-2018 vs 2019-2023)
#   era(민감도) = index_period_f (2년 구간 5범주; 기존 Model 4 의 era 변수)
#   지표: 3성분 점수(선형혼합, 센터 랜덤절편) · bilateral injection(로지스틱 GLMM)
###############################################################################

message("\n================ A-8: sex by era interaction ================")

a8_rows <- list()
addA8 <- function(x) a8_rows[[length(a8_rows) + 1]] <<- x

COV_H <- c(setdiff(COV_M4_NOCENTER, "index_period_f"), "index_half_f")   # era = 2구간
COV_P <- COV_M4_NOCENTER                                                  # era = 5구간(원 변수)
ERA_NOTE <- paste0(
  "Treatment era is already a Model 4 covariate (index treatment-era period). Here it is also ",
  "entered as an effect modifier. For the two-period analysis the five-category era covariate is ",
  "replaced by the dichotomous era variable so that the modifier and the covariate are the same ",
  "variable; for the five-category analysis the original Model 4 era variable is used. ",
  "Covariates otherwise as in Model 4 minus centre volume (absorbed by the centre term).")

era_models <- function(d, yvar, outcome_lab, eravar, covs, type = c("lmm", "glmm")) {
  type <- match.arg(type)
  dd <- droplevels(d[!is.na(d[[yvar]]) & !is.na(d[[eravar]]), , drop = FALSE]); dd$.y <- dd[[yvar]]
  covs_noera <- setdiff(covs, eravar)
  rhs_main <- paste(c("female", covs), collapse = " + ")
  rhs_int  <- paste(c(sprintf("female * %s", eravar), covs_noera), collapse = " + ")
  out <- list()
  fitter <- function(form, data, ml = FALSE) {
    if (type == "lmm")
      lme4::lmer(form, data = data, REML = !ml)
    else
      lme4::glmer(form, data = data, family = stats::binomial(), nAGQ = 0L,
                  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
  }
  f0 <- stats::as.formula(paste0(".y ~ ", rhs_main, " + (1 | center_anon)"))
  f1 <- stats::as.formula(paste0(".y ~ ", rhs_int,  " + (1 | center_anon)"))
  m0 <- fitter(f0, dd, ml = TRUE); m1 <- fitter(f1, dd, ml = TRUE)
  lrt <- stats::anova(m0, m1)
  p_int <- lrt$`Pr(>Chisq)`[2]; df_int <- lrt$Df[2]
  # era 층별 재적합
  for (lv in levels(dd[[eravar]])) {
    ds <- droplevels(dd[dd[[eravar]] == lv, , drop = FALSE])
    ms <- tryCatch(fitter(stats::as.formula(paste0(
      ".y ~ ", paste(c("female", covs_noera), collapse = " + "), " + (1 | center_anon)")), ds),
      error = function(e) NULL)
    if (is.null(ms)) next
    sm <- summary(ms)$coefficients
    b <- sm["female", "Estimate"]; se <- sm["female", "Std. Error"]
    p <- 2 * stats::pnorm(-abs(b / se))
    if (type == "lmm") {
      e <- b; l <- b - 1.96 * se; h <- b + 1.96 * se; ef <- fmt3(e, l, h)
      meas <- "adjusted mean difference in score (women - men)"
    } else {
      e <- exp(b); l <- exp(b - 1.96 * se); h <- exp(b + 1.96 * se); ef <- fmt2(e, l, h)
      meas <- "adjusted odds ratio (women vs men)"
    }
    out[[length(out) + 1]] <- rw("A-8", outcome_lab,
      sprintf("Era-specific estimate, %s (re-fitted within era)", lv),
      nrow(ds), sum(ds$female == 1), e, l, h, ef, fmt_p(p),
      paste0("Stratum-specific estimate; the formal test of effect modification is the ",
             "interaction row for this outcome (likelihood-ratio test). ", ERA_NOTE),
      measure = meas, era_variable = eravar, era_level = lv)
  }
  out[[length(out) + 1]] <- rw("A-8", outcome_lab,
    sprintf("Sex x era interaction, %s (%d df likelihood-ratio test)",
            ifelse(eravar == "index_half_f", "two eras (2014-2018 vs 2019-2023)",
                   "five 2-year eras"), df_int),
    nrow(dd), sum(dd$female == 1), NA_real_, NA_real_, NA_real_, NA_character_, fmt_p(p_int),
    paste0("Likelihood-ratio test of the sex by era interaction (models refitted by maximum ",
           "likelihood). A non-significant test means the data do not show that the sex ",
           "difference changed across the enrolment period; it does not establish that it was ",
           "constant. ", ERA_NOTE),
    measure = "interaction P value", era_variable = eravar, interaction_df = df_int,
    interaction_p = p_int)
  dplyr::bind_rows(out)
}

for (ev in c("index_half_f", "index_period_f")) {
  cvs <- if (ev == "index_half_f") COV_H else COV_P
  message(sprintf("  [A-8] 3-component score x %s ...", ev))
  addA8(era_models(a3, "pis3", "3-component procedural intensity score (0-3)", ev, cvs, "lmm"))
  message(sprintf("  [A-8] bilateral injection x %s ...", ev))
  addA8(era_models(an, "Bilati_inj", "Bilateral contrast injection", ev, cvs, "glmm"))
}

# --- era 별 성별 조율(crude rate) --------------------------------------------
for (ev in c("index_half_f", "index_period_f")) {
  for (i in seq_along(IC3)) {
    v <- IC3[i]
    for (lv in levels(co[[ev]])) {
      for (g in c("Male", "Female")) {
        x <- co[[v]][co[[ev]] == lv & co$sex_f == g]
        x <- x[!is.na(x)]
        addA8(rw("A-8(crude)", IC3_LAB[i],
          sprintf("Crude rate, %s, %s", lv, ifelse(g == "Male", "men", "women")),
          length(x), ifelse(g == "Female", length(x), 0L),
          100 * mean(x), NA_real_, NA_real_,
          sprintf("%d/%d (%.1f%%)", sum(x == 1), length(x), 100 * mean(x)), NA_character_,
          paste0("Crude (unadjusted) proportion within era and sex, full cohort; denominators ",
                 "exclude patients with the component not recorded. ", ERA_NOTE),
          measure = "crude proportion (%)", era_variable = ev, era_level = lv, sex = g))
      }
    }
  }
}

A8 <- dplyr::bind_rows(a8_rows)
write.csv(A8, file.path(DIR_TAB, "table31_A8_sex_by_era_interaction.csv"),
          row.names = FALSE, na = "")
logmsg("[A-8] wrote %d rows to table31_A8_sex_by_era_interaction.csv", nrow(A8))

###############################################################################
# A-7 (R2-M7) — 창별 모형 진단 정리표
#   창별(MACCE, Model 4 / 절약형): Schoenfeld global · female 항 p, 모수 개수,
#   사건/모수(EPP) + Model 4 VIF 요약(tableS8b) + 연령 선형성 LR 검정
###############################################################################

message("\n================ A-7: model diagnostics by landmark window ================")

a7_rows <- list()
addA7 <- function(x) a7_rows[[length(a7_rows) + 1]] <<- x

for (sp in list(list(cv = COV_M4, lb = "Model 4"),
                list(cv = COV_PARS, lb = "reduced (age, DM, CKD, AF)"))) {
  for (w in WINDOWS) {
    d <- cc_subset(window_data(an, "MACCE_pub", w$start, w$end), sp$cv)
    m <- cox_model(d, sp$cv)
    z <- survival::cox.zph(m)
    npar <- length(stats::coef(m)); ev <- m$nevent
    pf <- unname(z$table["female", "p"]); pg <- unname(z$table["GLOBAL", "p"])
    addA7(rw("A-7", OUT$MACCE_pub$lab,
      sprintf("%s, %s", w$lab, sp$lb), m$n, sum(d$female == 1),
      ev / npar, NA_real_, NA_real_, sprintf("%.2f", ev / npar), fmt_p(pf),
      paste0("Diagnostics for the window-specific Cox model (centre-clustered robust standard ",
             "errors). The reported P value is the scaled-Schoenfeld test for the sex term; the ",
             "global test is given separately. Events per parameter below 10 indicates a sparsely ",
             "estimated model, which is why the reduced model is the primary specification in the ",
             "late window."),
      measure = "events per parameter", window = w$lab, covariate_set = sp$lb,
      n_events = ev, n_parameters = npar, events_per_parameter = round(ev / npar, 2),
      zph_p_female = pf, zph_p_female_fmt = fmt_p(pf),
      zph_p_global = pg, zph_p_global_fmt = fmt_p(pg)))
  }
}

# --- 전기간 PH 검정 (창별 값과 혼동하지 않도록 함께 보고) -------------------
for (key in c("MACCE_pub", "DMS")) {
  d <- cc_subset(make_surv(an, key), COV_M4)
  m <- cox_model(d, COV_M4)
  z <- survival::cox.zph(m)
  pf <- unname(z$table["female", "p"]); pg <- unname(z$table["GLOBAL", "p"])
  npar <- length(stats::coef(m)); ev <- m$nevent
  addA7(rw("A-7", OUT[[key]]$lab, "Whole follow-up (no window split), Model 4",
    m$n, sum(d$female == 1), ev / npar, NA_real_, NA_real_, sprintf("%.2f", ev / npar), fmt_p(pf),
    paste0("Scaled-Schoenfeld test for the sex term over the whole follow-up, reproducing table8a. ",
           "The assumption is not satisfied for the primary composite (MACCE) and is satisfied for ",
           "the hard composite, which is why the analysis is partitioned into landmark windows and ",
           "no single whole-period hazard ratio is reported for MACCE. The global test is ",
           "significant in both cases because of the covariates, not the sex term."),
    measure = "events per parameter", window = "whole follow-up", covariate_set = "Model 4",
    n_events = ev, n_parameters = npar, events_per_parameter = round(ev / npar, 2),
    zph_p_female = pf, zph_p_female_fmt = fmt_p(pf),
    zph_p_global = pg, zph_p_global_fmt = fmt_p(pg)))
}

# --- Model 4 VIF 요약 (tableS8b 재사용) -------------------------------------
vifs <- utils::read.csv(file.path(DIR_TAB, "tableS8b_vif_model4.csv"), stringsAsFactors = FALSE)
addA7(rw("A-7", "Collinearity (Model 4 design matrix)", "Maximum variance inflation factor",
  nrow(an), sum(an$female == 1), max(vifs$VIF), NA_real_, NA_real_,
  sprintf("%.2f (%s)", max(vifs$VIF), vifs$term[which.max(vifs$VIF)]), NA_character_,
  paste0("Variance inflation factors for the Model 4 design matrix (generalized VIF squared for ",
         "multi-level factors), reproduced from tableS8b. No covariate reaches the conventional ",
         "threshold of 5. Source: scripts/04_primary_analysis.R."),
  measure = "variance inflation factor",
  n_covariates_vif_ge5 = sum(vifs$VIF >= 5)))
addA7(rw("A-7", "Collinearity (Model 4 design matrix)", "Variance inflation factor of the sex term",
  nrow(an), sum(an$female == 1), vifs$VIF[vifs$term == "female"], NA_real_, NA_real_,
  sprintf("%.2f", vifs$VIF[vifs$term == "female"]), NA_character_,
  "Reproduced from tableS8b; the exposure is nearly orthogonal to the adjustment set.",
  measure = "variance inflation factor"))

# --- 연령 선형성 LR 검정 (DMS = 기존 보고값 재현, MACCE = 추가) -------------
for (key in c("DMS", "MACCE_pub")) {
  d <- cc_subset(make_surv(an, key), COV_M4)
  f1 <- stats::as.formula(paste0("Surv(.time, .event) ~ ",
                                 paste(c("female", COV_M4), collapse = " + ")))
  f2 <- stats::as.formula(paste0("Surv(.time, .event) ~ ",
                                 paste(c("female", "splines::ns(Age, 3)", setdiff(COV_M4, "Age")),
                                       collapse = " + ")))
  m1 <- survival::coxph(f1, data = d, ties = "efron")
  m2 <- survival::coxph(f2, data = d, ties = "efron")
  a <- stats::anova(m1, m2); pv <- a$`Pr(>|Chi|)`[2]
  addA7(rw("A-7", OUT[[key]]$lab, "Age linearity: likelihood-ratio test, natural cubic spline (3 df) vs linear",
    m1$n, sum(d$female == 1), NA_real_, NA_real_, NA_real_, NA_character_, fmt_p(pv),
    paste0("Linearity in age is rejected for the primary composite; the spline-age specification is ",
           "reported as a sensitivity analysis (Model 4S) and the sex estimate is unchanged. ",
           "Non-robust likelihood-ratio test, as required for nested model comparison."),
    measure = "likelihood-ratio test P value"))
}

A7 <- dplyr::bind_rows(a7_rows)
write.csv(A7, file.path(DIR_TAB, "table32_A7_model_diagnostics_by_window.csv"),
          row.names = FALSE, na = "")
logmsg("[A-7] wrote %d rows to table32_A7_model_diagnostics_by_window.csv", nrow(A7))

###############################################################################
# A-5 (R2-M5) — >3년 창 IPTW 진단 (rebuttal 전용)
#   ESS(성별) · 가중치 백분위 · 사후 표준화평균차(SMD) · 1/99 절단 추정치
#   가중 방식은 scripts/13_v2_landmark_remedies.R (R1-IPTW) 와 동일해야 한다.
###############################################################################

message("\n================ A-5: late-window IPTW diagnostics ================")

a5_rows <- list()
addA5 <- function(x) a5_rows[[length(a5_rows) + 1]] <<- x

d5 <- cc_subset(window_data(an, "MACCE_pub", 1095, Inf), COV_M4)
ps5 <- stats::glm(stats::as.formula(paste("female ~", paste(COV_M4, collapse = " + "))),
                  data = d5, family = stats::binomial())
pp <- stats::fitted(ps5); pf <- mean(d5$female)
d5$w <- ifelse(d5$female == 1, pf / pp, (1 - pf) / (1 - pp))
qs <- stats::quantile(d5$w, c(.01, .99)); d5$wt <- pmin(pmax(d5$w, qs[1]), qs[2])
logmsg("[A-5] late window: n = %d (women %d), events %d", nrow(d5), sum(d5$female == 1),
       sum(d5$.event))

IPTW_NOTE <- paste0(
  "Stabilized inverse-probability weights for the >3-year window (landmark at 3 years), weight ",
  "model = the Model 4 covariates, as in the sensitivity analysis reported in the submitted ",
  "Supplement. Sex is not an intervention: these weights standardize the observed covariate ",
  "distribution and the result is not a causal contrast. Rebuttal material: not proposed for the ",
  "manuscript.")

for (wv in c("w", "wt")) {
  m <- survival::coxph(Surv(.time, .event) ~ female, data = d5, weights = d5[[wv]],
                       cluster = d5$center_anon, robust = TRUE, ties = "efron")
  s <- summary(m)
  addA5(rw("A-5", OUT$MACCE_pub$lab,
    ifelse(wv == "w", "IPTW, untrimmed (reproduction of the reported late-window estimate)",
           "IPTW, trimmed at the 1st and 99th weight percentiles"),
    nrow(d5), sum(d5$female == 1), unname(s$conf.int[1, 1]), unname(s$conf.int[1, 3]),
    unname(s$conf.int[1, 4]), fmt2(s$conf.int[1, 1], s$conf.int[1, 3], s$conf.int[1, 4]),
    fmt_p(unname(s$coefficients[1, "Pr(>|z|)"])), IPTW_NOTE,
    measure = "hazard ratio (women vs men)", n_events = sum(d5$.event)))
}
hr_un <- summary(survival::coxph(Surv(.time, .event) ~ female, data = d5, weights = d5$w,
                                 cluster = d5$center_anon, robust = TRUE))$conf.int[1, 1]
stopifnot(abs(hr_un - 0.83) < 0.015)      # table23 R1-IPTW 재현 assertion

# --- ESS / 가중치 분포 ------------------------------------------------------
ess <- function(w) sum(w)^2 / sum(w^2)
for (grp in list(list(f = NULL, lab = "all patients"), list(f = 0, lab = "men"),
                 list(f = 1, lab = "women"))) {
  w <- if (is.null(grp$f)) d5$w else d5$w[d5$female == grp$f]
  nn <- length(w)
  addA5(rw("A-5", "Weight diagnostics", sprintf("Effective sample size, %s", grp$lab),
    nn, sum(d5$female == 1), ess(w), NA_real_, NA_real_,
    sprintf("%.0f of %d (%.0f%%)", ess(w), nn, 100 * ess(w) / nn), NA_character_,
    paste0("Effective sample size = (sum of weights)^2 / sum of squared weights. ", IPTW_NOTE),
    measure = "effective sample size", group = grp$lab))
  qv <- stats::quantile(w, c(0, .01, .5, .99, 1))
  addA5(rw("A-5", "Weight diagnostics", sprintf("Weight percentiles (min, 1st, 50th, 99th, max), %s", grp$lab),
    nn, sum(d5$female == 1), unname(qv[5]), unname(qv[1]), unname(qv[5]),
    sprintf("%.2f / %.2f / %.2f / %.2f / %.2f", qv[1], qv[2], qv[3], qv[4], qv[5]), NA_character_,
    paste0("Untrimmed stabilized weights. ", IPTW_NOTE),
    measure = "weight percentiles", group = grp$lab,
    w_min = unname(qv[1]), w_p1 = unname(qv[2]), w_p50 = unname(qv[3]),
    w_p99 = unname(qv[4]), w_max = unname(qv[5])))
}

# --- 사후 SMD (Model 4 공변량 전개) ----------------------------------------
Xb <- stats::model.matrix(stats::as.formula(paste("~", paste(COV_M4, collapse = " + "))),
                          data = d5)[, -1, drop = FALSE]
smd_tab <- do.call(rbind, lapply(colnames(Xb), function(v) {
  x <- Xb[, v]; i1 <- d5$female == 1; i0 <- !i1
  sdp <- sqrt((stats::var(x[i1]) + stats::var(x[i0])) / 2)          # 미가중 pooled SD (분모)
  un <- (mean(x[i1]) - mean(x[i0])) / sdp
  wm <- function(ix, w) sum(x[ix] * w[ix]) / sum(w[ix])
  aw <- (wm(i1, d5$w) - wm(i0, d5$w)) / sdp
  tw <- (wm(i1, d5$wt) - wm(i0, d5$wt)) / sdp
  data.frame(variable = v, smd_unweighted = un, smd_weighted = aw, smd_trimmed = tw,
             stringsAsFactors = FALSE)
}))
for (i in seq_len(nrow(smd_tab))) {
  addA5(rw("A-5", "Covariate balance after weighting", smd_tab$variable[i],
    nrow(d5), sum(d5$female == 1), smd_tab$smd_weighted[i], NA_real_, NA_real_,
    sprintf("%+.3f (unweighted %+.3f; trimmed %+.3f)", smd_tab$smd_weighted[i],
            smd_tab$smd_unweighted[i], smd_tab$smd_trimmed[i]), NA_character_,
    paste0("Standardized mean difference (women minus men) after weighting, with the unweighted ",
           "pooled standard deviation as the denominator throughout so that weighted and ",
           "unweighted values are on one scale. ", IPTW_NOTE),
    measure = "standardized mean difference", smd_unweighted = smd_tab$smd_unweighted[i],
    smd_weighted = smd_tab$smd_weighted[i], smd_trimmed = smd_tab$smd_trimmed[i]))
}
imax <- which.max(abs(smd_tab$smd_weighted))
addA5(rw("A-5", "Covariate balance after weighting", "Maximum absolute standardized mean difference",
  nrow(d5), sum(d5$female == 1), abs(smd_tab$smd_weighted[imax]), NA_real_, NA_real_,
  sprintf("%.3f (%s)", abs(smd_tab$smd_weighted[imax]), smd_tab$variable[imax]), NA_character_,
  paste0("Largest absolute standardized mean difference across all ", nrow(smd_tab),
         " expanded Model 4 covariate terms after weighting (unweighted maximum ",
         sprintf("%.3f", max(abs(smd_tab$smd_unweighted))), "). ", IPTW_NOTE),
  measure = "standardized mean difference",
  n_covariate_terms = nrow(smd_tab),
  n_terms_smd_gt_0.1_weighted = sum(abs(smd_tab$smd_weighted) > 0.1)))
logmsg("[A-5] max |SMD| after weighting = %.3f (%s); terms >0.1 = %d of %d",
       abs(smd_tab$smd_weighted[imax]), smd_tab$variable[imax],
       sum(abs(smd_tab$smd_weighted) > 0.1), nrow(smd_tab))

A5 <- dplyr::bind_rows(a5_rows)
write.csv(A5, file.path(DIR_TAB, "table33_A5_iptw_diagnostics_late_window.csv"),
          row.names = FALSE, na = "")
logmsg("[A-5] wrote %d rows to table33_A5_iptw_diagnostics_late_window.csv", nrow(A5))

###############################################################################
# A-2 (R2-M2) 1·2 — 무추적 환자의 센터 분포 + 해당 센터 제외 재적합
#   (1) fudur2 == 0 (n = 779) 및 person-time 미기여(death_fudur == 0 & Death == 0,
#       n = 756) 의 센터 x 성별 교차표 (익명 센터 지표만 보고)
#   (2) 최대 기여 센터를 완전히 제외한 landmark 모형 재적합 (MACCE · DMS)
#   ⛔ (3) delta-adjusted tipping-point 는 연구자 승인 대기 — 미수행
###############################################################################

message("\n================ A-2: no-follow-up centre distribution and exclusion refit ================")

a2_rows <- list()
addA2 <- function(x) a2_rows[[length(a2_rows) + 1]] <<- x

co$no_fu_time <- as.numeric(co$death_fudur == 0 & co$Death == 0)
stopifnot(sum(co$fu_zero) == 779, sum(co$no_fu_time) == 756)

A2_NOTE <- paste0(
  "Two counts are distinguished throughout: 779 patients have no recorded follow-up duration ",
  "beyond the index admission (fudur2 = 0) and 756 of them contribute no person-time to the ",
  "time-to-event models (no recorded follow-up and no event). Centres are identified only by an ",
  "anonymised index; no centre name, code or patient identifier is reported anywhere.")

cen <- co %>% dplyr::group_by(center_anon) %>%
  dplyr::summarise(n = dplyr::n(), n_men = sum(female == 0), n_women = sum(female == 1),
                   z779 = sum(fu_zero), z779_men = sum(fu_zero == 1 & female == 0),
                   z779_women = sum(fu_zero == 1 & female == 1),
                   z756 = sum(no_fu_time), z756_men = sum(no_fu_time == 1 & female == 0),
                   z756_women = sum(no_fu_time == 1 & female == 1), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(z779))
TOPC <- cen$center_anon[1]
rank_by_size <- unname(rank(-cen$n, ties.method = "min")[which(cen$center_anon == TOPC)])
ordinal <- function(k) {
  sfx <- if (k %% 100 %in% 11:13) "th" else switch(as.character(k %% 10),
         "1" = "st", "2" = "nd", "3" = "rd", "th")
  paste0(k, sfx)
}
logmsg("[A-2(1)] centre %d holds %d of the %d patients without recorded follow-up (%.1f%%); that centre has n = %d (%s largest of %d centres) and a within-centre rate of %.1f%% (men %.1f%% / women %.1f%%)",
       TOPC, cen$z779[1], sum(co$fu_zero), 100 * cen$z779[1] / sum(co$fu_zero), cen$n[1],
       ordinal(rank_by_size), nrow(cen), 100 * cen$z779[1] / cen$n[1],
       100 * cen$z779_men[1] / cen$n_men[1], 100 * cen$z779_women[1] / cen$n_women[1])
addA2(rw("A-2(1)", sprintf("Largest contributing centre (anonymised index %d)", TOPC),
  sprintf("Share of all patients with absent follow-up held by one centre (%s largest of %d centres)",
          ordinal(rank_by_size), nrow(cen)),
  cen$n[1], cen$n_women[1], 100 * cen$z779[1] / sum(co$fu_zero), NA_real_, NA_real_,
  sprintf("%d of %d (%.1f%%) under the 779 definition; %d of %d (%.1f%%) under the 756 definition",
          cen$z779[1], sum(co$fu_zero), 100 * cen$z779[1] / sum(co$fu_zero),
          cen$z756[1], sum(co$no_fu_time), 100 * cen$z756[1] / sum(co$no_fu_time)), NA_character_,
  paste0("Within that centre, follow-up is absent for ",
         sprintf("%.1f%% of men and %.1f%% of women", 100 * cen$z779_men[1] / cen$n_men[1],
                 100 * cen$z779_women[1] / cen$n_women[1]),
         ", i.e. the centre reported almost no follow-up for either sex rather than ",
         "differentially by sex. ", A2_NOTE),
  measure = "share of all patients with absent follow-up (%)",
  centre_index = TOPC, centre_size_rank = rank_by_size, n_centres = nrow(cen),
  n_zero_fu_779 = cen$z779[1], n_no_person_time_756 = cen$z756[1],
  share_of_all_779_pct = round(100 * cen$z779[1] / sum(co$fu_zero), 1),
  share_of_all_756_pct = round(100 * cen$z756[1] / sum(co$no_fu_time), 1),
  pct_zero_fu_men = round(100 * cen$z779_men[1] / cen$n_men[1], 2),
  pct_zero_fu_women = round(100 * cen$z779_women[1] / cen$n_women[1], 2)))

for (i in seq_len(nrow(cen))) {
  r <- cen[i, ]
  addA2(rw("A-2(1)", sprintf("Centre %d (anonymised index)", r$center_anon),
    "Patients without recorded follow-up beyond the index admission (n = 779 definition)",
    r$n, r$n_women, 100 * r$z779 / r$n, NA_real_, NA_real_,
    sprintf("%d/%d (%.1f%%)", r$z779, r$n, 100 * r$z779 / r$n), NA_character_, A2_NOTE,
    measure = "percent of the centre's patients",
    centre_index = r$center_anon, n_men = r$n_men, n_women = r$n_women,
    n_zero_fu_779 = r$z779, n_zero_fu_779_men = r$z779_men, n_zero_fu_779_women = r$z779_women,
    pct_zero_fu_men = round(100 * r$z779_men / r$n_men, 2),
    pct_zero_fu_women = round(100 * r$z779_women / r$n_women, 2),
    n_no_person_time_756 = r$z756, n_no_person_time_756_men = r$z756_men,
    n_no_person_time_756_women = r$z756_women,
    share_of_all_779_pct = round(100 * r$z779 / sum(co$fu_zero), 1),
    share_of_all_756_pct = round(100 * r$z756 / sum(co$no_fu_time), 1)))
}

# 전체 / 해당 센터 제외 성별 비교
for (sc in list(list(d = co, lab = "all 33 centres"),
                list(d = co[co$center_anon != TOPC, ], lab = sprintf("excluding centre %d", TOPC)))) {
  for (vv in list(list(v = "fu_zero", lab = "no recorded follow-up (n = 779 definition)"),
                  list(v = "no_fu_time", lab = "no person-time contributed (n = 756 definition)"))) {
    d <- sc$d; x <- d[[vv$v]]
    tb <- table(d$sex_f, x)
    ct <- stats::chisq.test(tb)
    pm <- 100 * mean(x[d$female == 0]); pw <- 100 * mean(x[d$female == 1])
    addA2(rw("A-2(1)", vv$lab, sprintf("Sex comparison, %s", sc$lab),
      nrow(d), sum(d$female == 1), pw - pm, NA_real_, NA_real_,
      sprintf("men %.2f%% vs women %.2f%% (difference %+.2f pp)", pm, pw, pw - pm),
      fmt_p(ct$p.value),
      paste0("Chi-squared test of the sex difference in the proportion with absent follow-up. ",
             "The sex imbalance reported in the submitted manuscript is confined to one centre: ",
             "with that centre removed the proportions are nearly identical. ", A2_NOTE),
      measure = "difference in proportion (percentage points, women - men)",
      pct_men = round(pm, 2), pct_women = round(pw, 2),
      n_men = sum(d$female == 0), n_affected_men = sum(x[d$female == 0]),
      n_affected_women = sum(x[d$female == 1])))
  }
}

## --- (2) 해당 센터 제외 landmark 재적합 -----------------------------------
t22 <- utils::read.csv(file.path(DIR_TAB, "table22_landmark_partition.csv"), stringsAsFactors = FALSE)
map_lab <- c("unadjusted" = "unadjusted", "age only" = "age only",
             "reduced (age, DM, CKD, AF)" = "parsimonious (age, DM, CKD, AF)",
             "Model 4 (re-fitted within window)" = "Model 4 (re-fitted within window)")
an_x <- droplevels(an[an$center_anon != TOPC, , drop = FALSE])
logmsg("[A-2(2)] analysis sample excluding centre %d: n = %d (women %d) of %d",
       TOPC, nrow(an_x), sum(an_x$female), nrow(an))
for (key in c("MACCE_pub", "DMS")) {
  for (w in WINDOWS) {
    for (sp in SPECS) {
      d <- cc_subset(window_data(an_x, key, w$start, w$end), sp$covs)
      wmsg <- character(0)
      m <- withCallingHandlers(
        tryCatch(cox_model(d, sp$covs), error = function(e) NULL),
        warning = function(wn) { wmsg <<- c(wmsg, conditionMessage(wn)); invokeRestart("muffleWarning") })
      if (is.null(m)) next
      s <- summary(m); i <- which(rownames(s$coefficients) == "female")
      pub <- t22[t22$outcome_key == key & t22$window == w$lab &
                   t22$model == map_lab[[sp$lab]], ]
      addA2(rw("A-2(2)", OUT[[key]]$lab,
        sprintf("%s, %s, excluding centre %d", w$lab, sp$lab, TOPC),
        m$n, sum(d$female == 1), unname(s$conf.int[i, 1]), unname(s$conf.int[i, 3]),
        unname(s$conf.int[i, 4]), fmt2(s$conf.int[i, 1], s$conf.int[i, 3], s$conf.int[i, 4]),
        fmt_p(unname(s$coefficients[i, "Pr(>|z|)"])),
        paste0("Landmark model re-fitted after removing every patient of the one centre that ",
               "accounts for most of the absent follow-up. Published counterpart (all centres, ",
               "table22): ", if (nrow(pub)) sprintf("%s, P %s, n = %d", pub$HR_fmt[1], pub$p_fmt[1],
               pub$n[1]) else "not matched", ". Windows condition on being event-free and under ",
               "observation at the window start. ", A2_NOTE,
               if (length(wmsg)) paste0(" Fit warning(s): ",
                                        paste(trimws(unique(wmsg)), collapse = " | ")) else ""),
        measure = "hazard ratio (women vs men)", window = w$lab, covariate_set = sp$lab,
        n_events = m$nevent,
        published_HR_fmt = if (nrow(pub)) pub$HR_fmt[1] else NA_character_,
        published_p_fmt = if (nrow(pub)) pub$p_fmt[1] else NA_character_,
        published_n = if (nrow(pub)) pub$n[1] else NA_integer_,
        published_HR = if (nrow(pub)) pub$HR[1] else NA_real_))
    }
  }
}

A2 <- dplyr::bind_rows(a2_rows)
write.csv(A2, file.path(DIR_TAB, "table34_A2_nofollowup_center_and_refit.csv"),
          row.names = FALSE, na = "")
logmsg("[A-2] wrote %d rows to table34_A2_nofollowup_center_and_refit.csv", nrow(A2))

###############################################################################
# 로그 · sessionInfo
###############################################################################

writeLines(c(sprintf("17_phase7_reviewer_analyses.R  run %s", format(Sys.time())),
             sprintf("seed = %d | g-formula bootstrap reps = %d | Monte Carlo draws = %d",
                     SEED, B_BOOT, B_MC), "", LOG),
           file.path(DIR_REP, "phase7_analyses_log.txt"))
utils::capture.output(utils::sessionInfo(), file = file.path(DIR_REP, "phase7_sessionInfo.txt"))
message("\n== done ==")
