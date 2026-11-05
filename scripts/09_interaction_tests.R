###############################################################################
# 09_interaction_tests.R — A18: sex x lesion-complexity interaction tests
#
# 실행:  Rscript scripts/09_interaction_tests.R
# 입력:  data/processed/kcto_sex_analysis.rds  (via scripts/00_cohort_prep.R)
# 출력:  output/tables/table21_interaction_tests.csv
#        output/tables/table21_interaction_tests.docx
#        interim_report/phase2_stats/reports/interaction_09_log.txt
#        interim_report/phase2_stats/reports/interaction_09_sessionInfo.txt
#
# 코호트: flag_cto_pci & index_yr 2014-2023 -> N = 6,777 (CLAUDE.md §1)
# 분석표본: Model 4 complete case (04/05/06 스크립트와 동일)
#
# ⚠️ 관찰연구다. 모든 추정치는 연관성(association)이며 인과 해석 금지.
# ⛔ 사전 계획을 시사하는 표현 금지 (CLAUDE.md §1-2). 본 분석은 사전등록되지 않았다.
#
#-----------------------------------------------------------------------------
# 왜 필요한가
#-----------------------------------------------------------------------------
# 원고 Discussion 은 "복잡할수록 격차가 커진다"(complexity gradient) 를 해부학적
# 대안설명을 배제하는 결정적 근거("More decisively")로 사용한다. 그러나 지금까지
# 산출된 것은 table15c 의 **층별 개별 추정치를 나란히 제시한 것**뿐이며,
# 두 층의 차이가 우연 범위를 벗어나는지에 대한 검정은 원고에도 표에도 없다.
# 본 스크립트는 그 검정(sex x complexity interaction)을 정식 모형으로 수행한다.
#
#-----------------------------------------------------------------------------
# ★ 판정 기준 — 결과를 보기 전에 기록됨 (본 분석 요청서에 Lead 가 기재한 규칙을
#   그대로 옮긴 것이며, 모형 적합 이후 수정하지 않는다)
#-----------------------------------------------------------------------------
#   주 검정 = 3성분 procedural intensity score (bilateral + retrograde + IVUS)
#             x J-CTO score >=2 (이분), 선형혼합모형 + 센터 랜덤절편,
#             Model 4 공변량(센터 volume 제외).
#
#   (1) 상호작용 P < 0.05                      -> gradient 주장 유지
#   (2) 상호작용 P >= 0.05 이나 방향 일관       -> "directionally consistent but
#                                                 not formally supported" 로 하향
#   (3) 상호작용 방향 역전                      -> gradient 주장 철회
#
#   다중성: 본 스크립트의 검정은 하나의 질문(격차가 복잡도에 따라 달라지는가)에
#   대한 여러 조작화이다. 주 검정 1개를 위와 같이 지정하고, 나머지는 보조로
#   보고하며 별도의 다중비교 보정은 적용하지 않는다(모두 보고하므로 선택적
#   보고 문제는 발생하지 않는다).
#
#-----------------------------------------------------------------------------
# 구성
#   A18-1  주 검정: 3성분 점수 x J-CTO >=2 / x 폐색길이 >=20 mm
#   A18-2  연속형 복잡도: 3성분 점수 x J-CTO (연속) / 4범주 / 스플라인
#   A18-3  Bilateral injection (로지스틱 + 센터 랜덤절편)
#   A18-4  4성분 점수 (민감도)
#   A18-5  table15c 층별 추정치와의 정합성 확인
#   A18-6  gradient 이외의 해부학 반론 근거 2종 정량 평가
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(splines)
})

set.seed(42)

PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "00_cohort_prep.R"))

DIR_TAB  <- file.path(PROJ, "output", "tables")
DIR_REPX <- file.path(PROJ, "interim_report", "phase2_stats", "reports")
dir.create(DIR_TAB,  recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_REPX, recursive = TRUE, showWarnings = FALSE)

LOG <- character(0)
logmsg <- function(...) { s <- sprintf(...); LOG <<- c(LOG, s); message(s) }

DECISION_RULE <- paste0(
  "Decision rule recorded before these models were fitted and not modified afterwards. ",
  "Primary test = three-component procedural intensity score (bilateral injection + retrograde ",
  "attempt + IVUS) x J-CTO score >=2, linear mixed model with a centre random intercept and ",
  "Model 4 covariates (centre volume omitted, absorbed by centre). ",
  "(1) interaction P < 0.05 -> the complexity-gradient claim stands; ",
  "(2) interaction P >= 0.05 with a consistent direction -> the claim is downgraded to ",
  "'directionally consistent but not formally supported'; ",
  "(3) reversed direction -> the claim is withdrawn. ",
  "All other interaction tests are secondary operationalisations of the same question and are ",
  "reported without multiplicity adjustment because all of them are reported.")

CAVEAT <- paste0("Observational secondary analysis of a registry; estimates are associations. ",
                 "These analyses were not prospectively registered.")

###############################################################################
# 0. 코호트 · 점수 구성 (06_intensity_robustness.R 과 동일한 정의)
###############################################################################

co <- build_cohort()
stopifnot(nrow(co) == COHORT_N)

an <- co[stats::complete.cases(co[, COV_M4, drop = FALSE]), , drop = FALSE]
logmsg("[analysis set] Model 4 complete case: N = %d (F = %d)", nrow(an), sum(an$female))

# 성분 정의 — 06 과 완전히 동일해야 한다
IC4 <- c("Bilati_inj", "retro_attempt_any", "IVUS", "dapt_potent")  # 4성분 (민감도)
IC3 <- c("Bilati_inj", "retro_attempt_any", "IVUS")                 # 3성분 (주 지표)

mk_score <- function(dat, vars) {
  M <- as.data.frame(dat[, vars, drop = FALSE])
  ifelse(stats::complete.cases(M), rowSums(M), NA_real_)
}
an$score3 <- mk_score(an, IC3)
an$score4 <- mk_score(an, IC4)

# 복잡도 수정변수 -------------------------------------------------------------
# JCTO_bi 는 codebook 상 "J-CTO score as a binary variable (>=2)" 이며
# JCTO_re >= 2 와 동일함을 assertion 으로 확인한다.
an$jcto_ge2 <- as.numeric(an$JCTO_re >= 2)
stopifnot(all(an$jcto_ge2 == an$JCTO_bi, na.rm = TRUE))
stopifnot(!any(is.na(an$jcto_ge2)))          # JCTO_re 는 COV_M4 에 포함 -> 완전관측
stopifnot(!any(is.na(an$Length20)))
an$len20 <- as.numeric(an$Length20)

# 4범주 J-CTO — 계수명이 안전하도록 라벨을 영숫자로 (JCTO_g4_f 의 ">=3" 라벨은
# 모형 계수명에 부등호가 들어가 선형결합 지정이 취약해진다)
an$jcto_cat <- factor(an$JCTO_re.g4, levels = 0:3,
                      labels = c("c0", "c1", "c2", "c3plus"))
stopifnot(!any(is.na(an$jcto_cat)))

COV_FULL <- COV_M4_NOCENTER
LAB_FULL <- "Model 4 covariates minus centre volume (absorbed by centre)"

###############################################################################
# 1. helper
###############################################################################

.cc_subset <- function(dat, yvar, covs) {
  vars <- unique(c(yvar, "female", "center_anon", covs))
  vars <- intersect(vars, names(dat))
  droplevels(dat[stats::complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE])
}

# 선형결합 추정 (b'L, SE = sqrt(L' V L))
.lincom <- function(b, V, L) {
  est <- sum(b * L)
  se  <- sqrt(as.numeric(t(L) %*% V %*% L))
  p   <- 2 * stats::pnorm(-abs(est / se))
  c(est = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se, p = p)
}

.mkL <- function(nms, which) { L <- rep(0, length(nms)); L[match(which, nms)] <- 1; L }

# 80% 검정력에서 탐지 가능한 최소 상호작용 크기 (양측 alpha=0.05)
.mde80 <- function(se) (stats::qnorm(0.975) + stats::qnorm(0.80)) * se

res <- list()
addrow <- function(x) {
  if (is.null(x) || !nrow(x)) return(invisible(NULL))
  res[[length(res) + 1]] <<- as.data.frame(x, stringsAsFactors = FALSE)
  invisible(NULL)
}

new_row <- function(analysis, outcome, modifier, model_class, cov_lab, d,
                    quantity, measure, est, lo, hi, p, test, df = NA_real_,
                    est_fmt = NULL, note = "") {
  data.frame(
    analysis = analysis, outcome = outcome, modifier = modifier,
    model_class = model_class, covariate_set = cov_lab,
    n = if (is.null(d)) NA_integer_ else nrow(d),
    n_male = if (is.null(d)) NA_integer_ else sum(d$female == 0),
    n_female = if (is.null(d)) NA_integer_ else sum(d$female == 1),
    quantity = quantity, measure = measure,
    estimate = est, lo = lo, hi = hi,
    est_fmt = if (!is.null(est_fmt)) est_fmt else
      ifelse(is.na(est), NA_character_, sprintf("%.3f (%.3f to %.3f)", est, lo, hi)),
    p_value = p, p_fmt = fmt_p(p), test = test, df = df,
    note = note, stringsAsFactors = FALSE, row.names = NULL)
}

#-----------------------------------------------------------------------------
# 이분 수정변수 x female 상호작용 (선형혼합 / 로지스틱 GLMM 공용)
#   반환: 층별 female 효과 2행 + 상호작용 1행 + LRT 1행
#-----------------------------------------------------------------------------
fit_interaction <- function(dat, yvar, modvar, covs, analysis, outcome_lab,
                            modifier_lab, family = c("gaussian", "binomial"),
                            cov_lab = LAB_FULL, note = "",
                            mod_levels = c("modifier = 0", "modifier = 1")) {
  family <- match.arg(family)
  covs_use <- setdiff(covs, modvar)
  d <- .cc_subset(dat, yvar, c(modvar, covs_use))
  d$.y <- d[[yvar]]; d$.m <- d[[modvar]]
  if (nrow(d) < 100) return(NULL)

  rhs_cov <- if (length(covs_use)) paste0(" + ", paste(covs_use, collapse = " + ")) else ""
  f_int  <- stats::as.formula(paste0(".y ~ female * .m", rhs_cov, " + (1 | center_anon)"))
  f_main <- stats::as.formula(paste0(".y ~ female + .m", rhs_cov, " + (1 | center_anon)"))

  if (family == "gaussian") {
    m1 <- lme4::lmer(f_int,  data = d, REML = TRUE)
    m1ml <- lme4::lmer(f_int,  data = d, REML = FALSE)
    m0ml <- lme4::lmer(f_main, data = d, REML = FALSE)
    meas_main <- "adjusted mean difference in score (female - male)"
    meas_int  <- "difference of adjusted mean differences (high complexity - low complexity)"
    tf <- identity
  } else {
    ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5))
    m1 <- lme4::glmer(f_int,  data = d, family = stats::binomial(), nAGQ = 0L, control = ctrl)
    m1ml <- m1
    m0ml <- lme4::glmer(f_main, data = d, family = stats::binomial(), nAGQ = 0L, control = ctrl)
    meas_main <- "OR (female vs male)"
    meas_int  <- "ratio of ORs (high complexity / low complexity)"
    tf <- exp
  }

  b <- lme4::fixef(m1); V <- as.matrix(stats::vcov(m1)); nms <- names(b)
  i_f <- "female"; i_x <- if ("female:.m" %in% nms) "female:.m" else ".m:female"

  e0 <- .lincom(b, V, .mkL(nms, i_f))
  e1 <- .lincom(b, V, .mkL(nms, c(i_f, i_x)))
  ex <- .lincom(b, V, .mkL(nms, i_x))

  lr <- stats::anova(m0ml, m1ml)
  p_lrt <- lr$`Pr(>Chisq)`[2]; chi <- lr$Chisq[2]; ddf <- lr$Df[2]

  sig2 <- as.numeric(lme4::VarCorr(m1)$center_anon[1])
  nt <- paste0(note, sprintf(" | centre random-intercept variance = %.4f", sig2))

  fmt <- function(z) if (family == "gaussian")
    sprintf("%.3f (%.3f to %.3f)", z["est"], z["lo"], z["hi"]) else
    sprintf("%.2f (%.2f-%.2f)", exp(z["est"]), exp(z["lo"]), exp(z["hi"]))

  out <- list(
    new_row(analysis, outcome_lab, modifier_lab, mclass(family), cov_lab, d,
            paste0("female effect, ", mod_levels[1]), meas_main,
            tf(e0["est"]), tf(e0["lo"]), tf(e0["hi"]), e0["p"], "Wald z", 1,
            est_fmt = fmt(e0), note = nt),
    new_row(analysis, outcome_lab, modifier_lab, mclass(family), cov_lab, d,
            paste0("female effect, ", mod_levels[2]), meas_main,
            tf(e1["est"]), tf(e1["lo"]), tf(e1["hi"]), e1["p"], "Wald z", 1,
            est_fmt = fmt(e1), note = nt),
    new_row(analysis, outcome_lab, modifier_lab, mclass(family), cov_lab, d,
            "INTERACTION sex x complexity", meas_int,
            tf(ex["est"]), tf(ex["lo"]), tf(ex["hi"]), ex["p"], "Wald z", 1,
            est_fmt = fmt(ex),
            note = paste0(nt, sprintf(paste0(" | SE = %.4f | smallest interaction detectable with 80%% ",
                                             "power at two-sided alpha 0.05 = %.3f on this scale"),
                                      ex["se"], .mde80(ex["se"])))),
    new_row(analysis, outcome_lab, modifier_lab, mclass(family), cov_lab, d,
            "INTERACTION sex x complexity", "likelihood-ratio test",
            NA_real_, NA_real_, NA_real_, p_lrt, "LRT", ddf,
            est_fmt = sprintf("chi-square = %.2f, df = %d", chi, ddf), note = nt))
  dplyr::bind_rows(out)
}

mclass <- function(family) {
  if (family == "gaussian") "Linear mixed model, centre random intercept"
  else "Logistic GLMM, centre random intercept (nAGQ = 0)"
}

#-----------------------------------------------------------------------------
# 연속 / 다범주 수정변수 x female (LRT + 지정 지점의 female 효과)
#-----------------------------------------------------------------------------
fit_interaction_gen <- function(dat, yvar, mod_term, covs_drop, covs, analysis,
                                outcome_lab, modifier_lab, family = c("gaussian", "binomial"),
                                at = NULL, at_labs = NULL, cov_lab = LAB_FULL, note = "",
                                extra_vars = character(0)) {
  family <- match.arg(family)
  covs_use <- setdiff(covs, covs_drop)
  d <- .cc_subset(dat, yvar, c(covs_use, extra_vars))
  d$.y <- d[[yvar]]
  if (nrow(d) < 100) return(NULL)

  rhs_cov <- if (length(covs_use)) paste0(" + ", paste(covs_use, collapse = " + ")) else ""
  f_int  <- stats::as.formula(paste0(".y ~ female * ", mod_term, rhs_cov, " + (1 | center_anon)"))
  f_main <- stats::as.formula(paste0(".y ~ female + ", mod_term, rhs_cov, " + (1 | center_anon)"))

  if (family == "gaussian") {
    m1   <- lme4::lmer(f_int,  data = d, REML = TRUE)
    m1ml <- lme4::lmer(f_int,  data = d, REML = FALSE)
    m0ml <- lme4::lmer(f_main, data = d, REML = FALSE)
    meas <- "adjusted mean difference in score (female - male)"
    tf <- identity
    fmt <- function(z) sprintf("%.3f (%.3f to %.3f)", z["est"], z["lo"], z["hi"])
  } else {
    ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5))
    m1   <- lme4::glmer(f_int,  data = d, family = stats::binomial(), nAGQ = 0L, control = ctrl)
    m1ml <- m1
    m0ml <- lme4::glmer(f_main, data = d, family = stats::binomial(), nAGQ = 0L, control = ctrl)
    meas <- "OR (female vs male)"
    tf <- exp
    fmt <- function(z) sprintf("%.2f (%.2f-%.2f)", exp(z["est"]), exp(z["lo"]), exp(z["hi"]))
  }

  b <- lme4::fixef(m1); V <- as.matrix(stats::vcov(m1)); nms <- names(b)
  ix <- grep("^female:|:female$", nms, value = TRUE)
  lr <- stats::anova(m0ml, m1ml)

  out <- list()
  # 상호작용 항이 1개면 Wald 추정치도 제공
  if (length(ix) == 1) {
    ex <- .lincom(b, V, .mkL(nms, ix))
    out[[length(out) + 1]] <- new_row(
      analysis, outcome_lab, modifier_lab, mclass(family), cov_lab, d,
      "INTERACTION sex x complexity",
      if (family == "gaussian") "change in adjusted mean difference per unit of the modifier"
      else "ratio of ORs per unit of the modifier",
      tf(ex["est"]), tf(ex["lo"]), tf(ex["hi"]), ex["p"], "Wald z", 1,
      est_fmt = fmt(ex),
      note = paste0(note, sprintf(paste0(" | SE = %.4f | smallest per-unit interaction detectable ",
                                         "with 80%% power = %.3f on this scale"),
                                  ex["se"], .mde80(ex["se"]))))
  }
  out[[length(out) + 1]] <- new_row(
    analysis, outcome_lab, modifier_lab, mclass(family), cov_lab, d,
    "INTERACTION sex x complexity", "likelihood-ratio test",
    NA_real_, NA_real_, NA_real_, lr$`Pr(>Chisq)`[2], "LRT", lr$Df[2],
    est_fmt = sprintf("chi-square = %.2f, df = %d", lr$Chisq[2], lr$Df[2]), note = note)

  # 지정된 지점에서의 female 효과 (L 벡터 직접 구성)
  if (!is.null(at)) {
    for (k in seq_along(at)) {
      L <- rep(0, length(nms)); L[match("female", nms)] <- 1
      w <- at[[k]]
      for (nm in names(w)) if (nm %in% nms) L[match(nm, nms)] <- w[[nm]]
      e <- .lincom(b, V, L)
      out[[length(out) + 1]] <- new_row(
        analysis, outcome_lab, modifier_lab, mclass(family), cov_lab, d,
        paste0("female effect, ", at_labs[k]), meas,
        tf(e["est"]), tf(e["lo"]), tf(e["hi"]), e["p"], "Wald z", 1,
        est_fmt = fmt(e), note = note)
    }
  }
  list(rows = dplyr::bind_rows(out), m1 = m1, m1ml = m1ml, d = d)
}

#-----------------------------------------------------------------------------
# 층별 재적합 (table15c 재현) + 두 층 추정치의 naive 대비 검정
#-----------------------------------------------------------------------------
fit_stratum <- function(dat, yvar, covs, family = c("gaussian", "binomial")) {
  family <- match.arg(family)
  d <- .cc_subset(dat, yvar, covs); d$.y <- d[[yvar]]
  rhs <- paste(c("female", covs), collapse = " + ")
  f <- stats::as.formula(paste0(".y ~ ", rhs, " + (1 | center_anon)"))
  m <- if (family == "gaussian") lme4::lmer(f, data = d, REML = TRUE) else
    lme4::glmer(f, data = d, family = stats::binomial(), nAGQ = 0L,
                control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
  s <- summary(m)
  b <- s$coefficients["female", "Estimate"]; se <- s$coefficients["female", "Std. Error"]
  list(b = b, se = se, lo = b - 1.96 * se, hi = b + 1.96 * se,
       p = 2 * stats::pnorm(-abs(b / se)), n = nrow(d),
       n_male = sum(d$female == 0), n_female = sum(d$female == 1), d = d)
}

###############################################################################
# A18-0. 참조 재현 — 3성분 점수 전체 코호트 추정치가 table15c 와 일치하는지
###############################################################################

message("== A18-0: reference reproduction ==")

ref3 <- fit_stratum(an, "score3", COV_FULL, "gaussian")
logmsg("[A18-0] three-component score, full analysis sample: n = %d (F %d), estimate = %.3f (%.3f to %.3f), P = %s",
       ref3$n, ref3$n_female, ref3$b, ref3$lo, ref3$hi, fmt_p(ref3$p))
stopifnot(ref3$n == 5879)
stopifnot(abs(ref3$b - (-0.088)) < 0.002)     # table15c A16-1/A16-5 값 재현
addrow(new_row("A18-0", "Three-component procedural intensity score (0-3)", "none (reference)",
               mclass("gaussian"), LAB_FULL, ref3$d,
               "female effect, whole analysis sample",
               "adjusted mean difference in score (female - male)",
               ref3$b, ref3$lo, ref3$hi, ref3$p, "Wald z", 1,
               note = "reference estimate reproduced from table15c (A16-1 / A16-5)"))

ref4 <- fit_stratum(an, "score4", COV_FULL, "gaussian")
logmsg("[A18-0] four-component score, full analysis sample: n = %d (F %d), estimate = %.3f (%.3f to %.3f), P = %s",
       ref4$n, ref4$n_female, ref4$b, ref4$lo, ref4$hi, fmt_p(ref4$p))
stopifnot(abs(ref4$b - (-0.120)) < 0.002)
addrow(new_row("A18-0", "Four-component procedural intensity score (0-4)", "none (reference)",
               mclass("gaussian"), LAB_FULL, ref4$d,
               "female effect, whole analysis sample",
               "adjusted mean difference in score (female - male)",
               ref4$b, ref4$lo, ref4$hi, ref4$p, "Wald z", 1,
               note = "reference estimate reproduced from table15c (A16-0)"))

refb <- fit_stratum(an, "Bilati_inj", COV_FULL, "binomial")
logmsg("[A18-0] bilateral injection, full analysis sample: n = %d (F %d), OR = %.2f (%.2f-%.2f), P = %s",
       refb$n, refb$n_female, exp(refb$b), exp(refb$lo), exp(refb$hi), fmt_p(refb$p))
addrow(new_row("A18-0", "Bilateral contrast injection", "none (reference)",
               mclass("binomial"), LAB_FULL, refb$d,
               "female effect, whole analysis sample", "OR (female vs male)",
               exp(refb$b), exp(refb$lo), exp(refb$hi), refb$p, "Wald z", 1,
               est_fmt = sprintf("%.2f (%.2f-%.2f)", exp(refb$b), exp(refb$lo), exp(refb$hi)),
               note = "reference estimate reproduced from table14 / table18 (within-centre)"))

###############################################################################
# A18-1. 주 검정 — 3성분 점수 x 이분 복잡도
###############################################################################

message("== A18-1: primary interaction tests (binary complexity) ==")

addrow(fit_interaction(an, "score3", "jcto_ge2", COV_FULL, "A18-1",
                       "Three-component procedural intensity score (0-3)",
                       "J-CTO score >=2 (binary)", "gaussian",
                       note = "PRIMARY TEST of the complexity gradient",
                       mod_levels = c("J-CTO score <2", "J-CTO score >=2")))

addrow(fit_interaction(an, "score3", "len20", COV_FULL, "A18-1",
                       "Three-component procedural intensity score (0-3)",
                       "Occlusion length >=20 mm (binary)", "gaussian",
                       note = "second binary operationalisation of complexity",
                       mod_levels = c("occlusion length <20 mm", "occlusion length >=20 mm")))

###############################################################################
# A18-2. 연속형 복잡도 (검정력 확보) + 다범주 + 스플라인
###############################################################################

message("== A18-2: continuous / categorical / spline complexity ==")

# (a) 연속 J-CTO (0-5), 선형 상호작용
g_lin <- fit_interaction_gen(an, "score3", "JCTO_re", "JCTO_re", COV_FULL, "A18-2",
                             "Three-component procedural intensity score (0-3)",
                             "J-CTO score (continuous, 0-5)", "gaussian",
                             at = list(c(`female:JCTO_re` = 0), c(`female:JCTO_re` = 1),
                                       c(`female:JCTO_re` = 2), c(`female:JCTO_re` = 3),
                                       c(`female:JCTO_re` = 4)),
                             at_labs = c("J-CTO 0", "J-CTO 1", "J-CTO 2", "J-CTO 3", "J-CTO 4"),
                             extra_vars = "JCTO_re",
                             note = "continuous modifier retains more information than dichotomisation")
addrow(g_lin$rows)

# (b) 4범주 J-CTO (0 / 1 / 2 / >=3) — 용량-반응 형태 확인
g_cat <- fit_interaction_gen(an, "score3", "jcto_cat", c("JCTO_re"), COV_FULL, "A18-2",
                             "Three-component procedural intensity score (0-3)",
                             "J-CTO score, 4 categories (0 / 1 / 2 / >=3)", "gaussian",
                             at = list(NULL, c(`female:jcto_catc1` = 1), c(`female:jcto_catc2` = 1),
                                       c(`female:jcto_catc3plus` = 1)),
                             at_labs = c("J-CTO 0 (reference)", "J-CTO 1", "J-CTO 2", "J-CTO >=3"),
                             extra_vars = "jcto_cat",
                             note = "dose-response form of the modifier; 3 df global test")
addrow(g_cat$rows)
stopifnot(all(c("female:jcto_catc1", "female:jcto_catc2", "female:jcto_catc3plus") %in%
                names(lme4::fixef(g_cat$m1))))

# (c) 스플라인 상호작용 (선형성 가정 완화)
g_spl <- fit_interaction_gen(an, "score3", "splines::ns(JCTO_re, df = 3)", "JCTO_re", COV_FULL,
                             "A18-2", "Three-component procedural intensity score (0-3)",
                             "J-CTO score, natural cubic spline (3 df)", "gaussian",
                             extra_vars = "JCTO_re",
                             note = "relaxes linearity of the modifier and of the interaction")
addrow(g_spl$rows)

# (d) 비선형성의 소재 분해 — 주효과인가 상호작용인가
#     (c) 의 3 df 검정은 "비선형 주효과를 허용한 상태에서의 상호작용" 이다.
#     여기서는 주효과 자체의 비선형성(2 df)과 선형 대비 스플라인 전체(4 df)를 함께 낸다.
d_sp <- .cc_subset(an, "score3", COV_FULL); d_sp$.y <- d_sp$score3
cov_nj <- setdiff(COV_FULL, "JCTO_re")
.f <- function(s) stats::as.formula(paste0(".y ~ ", s, " + ", paste(cov_nj, collapse = " + "),
                                           " + (1 | center_anon)"))
m_ml_lin_main <- lme4::lmer(.f("female + JCTO_re"), data = d_sp, REML = FALSE)
m_ml_spl_main <- lme4::lmer(.f("female + splines::ns(JCTO_re, df = 3)"), data = d_sp, REML = FALSE)
m_ml_lin_int  <- lme4::lmer(.f("female * JCTO_re"), data = d_sp, REML = FALSE)
m_ml_spl_int  <- lme4::lmer(.f("female * splines::ns(JCTO_re, df = 3)"), data = d_sp, REML = FALSE)

lr_main <- stats::anova(m_ml_lin_main, m_ml_spl_main)
addrow(new_row("A18-2", "Three-component procedural intensity score (0-3)",
               "J-CTO score, main effect only", mclass("gaussian"), LAB_FULL, d_sp,
               "departure of the MAIN effect of J-CTO from linearity", "likelihood-ratio test",
               NA_real_, NA_real_, NA_real_, lr_main$`Pr(>Chisq)`[2], "LRT", lr_main$Df[2],
               est_fmt = sprintf("chi-square = %.2f, df = %d", lr_main$Chisq[2], lr_main$Df[2]),
               note = paste0("no sex interaction in either model; locates any non-linearity in the ",
                             "association of complexity with intensity itself")))

lr_sp <- stats::anova(m_ml_lin_int, m_ml_spl_int)
addrow(new_row("A18-2", "Three-component procedural intensity score (0-3)",
               "J-CTO score, spline vs linear specification", mclass("gaussian"), LAB_FULL, d_sp,
               "joint departure from linearity (main effect and interaction)", "likelihood-ratio test",
               NA_real_, NA_real_, NA_real_, lr_sp$`Pr(>Chisq)`[2], "LRT", lr_sp$Df[2],
               est_fmt = sprintf("chi-square = %.2f, df = %d", lr_sp$Chisq[2], lr_sp$Df[2]),
               note = paste0("this test is NOT specific to the interaction; the interaction-specific ",
                             "component is the 3 df spline interaction test reported above")))

###############################################################################
# A18-3. Bilateral injection — 유일하게 살아남은 개별 지표
###############################################################################

message("== A18-3: bilateral injection ==")

addrow(fit_interaction(an, "Bilati_inj", "jcto_ge2", COV_FULL, "A18-3",
                       "Bilateral contrast injection", "J-CTO score >=2 (binary)", "binomial",
                       note = "the only individual indicator that survived the earlier challenges",
                       mod_levels = c("J-CTO score <2", "J-CTO score >=2")))

addrow(fit_interaction(an, "Bilati_inj", "len20", COV_FULL, "A18-3",
                       "Bilateral contrast injection", "Occlusion length >=20 mm (binary)", "binomial",
                       note = "second binary operationalisation of complexity",
                       mod_levels = c("occlusion length <20 mm", "occlusion length >=20 mm")))

g_bil <- fit_interaction_gen(an, "Bilati_inj", "JCTO_re", "JCTO_re", COV_FULL, "A18-3",
                             "Bilateral contrast injection", "J-CTO score (continuous, 0-5)",
                             "binomial",
                             at = list(c(`female:JCTO_re` = 0), c(`female:JCTO_re` = 1),
                                       c(`female:JCTO_re` = 2), c(`female:JCTO_re` = 3)),
                             at_labs = c("J-CTO 0", "J-CTO 1", "J-CTO 2", "J-CTO 3"),
                             extra_vars = "JCTO_re",
                             note = "continuous modifier")
addrow(g_bil$rows)

###############################################################################
# A18-4. 4성분 점수 (민감도)
###############################################################################

message("== A18-4: four-component score (sensitivity) ==")

addrow(fit_interaction(an, "score4", "jcto_ge2", COV_FULL, "A18-4",
                       "Four-component procedural intensity score (0-4)",
                       "J-CTO score >=2 (binary)", "gaussian",
                       note = "sensitivity: adds potent P2Y12 at discharge, which is not a procedural step",
                       mod_levels = c("J-CTO score <2", "J-CTO score >=2")))

addrow(fit_interaction(an, "score4", "len20", COV_FULL, "A18-4",
                       "Four-component procedural intensity score (0-4)",
                       "Occlusion length >=20 mm (binary)", "gaussian",
                       note = "sensitivity",
                       mod_levels = c("occlusion length <20 mm", "occlusion length >=20 mm")))

g_lin4 <- fit_interaction_gen(an, "score4", "JCTO_re", "JCTO_re", COV_FULL, "A18-4",
                              "Four-component procedural intensity score (0-4)",
                              "J-CTO score (continuous, 0-5)", "gaussian",
                              extra_vars = "JCTO_re",
                              note = "sensitivity, continuous modifier")
addrow(g_lin4$rows)

###############################################################################
# A18-5. table15c 층별 추정치와의 정합성
#   (i) 층별 모형을 재적합해 table15c 값을 재현
#   (ii) 두 층 추정치의 naive 대비 검정 (모든 공변량 효과가 층별로 다른 모형과 동치)
#   (iii) 상호작용 모형이 함의하는 층별 값과 비교
###############################################################################

message("== A18-5: consistency with the stratified estimates ==")

consist <- list()
for (cfg in list(
  # ref = table15c / table18 이 보고한 층별 값 (재현 검증용)
  list(sc = "score3", sclab = "Three-component procedural intensity score (0-3)",
       fam = "gaussian", ref_lo = -0.066, ref_hi = -0.107, tol = 0.002),
  list(sc = "score4", sclab = "Four-component procedural intensity score (0-4)",
       fam = "gaussian", ref_lo = -0.078, ref_hi = -0.154, tol = 0.002),
  list(sc = "Bilati_inj", sclab = "Bilateral contrast injection",
       fam = "binomial", ref_lo = NA_real_, ref_hi = log(0.71), tol = 0.03))) {

  s_lo <- fit_stratum(an[an$jcto_ge2 == 0, , drop = FALSE], cfg$sc, COV_FULL, cfg$fam)
  s_hi <- fit_stratum(an[an$jcto_ge2 == 1, , drop = FALSE], cfg$sc, COV_FULL, cfg$fam)
  # 재현 검증 — 층별 추정치가 기존 표와 일치해야 한다
  if (!is.na(cfg$ref_lo)) stopifnot(abs(s_lo$b - cfg$ref_lo) < cfg$tol)
  stopifnot(abs(s_hi$b - cfg$ref_hi) < cfg$tol)

  dd  <- s_hi$b - s_lo$b
  sed <- sqrt(s_lo$se^2 + s_hi$se^2)
  pz  <- 2 * stats::pnorm(-abs(dd / sed))

  tf  <- if (cfg$fam == "gaussian") identity else exp
  fm  <- if (cfg$fam == "gaussian")
    function(e, l, h) sprintf("%.3f (%.3f to %.3f)", e, l, h) else
    function(e, l, h) sprintf("%.2f (%.2f-%.2f)", exp(e), exp(l), exp(h))
  meas_main <- if (cfg$fam == "gaussian")
    "adjusted mean difference in score (female - male)" else "OR (female vs male)"

  for (z in list(list(s = s_lo, lab = "J-CTO score <2"), list(s = s_hi, lab = "J-CTO score >=2"))) {
    addrow(new_row("A18-5", cfg$sclab, "J-CTO score >=2 (binary)", mclass(cfg$fam), LAB_FULL, NULL,
                   paste0("STRATIFIED refit, female effect, ", z$lab), meas_main,
                   tf(z$s$b), tf(z$s$lo), tf(z$s$hi), z$s$p, "Wald z", 1,
                   est_fmt = fm(z$s$b, z$s$lo, z$s$hi),
                   note = paste0("separate model fitted within the stratum; all covariate effects, ",
                                 "the residual variance and the centre random-intercept variance are ",
                                 "free to differ between strata (this is what table15c reports)")))
    res[[length(res)]]$n <- z$s$n
    res[[length(res)]]$n_male <- z$s$n_male
    res[[length(res)]]$n_female <- z$s$n_female
  }
  addrow(new_row("A18-5", cfg$sclab, "J-CTO score >=2 (binary)", "Contrast of two stratified models",
                 LAB_FULL, NULL, "INTERACTION sex x complexity (fully stratified contrast)",
                 if (cfg$fam == "gaussian")
                   "difference of adjusted mean differences (high complexity - low complexity)" else
                   "ratio of ORs (high complexity / low complexity)",
                 tf(dd), tf(dd - 1.96 * sed), tf(dd + 1.96 * sed), pz, "Wald z", 1,
                 est_fmt = fm(dd, dd - 1.96 * sed, dd + 1.96 * sed),
                 note = sprintf(paste0("z = %.2f | SE = %.4f | equivalent to a model in which every ",
                                       "covariate is allowed to interact with the stratum; compare with ",
                                       "the single-interaction-term model in A18-1"), dd / sed, sed)))
  res[[length(res)]]$n <- s_lo$n + s_hi$n
  res[[length(res)]]$n_male <- s_lo$n_male + s_hi$n_male
  res[[length(res)]]$n_female <- s_lo$n_female + s_hi$n_female

  consist[[cfg$sc]] <- list(lo = s_lo, hi = s_hi, diff = dd, se = sed, p = pz)
  logmsg("[A18-5] %s : <2 = %.3f (SE %.4f, n %d) | >=2 = %.3f (SE %.4f, n %d) | difference %.3f (SE %.4f), z = %.2f, P = %s",
         cfg$sc, s_lo$b, s_lo$se, s_lo$n, s_hi$b, s_hi$se, s_hi$n, dd, sed, dd / sed, fmt_p(pz))
}

# 폐색길이 층 (table15c 는 >=20 mm 층만 보고했으므로 <20 mm 대비층을 새로 적합)
for (cfg in list(
  list(sc = "score3", sclab = "Three-component procedural intensity score (0-3)", fam = "gaussian"),
  list(sc = "Bilati_inj", sclab = "Bilateral contrast injection", fam = "binomial"))) {
  s_lo <- fit_stratum(an[an$len20 == 0, , drop = FALSE], cfg$sc, COV_FULL, cfg$fam)
  s_hi <- fit_stratum(an[an$len20 == 1, , drop = FALSE], cfg$sc, COV_FULL, cfg$fam)
  dd <- s_hi$b - s_lo$b; sed <- sqrt(s_lo$se^2 + s_hi$se^2); pz <- 2 * stats::pnorm(-abs(dd / sed))
  tf <- if (cfg$fam == "gaussian") identity else exp
  fm <- if (cfg$fam == "gaussian")
    function(e, l, h) sprintf("%.3f (%.3f to %.3f)", e, l, h) else
    function(e, l, h) sprintf("%.2f (%.2f-%.2f)", exp(e), exp(l), exp(h))
  meas_main <- if (cfg$fam == "gaussian")
    "adjusted mean difference in score (female - male)" else "OR (female vs male)"
  for (z in list(list(s = s_lo, lab = "occlusion length <20 mm"),
                 list(s = s_hi, lab = "occlusion length >=20 mm"))) {
    addrow(new_row("A18-5", cfg$sclab, "Occlusion length >=20 mm (binary)", mclass(cfg$fam),
                   LAB_FULL, NULL, paste0("STRATIFIED refit, female effect, ", z$lab), meas_main,
                   tf(z$s$b), tf(z$s$lo), tf(z$s$hi), z$s$p, "Wald z", 1,
                   est_fmt = fm(z$s$b, z$s$lo, z$s$hi),
                   note = paste0("the <20 mm contrast stratum was not reported in table15c and is ",
                                 "fitted here so that the two strata can be compared")))
    res[[length(res)]]$n <- z$s$n
    res[[length(res)]]$n_male <- z$s$n_male
    res[[length(res)]]$n_female <- z$s$n_female
  }
  addrow(new_row("A18-5", cfg$sclab, "Occlusion length >=20 mm (binary)",
                 "Contrast of two stratified models", LAB_FULL, NULL,
                 "INTERACTION sex x complexity (fully stratified contrast)",
                 if (cfg$fam == "gaussian")
                   "difference of adjusted mean differences (high complexity - low complexity)" else
                   "ratio of ORs (high complexity / low complexity)",
                 tf(dd), tf(dd - 1.96 * sed), tf(dd + 1.96 * sed), pz, "Wald z", 1,
                 est_fmt = fm(dd, dd - 1.96 * sed, dd + 1.96 * sed),
                 note = sprintf("z = %.2f | SE = %.4f", dd / sed, sed)))
  res[[length(res)]]$n <- s_lo$n + s_hi$n
  res[[length(res)]]$n_male <- s_lo$n_male + s_hi$n_male
  res[[length(res)]]$n_female <- s_lo$n_female + s_hi$n_female
  logmsg("[A18-5] %s (length >=20 mm) : <20 = %.3f (SE %.4f) | >=20 = %.3f (SE %.4f) | difference %.3f, z = %.2f, P = %s",
         cfg$sc, s_lo$b, s_lo$se, s_hi$b, s_hi$se, dd, dd / sed, fmt_p(pz))
}

# 미보정 층별 평균 (해석 보조) -----------------------------------------------
desc <- an[!is.na(an$score3), ]
for (g in 0:1) {
  dg <- desc[desc$jcto_ge2 == g, ]
  logmsg("[A18-5] crude mean three-component score, J-CTO %s : men %.3f (n %d) / women %.3f (n %d), crude difference %.3f",
         ifelse(g == 1, ">=2", "<2"), mean(dg$score3[dg$female == 0]), sum(dg$female == 0),
         mean(dg$score3[dg$female == 1]), sum(dg$female == 1),
         mean(dg$score3[dg$female == 1]) - mean(dg$score3[dg$female == 0]))
  addrow(new_row("A18-5", "Three-component procedural intensity score (0-3)",
                 "J-CTO score >=2 (binary)", "Crude means (no adjustment)", "none", dg,
                 paste0("crude mean score difference, J-CTO ", ifelse(g == 1, ">=2", "<2")),
                 "crude mean difference in score (female - male)",
                 mean(dg$score3[dg$female == 1]) - mean(dg$score3[dg$female == 0]),
                 NA_real_, NA_real_, NA_real_, "descriptive", NA_real_,
                 est_fmt = sprintf("men %.3f, women %.3f, difference %.3f",
                                   mean(dg$score3[dg$female == 0]), mean(dg$score3[dg$female == 1]),
                                   mean(dg$score3[dg$female == 1]) - mean(dg$score3[dg$female == 0])),
                 note = "descriptive only"))
}

###############################################################################
# A18-5b. 불일치의 원인 규명 — 공변량 제약 사다리
#
#   A18-1(단일 상호작용항, 공변량 효과를 층 간 공통으로 제약)과
#   A18-5(층별 개별 적합의 대비, 모든 공변량 효과가 층별로 자유)는 같은 질문에
#   대해 서로 다른 답을 준다. 두 극단 사이를 단계적으로 이어 원인을 규명한다.
#     S1 공변량 효과 공통                      (= A18-1)
#     S2 + 수정변수 x 병변 공변량 상호작용
#     S3 + 수정변수 x 모든 Model 4 공변량      (= A18-5 의 층별 대비와 사실상 동치)
#     S4 (폐색길이 전용) 공변량에서 J-CTO 제거 — 폐색길이는 J-CTO 의 구성요소이므로
#        같은 정보를 두 번 넣은 설정이 아닌지 확인
###############################################################################

message("== A18-5b: covariate-specification ladder ==")

spec_interaction <- function(dat, yvar, modvar, covs, mod_x, analysis, outcome_lab,
                             modifier_lab, spec_lab, family = c("gaussian", "binomial"),
                             note = "") {
  family <- match.arg(family)
  covs_use <- setdiff(covs, modvar)
  d <- .cc_subset(dat, yvar, c(modvar, covs_use)); d$.y <- d[[yvar]]; d$.m <- d[[modvar]]
  x_terms <- setdiff(mod_x, modvar)
  base <- paste0(".m + ", paste(covs_use, collapse = " + "))
  if (length(x_terms)) base <- paste0(base, " + .m:(", paste(x_terms, collapse = " + "), ")")
  f_int  <- stats::as.formula(paste0(".y ~ female + ", base, " + female:.m + (1 | center_anon)"))
  f_main <- stats::as.formula(paste0(".y ~ female + ", base, " + (1 | center_anon)"))
  if (family == "gaussian") {
    m1 <- lme4::lmer(f_int, data = d, REML = TRUE)
    lr <- stats::anova(lme4::lmer(f_main, data = d, REML = FALSE),
                       lme4::lmer(f_int,  data = d, REML = FALSE))
    tf <- identity
    fm <- function(z) sprintf("%.3f (%.3f to %.3f)", z["est"], z["lo"], z["hi"])
    meas <- "difference of adjusted mean differences (high complexity - low complexity)"
  } else {
    ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5))
    m1 <- lme4::glmer(f_int, data = d, family = stats::binomial(), nAGQ = 0L, control = ctrl)
    lr <- stats::anova(lme4::glmer(f_main, data = d, family = stats::binomial(), nAGQ = 0L,
                                   control = ctrl), m1)
    tf <- exp
    fm <- function(z) sprintf("%.2f (%.2f-%.2f)", exp(z["est"]), exp(z["lo"]), exp(z["hi"]))
    meas <- "ratio of ORs (high complexity / low complexity)"
  }
  b <- lme4::fixef(m1); V <- as.matrix(stats::vcov(m1)); nms <- names(b)
  ix <- if ("female:.m" %in% nms) "female:.m" else ".m:female"
  ex <- .lincom(b, V, .mkL(nms, ix))
  logmsg("[A18-5b] %s | %s | %s : interaction %s, P = %s (LRT P = %s)",
         yvar, modifier_lab, spec_lab, fm(ex), fmt_p(ex["p"]), fmt_p(lr$`Pr(>Chisq)`[2]))
  new_row(analysis, outcome_lab, modifier_lab,
          paste0(mclass(family), " | ", spec_lab), LAB_FULL, d,
          "INTERACTION sex x complexity", meas,
          tf(ex["est"]), tf(ex["lo"]), tf(ex["hi"]), ex["p"], "Wald z", 1,
          est_fmt = fm(ex),
          note = paste0(note, sprintf(" | likelihood-ratio P = %s | number of fixed effects = %d",
                                      fmt_p(lr$`Pr(>Chisq)`[2]), length(b))))
}

SPEC <- list(
  list(lab = "S1: covariate effects common to both strata", mx = character(0),
       nt = "identical to the model reported in A18-1 / A18-3 / A18-4"),
  list(lab = "S2: modifier also interacted with the lesion covariates", mx = COV_LES,
       nt = "allows the lesion covariates to act differently in the two strata"),
  list(lab = "S3: modifier interacted with every Model 4 covariate", mx = COV_M4_NOCENTER,
       nt = paste0("equivalent to fitting the two strata separately except that the residual and ",
                   "centre variance components remain common; converges on the stratified contrast ",
                   "reported in A18-5")))

for (cfg in list(
  list(y = "score3", ylab = "Three-component procedural intensity score (0-3)", fam = "gaussian"),
  list(y = "score4", ylab = "Four-component procedural intensity score (0-4)", fam = "gaussian"),
  list(y = "Bilati_inj", ylab = "Bilateral contrast injection", fam = "binomial"))) {
  for (md in list(list(v = "jcto_ge2", lab = "J-CTO score >=2 (binary)"),
                  list(v = "len20",    lab = "Occlusion length >=20 mm (binary)"))) {
    for (sp in SPEC) {
      addrow(spec_interaction(an, cfg$y, md$v, COV_FULL, sp$mx, "A18-5b", cfg$ylab,
                              md$lab, sp$lab, cfg$fam, sp$nt))
    }
  }
}

# S4: 폐색길이 분석에서 J-CTO 를 공변량에서 제거 (구성요소 이중 투입 확인)
for (cfg in list(
  list(y = "score3", ylab = "Three-component procedural intensity score (0-3)", fam = "gaussian"),
  list(y = "Bilati_inj", ylab = "Bilateral contrast injection", fam = "binomial"))) {
  addrow(spec_interaction(an, cfg$y, "len20", setdiff(COV_FULL, "JCTO_re"), character(0),
                          "A18-5b", cfg$ylab, "Occlusion length >=20 mm (binary)",
                          "S4: J-CTO score removed from the covariates", cfg$fam,
                          paste0("occlusion length >=20 mm is one of the components of the J-CTO ",
                                 "score, so the S1 specification adjusts for the modifier twice; ",
                                 "this specification removes that overlap")))
}

###############################################################################
# A18-6. gradient 이외의 해부학 반론 근거 2종 — 정량 평가
#   근거 1: 신장·체중 공변량 추가 시 3성분 점수 추정치 불변 (-0.088 -> -0.088)
#   근거 2: 신장 공통지지 구간 제한 시 bilateral 0.63-0.68 유지
###############################################################################

message("== A18-6: strength of the two remaining anatomical counter-arguments ==")

# --- 근거 1: 체격 공변량의 정보량 자체를 정량화 ------------------------------
d_bs <- .cc_subset(an, "score3", c(COV_FULL, "Height", "Weight"))
d_bs$.y <- d_bs$score3
m_bs <- lme4::lmer(stats::as.formula(paste0(
  ".y ~ female + ", paste(c(COV_FULL, "Height", "Weight"), collapse = " + "),
  " + (1 | center_anon)")), data = d_bs, REML = TRUE)
s_bs <- summary(m_bs)$coefficients
m_nb <- lme4::lmer(stats::as.formula(paste0(
  ".y ~ female + ", paste(COV_FULL, collapse = " + "), " + (1 | center_anon)")),
  data = d_bs, REML = TRUE)
s_nb <- summary(m_nb)$coefficients

b_fem_nb <- s_nb["female", "Estimate"]; se_fem_nb <- s_nb["female", "Std. Error"]
b_fem_bs <- s_bs["female", "Estimate"]; se_fem_bs <- s_bs["female", "Std. Error"]
b_h <- s_bs["Height", "Estimate"]; se_h <- s_bs["Height", "Std. Error"]
b_w <- s_bs["Weight", "Estimate"]; se_w <- s_bs["Weight", "Std. Error"]
dh <- mean(d_bs$Height[d_bs$female == 1]) - mean(d_bs$Height[d_bs$female == 0])
dw <- mean(d_bs$Weight[d_bs$female == 1]) - mean(d_bs$Weight[d_bs$female == 0])
contrib <- b_h * dh + b_w * dw

logmsg(paste0("[A18-6] body size: height coefficient %.5f per cm (SE %.5f, P = %s), weight ",
              "coefficient %.5f per kg (SE %.5f, P = %s); mean female-male difference %.1f cm / %.1f kg; ",
              "the score difference attributable to measured body size is %.4f points, against an ",
              "adjusted sex difference of %.3f"),
       b_h, se_h, fmt_p(2 * stats::pnorm(-abs(b_h / se_h))), b_w, se_w,
       fmt_p(2 * stats::pnorm(-abs(b_w / se_w))), dh, dw, contrib, b_fem_bs)

addrow(new_row("A18-6", "Three-component procedural intensity score (0-3)", "body size",
               mclass("gaussian"), paste0(LAB_FULL, " + height + weight"), d_bs,
               "partial association of height with the score", "change in score per cm of height",
               b_h, b_h - 1.96 * se_h, b_h + 1.96 * se_h,
               2 * stats::pnorm(-abs(b_h / se_h)), "Wald z", 1,
               est_fmt = sprintf("%.5f (%.5f to %.5f)", b_h, b_h - 1.96 * se_h, b_h + 1.96 * se_h),
               note = paste0("quantifies how much information measured body size carries about ",
                             "procedural intensity once the other covariates are in the model")))
addrow(new_row("A18-6", "Three-component procedural intensity score (0-3)", "body size",
               mclass("gaussian"), paste0(LAB_FULL, " + height + weight"), d_bs,
               "partial association of weight with the score", "change in score per kg of weight",
               b_w, b_w - 1.96 * se_w, b_w + 1.96 * se_w,
               2 * stats::pnorm(-abs(b_w / se_w)), "Wald z", 1,
               est_fmt = sprintf("%.5f (%.5f to %.5f)", b_w, b_w - 1.96 * se_w, b_w + 1.96 * se_w),
               note = ""))
addrow(new_row("A18-6", "Three-component procedural intensity score (0-3)", "body size",
               "Arithmetic from the fitted model", paste0(LAB_FULL, " + height + weight"), d_bs,
               "score difference attributable to measured body size",
               "score points", contrib, NA_real_, NA_real_, NA_real_, "derived", NA_real_,
               est_fmt = sprintf("%.4f points (height %.1f cm x %.5f + weight %.1f kg x %.5f)",
                                 contrib, dh, b_h, dw, b_w),
               note = sprintf(paste0("Adjusted sex difference without body size = %.3f (SE %.4f); ",
                                     "with body size = %.3f (SE %.4f); attenuation = %.1f%%; SE ",
                                     "inflation from adding two variables that are almost collinear ",
                                     "with sex = %.0f%%. The contribution of measured body size has ",
                                     "the OPPOSITE sign to the sex difference, so the absence of ",
                                     "attenuation is what this model implies rather than evidence ",
                                     "that body size is irrelevant. Multiplying the confidence limits ",
                                     "of the height coefficient by the mean height difference of ",
                                     "%.1f cm gives a height-attributable score difference anywhere ",
                                     "between %.3f and %.3f points, a range that includes values ",
                                     "several times larger than the observed sex difference; height ",
                                     "and sex cannot be separated in this cohort (height |SMD| = %.2f)."),
                              b_fem_nb, se_fem_nb, b_fem_bs, se_fem_bs,
                              100 * (1 - b_fem_bs / b_fem_nb), 100 * (se_fem_bs / se_fem_nb - 1),
                              dh, min((b_h - 1.96 * se_h) * dh, (b_h + 1.96 * se_h) * dh),
                              max((b_h - 1.96 * se_h) * dh, (b_h + 1.96 * se_h) * dh),
                              abs(smd_cont_h <- (mean(d_bs$Height[d_bs$female == 1]) -
                                                   mean(d_bs$Height[d_bs$female == 0])) /
                                    sqrt((stats::var(d_bs$Height[d_bs$female == 1]) +
                                            stats::var(d_bs$Height[d_bs$female == 0])) / 2)))))

# --- 근거 2: 신장 공통지지 구간 제한의 실제 효과 -----------------------------
h <- co$Height
CS_STRICT <- c(unname(stats::quantile(h[co$female == 0], .05,  na.rm = TRUE)),
               unname(stats::quantile(h[co$female == 1], .95,  na.rm = TRUE)))
CS_WIDE   <- c(unname(stats::quantile(h[co$female == 0], .025, na.rm = TRUE)),
               unname(stats::quantile(h[co$female == 1], .975, na.rm = TRUE)))
in_cs <- function(d, rng) !is.na(d$Height) & d$Height >= rng[1] & d$Height <= rng[2]
smd_cont <- function(x, g) {
  x1 <- x[g == 1]; x0 <- x[g == 0]
  (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) /
    sqrt((stats::var(x1, na.rm = TRUE) + stats::var(x0, na.rm = TRUE)) / 2)
}

for (cs in list(list(tag = "no restriction", rng = c(-Inf, Inf)),
                list(tag = sprintf("common support, strict (%.0f-%.0f cm)", CS_STRICT[1], CS_STRICT[2]),
                     rng = CS_STRICT),
                list(tag = sprintf("common support, wide (%.0f-%.0f cm)", CS_WIDE[1], CS_WIDE[2]),
                     rng = CS_WIDE))) {
  dsub <- an[in_cs(an, cs$rng) | is.infinite(cs$rng[1]), , drop = FALSE]
  if (is.infinite(cs$rng[1])) dsub <- an[!is.na(an$Height), , drop = FALSE]
  sb <- fit_stratum(dsub, "Bilati_inj", COV_FULL, "binomial")
  smd_h <- smd_cont(sb$d$Height, sb$d$female)
  frac_f <- sum(dsub$female == 1) / sum(an$female == 1)
  frac_m <- sum(dsub$female == 0) / sum(an$female == 0)
  logmsg(paste0("[A18-6] bilateral, %s: n = %d (M %d / F %d; %.0f%% of men and %.0f%% of women ",
                "retained), OR = %.2f (%.2f-%.2f), P = %s, residual height SMD = %.2f"),
         cs$tag, sb$n, sb$n_male, sb$n_female, 100 * frac_m, 100 * frac_f,
         exp(sb$b), exp(sb$lo), exp(sb$hi), fmt_p(sb$p), smd_h)
  addrow(new_row("A18-6", "Bilateral contrast injection", "height common support",
                 mclass("binomial"), LAB_FULL, NULL,
                 paste0("female effect, ", cs$tag), "OR (female vs male)",
                 exp(sb$b), exp(sb$lo), exp(sb$hi), sb$p, "Wald z", 1,
                 est_fmt = sprintf("%.2f (%.2f-%.2f)", exp(sb$b), exp(sb$lo), exp(sb$hi)),
                 note = sprintf(paste0("residual height SMD within the window = %.2f (a value above ",
                                       "0.25 indicates the window does not achieve balance); %.0f%% of ",
                                       "men and %.0f%% of women retained; 95%% CI width on the log-odds ",
                                       "scale = %.2f"),
                                smd_h, 100 * frac_m, 100 * frac_f, sb$hi - sb$lo)))
  res[[length(res)]]$n <- sb$n
  res[[length(res)]]$n_male <- sb$n_male
  res[[length(res)]]$n_female <- sb$n_female
}

###############################################################################
# 저장
###############################################################################

tab <- dplyr::bind_rows(res)
tab$decision_rule <- DECISION_RULE
tab$caveat <- CAVEAT
write.csv(tab, file.path(DIR_TAB, "table21_interaction_tests.csv"), row.names = FALSE, na = "")
logmsg("[A18] wrote %d rows to table21_interaction_tests.csv", nrow(tab))

# --- docx (핵심 요약만) ------------------------------------------------------
key <- tab[grepl("^INTERACTION", tab$quantity) | tab$analysis == "A18-0", ]
key_out <- data.frame(
  Analysis = key$analysis, Outcome = key$outcome, Modifier = key$modifier,
  Model = key$model_class, N = key$n, Quantity = key$quantity,
  Estimate = key$est_fmt, P = key$p_fmt, Test = key$test, stringsAsFactors = FALSE)
ft <- flextable::flextable(key_out)
ft <- flextable::autofit(flextable::fontsize(ft, size = 7, part = "all"))
ft <- flextable::add_footer_lines(ft, DECISION_RULE)
ft <- flextable::add_footer_lines(ft, CAVEAT)
doc <- officer::read_docx()
doc <- officer::body_add_par(doc, "Table 21. Tests of interaction between sex and lesion complexity",
                             style = "heading 1")
doc <- flextable::body_add_flextable(doc, ft)
print(doc, target = file.path(DIR_TAB, "table21_interaction_tests.docx"))

writeLines(LOG, file.path(DIR_REPX, "interaction_09_log.txt"))
capture.output(utils::sessionInfo(), file = file.path(DIR_REPX, "interaction_09_sessionInfo.txt"))
message("== 09_interaction_tests.R done ==")
