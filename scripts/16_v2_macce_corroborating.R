# =============================================================================
# 16_v2_macce_corroborating.R — Phase 5(v2) QA V2-03(b): primary(MACCE) 보강 배터리
#
# 창별(0-1 / 1-3 / >3년, 경계 365/1095) MACCE_pub 에 대해:
#   1) 엔트로피 균형 (WeightIt ebal ATE; outcome 모형은 female 1모수)
#   2) 1:1 PSM (MatchIt nearest glm, caliper 0.2 SD)
#   3) 다중대입 m=20 (mice, Nelson-Aalen 보조; 04 스크립트 사양) — 창별 적합·pool
#   4) day-0 사건 제외 (0-1년 창만 해당) / 5년 행정 절단 (>3년 창만 해당)
#   5) 중복행 제거 2종 (05 A7 규칙 재현: CO_I / CO_II)
# 분석 모형 원칙(Table 3 과 동일): 0-1·1-3년 = Model 4 / >3년 = 절약형.
# 출력: output/tables/table27_macce_corroborating.csv → Supplement S11(b)
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(survival) })
set.seed(42)
PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "00_cohort_prep.R"))

co <- build_cohort(); stopifnot(nrow(co) == 6777)
an <- co[stats::complete.cases(co[, COV_M4, drop = FALSE]), , drop = FALSE]

window_data <- function(dat, start, end = Inf, day0 = "half_day") {
  d <- make_surv(dat, "MACCE_pub", day0 = day0)
  d <- d[d$.time > start, , drop = FALSE]
  d$.time <- d$.time - start
  if (is.finite(end)) { w <- end - start
    d$.event <- ifelse(d$.time <= w & d$.event == 1, 1L, 0L); d$.time <- pmin(d$.time, w) }
  d
}
WINS <- list(list(start = 0, end = 365, lab = "0-1 year", covs = quote(COV_M4), covlab = "Model 4"),
             list(start = 365, end = 1095, lab = "1-3 years", covs = quote(COV_M4), covlab = "Model 4"),
             list(start = 1095, end = Inf, lab = ">3 years", covs = quote(COV_PARS), covlab = "parsimonious"))
cc <- function(d, covs) {
  v <- intersect(unique(c(".time",".event","female","center_anon",covs)), names(d))
  droplevels(d[stats::complete.cases(d[, v, drop = FALSE]), , drop = FALSE])
}
rows <- list()
addr <- function(win, method, m, note = "") {
  s <- summary(m); i <- which(rownames(s$coefficients) == "female")
  se_col <- if ("robust se" %in% colnames(s$coefficients)) "robust se" else "se(coef)"
  rows[[length(rows)+1]] <<- data.frame(
    window = win, method = method, n = m$n, n_events = m$nevent,
    HR = unname(s$conf.int[i,1]), lo = unname(s$conf.int[i,3]), hi = unname(s$conf.int[i,4]),
    HR_fmt = sprintf("%.2f (%.2f-%.2f)", s$conf.int[i,1], s$conf.int[i,3], s$conf.int[i,4]),
    p_value = unname(s$coefficients[i, "Pr(>|z|)"]),
    p_fmt = fmt_p(unname(s$coefficients[i, "Pr(>|z|)"])), note = note, stringsAsFactors = FALSE)
}

## ---- A7 중복 플래그 재현 (05_qa_sensitivity.R 과 동일 규칙; PHI 출력 금지) ----
dup_flags <- local({
  raw <- readRDS(file.path(PROJ, "data", "processed", "kcto_sex_analysis.rds"))
  raw$index_yr <- suppressWarnings(as.numeric(as.character(raw$index_yr)))
  sel <- raw$flag_cto_pci %in% TRUE & !is.na(raw$index_yr) & raw$index_yr >= 2014 & raw$index_yr <= 2023
  stopifnot(sum(sel) == nrow(co))
  sub <- raw[sel, , drop = FALSE]
  key <- paste(sub$id_code, sub$id, sep = "|")
  dup_entry <- duplicated(paste(key, as.character(sub$CTO_intervention_date), sub$Age, sub$Sex, sep = "||"))
  drop_repeat <- rep(FALSE, nrow(sub))
  keep1 <- which(!dup_entry)
  for (g in split(keep1, key[keep1])) {
    if (length(g) < 2) next
    sx <- sub$Sex[g]; ag <- suppressWarnings(as.numeric(sub$Age[g])); yr <- sub$index_yr[g]
    if (length(unique(sx[!is.na(sx)])) <= 1 && !any(is.na(ag)) && !any(is.na(yr)) &&
        abs((max(ag) - min(ag)) - (max(yr) - min(yr))) <= 1)
      drop_repeat[g[order(sub$index_yr[g])][-1]] <- TRUE
  }
  list(dup_entry = dup_entry, drop_repeat = drop_repeat)
})
CO_I  <- co[!dup_flags$dup_entry, , drop = FALSE]
CO_II <- co[!(dup_flags$dup_entry | dup_flags$drop_repeat), , drop = FALSE]
AN_I  <- CO_I [stats::complete.cases(CO_I [, COV_M4, drop = FALSE]), , drop = FALSE]
AN_II <- CO_II[stats::complete.cases(CO_II[, COV_M4, drop = FALSE]), , drop = FALSE]

for (w in WINS) {
  covs <- eval(w$covs)
  base <- cc(window_data(an, w$start, w$end), covs)
  message(sprintf("== window %s (n=%d, events=%d) ==", w$lab, nrow(base), sum(base$.event)))

  ## 1) 엔트로피 균형
  ps_form <- stats::as.formula(paste("female ~", paste(COV_M4, collapse = " + ")))
  eb <- tryCatch(WeightIt::weightit(ps_form, data = base, method = "ebal", estimand = "ATE"),
                 error = function(e) NULL)
  if (!is.null(eb)) {
    base$w_eb <- eb$weights
    m <- coxph(Surv(.time, .event) ~ female, data = base, weights = base$w_eb,
               cluster = base$center_anon, robust = TRUE, ties = "efron")
    addr(w$lab, "Entropy balancing (ATE)", m, sprintf("max weight %.1f", max(base$w_eb)))
  }

  ## 2) 1:1 PSM
  M <- tryCatch(MatchIt::matchit(ps_form, data = base, method = "nearest",
                                 distance = "glm", caliper = 0.2),
                error = function(e) NULL)
  if (!is.null(M)) {
    md <- MatchIt::match.data(M)
    m <- coxph(Surv(.time, .event) ~ female, data = md,
               cluster = md$center_anon, robust = TRUE, ties = "efron")
    addr(w$lab, "1:1 propensity-score matching (caliper 0.2 SD)", m,
         sprintf("matched pairs %d", sum(md$female == 1)))
  }

  ## 4) day-0 제외 / 5년 절단 — 해당 창만
  if (w$start == 0) {
    d0 <- cc(window_data(an, w$start, w$end, day0 = "exclude"), covs)
    addr(w$lab, "Day-0 events excluded", cox_model(d0, covs), w$covlab)
  }
  if (!is.finite(w$end)) {
    d5 <- base; d5$.event <- ifelse(d5$.time <= (1825 - w$start) & d5$.event == 1, 1L, 0L)
    d5$.time <- pmin(d5$.time, 1825 - w$start)
    addr(w$lab, "Follow-up censored at 5 years", cox_model(d5, covs), w$covlab)
  }

  ## 5) 중복 제거 2종
  for (cfg in list(list(d = AN_I,  tag = sprintf("Duplicate data-entry rows removed (N=%d)", nrow(CO_I))),
                   list(d = AN_II, tag = sprintf("Duplicates removed + index retained (N=%d)", nrow(CO_II))))) {
    dd <- cc(window_data(cfg$d, w$start, w$end), covs)
    addr(w$lab, cfg$tag, cox_model(dd, covs), w$covlab)
  }

  ## 3) 다중대입 m=20 (분석모형: 창 원칙 그대로 / 대입모형: 전 공변량 + NA 위험)
  mi_vars <- unique(c(COV_M4, "preEF"))
  full <- window_data(co, w$start, w$end)
  mi_df <- full[, c(".time", ".event", "female", "center_anon", mi_vars)]
  na_est <- survival::basehaz(coxph(Surv(.time, .event) ~ 1, data = mi_df), centered = FALSE)
  mi_df$.nelsonaalen <- stats::approx(na_est$time, na_est$hazard, xout = mi_df$.time,
                                      method = "constant", rule = 2)$y
  pred <- mice::make.predictorMatrix(mi_df)
  pred[, ".time"] <- 0; pred[".time", ] <- 0; pred[".event", ] <- 0; pred[".nelsonaalen", ] <- 0
  pred[, "center_anon"] <- 0; pred["center_anon", ] <- 0
  meth <- mice::make.method(mi_df); meth[c(".time", ".event", ".nelsonaalen", "center_anon")] <- ""
  imp <- mice::mice(mi_df, m = 20, maxit = 10, printFlag = FALSE, seed = 42,
                    predictorMatrix = pred, method = meth)
  f_an <- stats::as.formula(paste("Surv(.time, .event) ~ female +", paste(covs, collapse = " + ")))
  fits <- lapply(seq_len(imp$m), function(i) coxph(f_an, data = mice::complete(imp, i)))
  pl <- summary(mice::pool(mice::as.mira(fits)), conf.int = TRUE)
  pr <- pl[pl$term == "female", ]
  rows[[length(rows)+1]] <- data.frame(
    window = w$lab, method = "Multiple imputation (m=20)", n = nrow(mi_df),
    n_events = sum(mi_df$.event),
    HR = exp(pr$estimate), lo = exp(pr$`2.5 %`), hi = exp(pr$`97.5 %`),
    HR_fmt = sprintf("%.2f (%.2f-%.2f)", exp(pr$estimate), exp(pr$`2.5 %`), exp(pr$`97.5 %`)),
    p_value = pr$p.value, p_fmt = fmt_p(pr$p.value),
    note = paste0(w$covlab, "; SE not clustered (mice pooling)"), stringsAsFactors = FALSE)
  message(sprintf("  MI done: %s", tail(rows,1)[[1]]$HR_fmt))
}

tab <- dplyr::bind_rows(rows)
tab$caveat <- paste0("Analysis model follows the window principle of Table 3 (Model 4 in the first two ",
                     "windows, parsimonious beyond 3 years); balancing/matching estimate the propensity ",
                     "score on the full covariate set with a single-parameter outcome model. Not prospectively registered.")
write.csv(tab, file.path(PROJ, "output", "tables", "table27_macce_corroborating.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
print(tab[, c("window","method","n","n_events","HR_fmt","p_fmt")], row.names = FALSE)
cat("saved: output/tables/table27_macce_corroborating.csv\n")
