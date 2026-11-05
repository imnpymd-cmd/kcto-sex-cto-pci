# =============================================================================
# 14_v2_landmark_extras.R — v2 Phase 2 보완: §1-2v2 항목 9
#   (E1) landmark 시점 이동 민감도 — 후기 창 시작을 2년(730)·4년(1460)으로 이동, 절약형
#   (E2) 창별 절대위험 — 창 위험집합 내 KM 누적발생률(창 종료 시점), 성별
# 출력: interim_report/v2_nace_primary/table24_landmark_shift.csv
#       interim_report/v2_nace_primary/table25_window_absolute_risk.csv
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(survival) })
set.seed(42)
PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "00_cohort_prep.R"))
DIR_V2 <- file.path(PROJ, "output", "tables")   # v2 채택(2026-08-14)으로 정식 표 위치에 출력

co <- build_cohort(); stopifnot(nrow(co) == 6777)
an <- co[stats::complete.cases(co[, COV_M4, drop = FALSE]), , drop = FALSE]

window_data <- function(key, start, end = Inf) {
  d <- make_surv(an, key); d <- d[d$.time > start, , drop = FALSE]
  d$.time <- d$.time - start
  if (is.finite(end)) { w <- end - start
    d$.event <- ifelse(d$.time <= w & d$.event == 1, 1L, 0L); d$.time <- pmin(d$.time, w) }
  d
}
cc <- function(dat, covs) {
  vars <- intersect(unique(c(".time",".event","female","center_anon",covs)), names(dat))
  droplevels(dat[stats::complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE])
}

## E1: 후기 landmark 이동 (MACCE, 절약형 + 연령만)
cat("== E1: landmark shift (MACCE, late window) ==\n")
r1 <- list()
for (lm in c(730, 1095, 1460)) {
  d0 <- window_data("MACCE_pub", lm, Inf)
  for (sp in list(list(covs = COV_AGE, lab = "age only"),
                  list(covs = COV_PARS, lab = "parsimonious"))) {
    d <- cc(d0, sp$covs); m <- cox_model(d, sp$covs); s <- summary(m)
    i <- which(rownames(s$coefficients) == "female")
    r1[[length(r1)+1]] <- data.frame(
      landmark_days = lm, landmark_lab = sprintf("%.0f years", lm/365),
      model = sp$lab, n = m$n, n_events = m$nevent,
      HR_fmt = sprintf("%.2f (%.2f-%.2f)", s$conf.int[i,1], s$conf.int[i,3], s$conf.int[i,4]),
      p_fmt = fmt_p(s$coefficients[i, "Pr(>|z|)"]), stringsAsFactors = FALSE)
    cat(sprintf("  landmark %4.0fd %-12s HR %s p %s (events %d)\n", lm, sp$lab,
                tail(r1,1)[[1]]$HR_fmt, tail(r1,1)[[1]]$p_fmt, m$nevent))
  }
}
write.csv(bind_rows(r1), file.path(DIR_V2, "table24_landmark_shift.csv"), row.names = FALSE)

## E2: 창별 절대위험 (KM, 창 종료 시점; >3년 창은 landmark 후 2년 = 원 시간축 5년)
cat("== E2: window absolute risk (KM) ==\n")
WINDOWS <- list(list(start=0, end=365, at=365, lab="0-1 year"),
                list(start=365, end=1095, at=730, lab="1-3 years (landmark 1y)"),
                list(start=1095, end=Inf, at=730, lab=">3 years (landmark 3y, incidence at +2y)"))
r2 <- list()
for (key in c("MACCE_pub", "DMS")) {
  for (w in WINDOWS) {
    d <- window_data(key, w$start, w$end)
    f <- survfit(Surv(.time, .event) ~ female, data = d, conf.type = "log-log")
    sm <- summary(f, times = w$at, extend = TRUE)
    for (j in seq_along(sm$strata)) {
      sex <- ifelse(grepl("female=1", sm$strata[j]), "Women", "Men")
      r2[[length(r2)+1]] <- data.frame(
        outcome_key = key, window = w$lab, sex = sex,
        n_at_window_start = sum(d$female == (sex=="Women")),
        events_in_window = sum(d$.event[d$female == (sex=="Women")]),
        cuminc_pct = round(100*(1 - sm$surv[j]), 1),
        ci = sprintf("%.1f-%.1f", 100*(1-sm$upper[j]), 100*(1-sm$lower[j])),
        stringsAsFactors = FALSE)
      cat(sprintf("  %-10s %-38s %-6s %5.1f%% (%s) [n=%d, ev=%d]\n", key, w$lab, sex,
                  tail(r2,1)[[1]]$cuminc_pct, tail(r2,1)[[1]]$ci,
                  tail(r2,1)[[1]]$n_at_window_start, tail(r2,1)[[1]]$events_in_window))
    }
  }
}
write.csv(bind_rows(r2), file.path(DIR_V2, "table25_window_absolute_risk.csv"), row.names = FALSE)
cat("saved: table24_landmark_shift.csv / table25_window_absolute_risk.csv\n")
