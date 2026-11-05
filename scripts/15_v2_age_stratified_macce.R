# =============================================================================
# 15_v2_age_stratified_macce.R — v2 Figure 2 원천표: 연령대별 MACCE
#
# 연구자 결정(2026-08-14): Figure 2 를 DMS 기반(table6f)에서 MACCE 기반으로 교체.
# 산출: 연령대(<65 / 65-74 / >=75) × 성별
#   - 5년 KM 누적발생률 (95% CI, log-log)
#   - 연령대내 미보정 HR (여 vs 남, 센터 clustered SE) + 대내 PH 진단(zph)
# 출력: output/tables/table26_age_stratified_MACCE.csv
# 주의: MACCE 는 전기간 PH 위반이므로 대내 zph p 를 반드시 함께 인쇄하고,
#       figure legend 에 대내 HR 이 전기간 요약임을 명시할 것.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(survival) })
set.seed(42)
PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "00_cohort_prep.R"))

co <- build_cohort(); stopifnot(nrow(co) == 6777)
d0 <- make_surv(co, "MACCE_pub")
d0$band <- cut(d0$Age, breaks = c(-Inf, 64.999, 74.999, Inf),
               labels = c("<65 years", "65-74 years", ">=75 years"))
stopifnot(!any(is.na(d0$band)))

rows <- list()
for (b in levels(d0$band)) {
  d <- d0[d0$band == b, , drop = FALSE]
  f <- survfit(Surv(.time/365.25, .event) ~ female, data = d, conf.type = "log-log")
  sm <- summary(f, times = 5, extend = TRUE)
  m <- coxph(Surv(.time, .event) ~ female, data = d, cluster = d$center_anon,
             robust = TRUE, ties = "efron")
  s <- summary(m)
  z <- tryCatch(cox.zph(m)$table["female", "p"], error = function(e) NA_real_)
  for (j in seq_along(sm$strata)) {
    sex <- ifelse(grepl("female=1", sm$strata[j]), "Female", "Male")
    rows[[length(rows)+1]] <- data.frame(
      outcome_key = "MACCE_pub", age_stratum = b, sex = sex,
      n = sum(d$female == (sex=="Female")),
      events = sum(d$.event[d$female == (sex=="Female")]),
      cum_inc_5yr_pct = round(100*(1-sm$surv[j]), 1),
      ci_lo = round(100*(1-sm$upper[j]), 1), ci_hi = round(100*(1-sm$lower[j]), 1),
      within_stratum_HR = sprintf("%.2f (%.2f-%.2f)", s$conf.int[1,1], s$conf.int[1,3], s$conf.int[1,4]),
      within_stratum_p = sprintf("%.3f", s$coefficients[1, "Pr(>|z|)"]),
      within_stratum_zph_p = sprintf("%.3f", z),
      stringsAsFactors = FALSE)
  }
}
tab <- bind_rows(rows)
write.csv(tab, file.path(PROJ, "output", "tables", "table26_age_stratified_MACCE.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
print(tab[, c("age_stratum","sex","n","events","cum_inc_5yr_pct","ci_lo","ci_hi",
              "within_stratum_HR","within_stratum_p","within_stratum_zph_p")], row.names = FALSE)
cat("saved: output/tables/table26_age_stratified_MACCE.csv\n")
