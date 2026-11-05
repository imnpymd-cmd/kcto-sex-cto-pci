# =============================================================================
# 12_v2_macce_landmark.R — v2 검토용: landmark 분할 프레임워크 (연구자 질의 2026-08-14)
#
# 질의: "MACCE + 1년, 3년 landmark analysis를 하면?"
# 구조: 추적을 0–1년 / 1–3년(1년 landmark) / >3년(3년 landmark) 세 창으로 분할하고
#       각 창에서 모형을 **재적합**한다 (table8b 의 구간별 HR 은 공변량 효과를 전 구간
#       공유하는 time-varying 모형; 여기서는 A17(table19b)과 같은 창내 재적합).
# 경계: 365 / 1095일 — 04_primary_analysis.R 의 구간 절단(c(365,1095)) 및 A17(LM=1095)과 동일.
# 모형: 미보정 / 연령만 / 절약형(age·DM·CKD·AF) / Model 4 — A17 사양 그대로. EPP 병기.
# 출력: interim_report/v2_nace_primary/table22_landmark_partition.csv
#       (v1 산출물 동결 원칙에 따라 output/tables 가 아니라 v2 폴더에 저장)
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(survival) })
set.seed(42)

PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "00_cohort_prep.R"))
DIR_V2 <- file.path(PROJ, "output", "tables")   # v2 채택(2026-08-14)으로 정식 표 위치에 출력

co <- build_cohort(); stopifnot(nrow(co) == 6777)
an <- co[stats::complete.cases(co[, COV_M4, drop = FALSE]), , drop = FALSE]
message(sprintf("[analysis set] N = %d (F = %d)", nrow(an), sum(an$female)))

cc_subset <- function(dat, covs) {
  vars <- intersect(unique(c(".time", ".event", "female", "center_anon", covs)), names(dat))
  droplevels(dat[stats::complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE])
}

# 창 정의: start 시점 event-free & 관찰 중인 환자만 남기고, 시간 재설정 후 end 에서 절단
window_data <- function(key, start, end = Inf) {
  d <- make_surv(an, key)
  d <- d[d$.time > start, , drop = FALSE]
  d$.time <- d$.time - start
  if (is.finite(end)) {
    w <- end - start
    d$.event <- ifelse(d$.time <= w & d$.event == 1, 1L, 0L)
    d$.time  <- pmin(d$.time, w)
  }
  d
}

SPECS <- list(
  list(covs = character(0), lab = "unadjusted"),
  list(covs = COV_AGE,      lab = "age only"),
  list(covs = COV_PARS,     lab = "parsimonious (age, DM, CKD, AF)"),
  list(covs = COV_M4,       lab = "Model 4 (re-fitted within window)")
)
WINDOWS <- list(
  list(start = 0,    end = 365,  lab = "0-1 year (all patients, censored at 1 year)"),
  list(start = 365,  end = 1095, lab = "1-3 years (landmark at 1 year)"),
  list(start = 1095, end = Inf,  lab = ">3 years (landmark at 3 years)")
)

rows <- list()
for (key in c("MACCE_pub", "DMS", "NACE_pub", "MACE_pub", "Death", "AnyRevasc")) {
  for (w in WINDOWS) {
    dw <- window_data(key, w$start, w$end)
    for (sp in SPECS) {
      d <- cc_subset(dw, sp$covs)
      m <- tryCatch(cox_model(d, sp$covs), error = function(e) NULL)
      if (is.null(m)) next
      s <- summary(m); i <- which(rownames(s$coefficients) == "female")
      npar <- length(stats::coef(m))
      rows[[length(rows) + 1]] <- data.frame(
        analysis = "V2-LM", outcome = OUT[[key]]$lab, outcome_key = key,
        window = w$lab, model = sp$lab,
        n = m$n, n_events = m$nevent,
        n_female = sum(d$female == 1), events_female = sum(d$.event[d$female == 1] == 1),
        n_parameters = npar, events_per_parameter = round(m$nevent / npar, 2),
        HR = unname(s$conf.int[i, 1]), lo = unname(s$conf.int[i, 3]), hi = unname(s$conf.int[i, 4]),
        HR_fmt = sprintf("%.2f (%.2f-%.2f)", s$conf.int[i, 1], s$conf.int[i, 3], s$conf.int[i, 4]),
        p_value = unname(s$coefficients[i, "Pr(>|z|)"]),
        p_fmt = fmt_p(unname(s$coefficients[i, "Pr(>|z|)"])),
        stringsAsFactors = FALSE)
    }
  }
}
tab <- dplyr::bind_rows(rows)
tab$caveat <- paste0(
  "Landmark windows condition on being event-free and under observation at the window start ",
  "(depletion of susceptibles); estimates are associations, not effects. Exploratory v2 analysis; ",
  "not part of the frozen v1 output set. Not prospectively registered.")
write.csv(tab, file.path(DIR_V2, "table22_landmark_partition.csv"), row.names = FALSE, fileEncoding = "UTF-8")

for (k in unique(tab$outcome_key)) {
  cat("\n==", k, "==\n")
  print(tab[tab$outcome_key == k, c("window", "model", "n", "n_events", "events_per_parameter", "HR_fmt", "p_fmt")],
        row.names = FALSE)
}
cat("\nsaved: interim_report/v2_nace_primary/table22_landmark_partition.csv\n")
