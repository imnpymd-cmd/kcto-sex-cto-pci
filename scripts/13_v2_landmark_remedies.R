# =============================================================================
# 13_v2_landmark_remedies.R — v2 검토(3): "MACCE + 1·3년 landmark" 한계 극복 검증
#
# 연구자 질의(2026-08-14): "남는 한계를 극복할 수 있는 방안이 있나?"
# 검증 3종 (MACCE_pub 중심, 창 경계 365/1095 = scripts/12 와 동일):
#   (R1) EPP 부족(>3년 창 Model4 EPP 6.8) → 창내 IPTW(성향점수는 공변량 26개로 추정하되
#        outcome 모형은 female 1개 모수) 로 우회 — 절약형·Model4 와 수렴하는지 확인
#   (R2) 창별 PH 진단 — 각 창 내부에서 cox.zph (절약형 모형) : 창 분할이 PH 문제를
#        실제로 해소했는지 형식 확인
#   (R3) 창별 첫 사건 성분 분해 — 각 창의 MACCE 첫 사건이 어느 성분(사망/뇌졸중/MI/재관류)
#        인지 성별로 집계 : 재관류 혼입(순환논리) 논쟁의 실증 근거
# 출력: interim_report/v2_nace_primary/table23_landmark_remedies.csv (+ 콘솔 요약)
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(survival) })
set.seed(42)

PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "00_cohort_prep.R"))
DIR_V2 <- file.path(PROJ, "output", "tables")   # v2 채택(2026-08-14)으로 정식 표 위치에 출력

co <- build_cohort(); stopifnot(nrow(co) == 6777)
an <- co[stats::complete.cases(co[, COV_M4, drop = FALSE]), , drop = FALSE]

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
WINDOWS <- list(
  list(start = 0,    end = 365,  lab = "0-1 year"),
  list(start = 365,  end = 1095, lab = "1-3 years (landmark 1y)"),
  list(start = 1095, end = Inf,  lab = ">3 years (landmark 3y)")
)

out_rows <- list()
add <- function(...) out_rows[[length(out_rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

## ---- R1: 창내 IPTW (stabilized ATE) ----------------------------------------
cat("== R1: window IPTW (MACCE_pub) ==\n")
for (w in WINDOWS) {
  d <- window_data("MACCE_pub", w$start, w$end)
  d <- droplevels(d[stats::complete.cases(d[, COV_M4, drop = FALSE]), , drop = FALSE])
  ps <- stats::glm(stats::as.formula(paste("female ~", paste(COV_M4, collapse = " + "))),
                   data = d, family = stats::binomial())
  p  <- stats::fitted(ps); pf <- mean(d$female)
  d$w_iptw <- ifelse(d$female == 1, pf / p, (1 - pf) / (1 - p))
  # 극단 가중치 1/99 백분위 절단 (v1 IPTW 민감도와 동일 관행)
  qs <- stats::quantile(d$w_iptw, c(.01, .99)); d$w_trim <- pmin(pmax(d$w_iptw, qs[1]), qs[2])
  for (wt in c("w_iptw", "w_trim")) {
    m <- survival::coxph(Surv(.time, .event) ~ female, data = d, weights = d[[wt]],
                         cluster = d$center_anon, robust = TRUE, ties = "efron")
    s <- summary(m)
    add(analysis = "R1-IPTW", window = w$lab, model = ifelse(wt == "w_iptw", "IPTW", "IPTW trimmed 1/99"),
        n = nrow(d), n_events = sum(d$.event), n_parameters = 1,
        HR_fmt = sprintf("%.2f (%.2f-%.2f)", s$conf.int[1,1], s$conf.int[1,3], s$conf.int[1,4]),
        p_fmt = fmt_p(s$coefficients[1, "Pr(>|z|)"]),
        note = sprintf("PS model events(female)=%d, PS EPP=%.1f; weight max=%.1f",
                       sum(d$female), sum(d$female)/length(coef(ps)), max(d$w_iptw)))
    cat(sprintf("  %-26s %-18s HR %s p %s\n", w$lab, ifelse(wt=="w_iptw","IPTW","IPTW trim"),
                tail(out_rows,1)[[1]]$HR_fmt, tail(out_rows,1)[[1]]$p_fmt))
  }
}

## ---- R2: 창별 PH 진단 (절약형 모형) -----------------------------------------
cat("== R2: within-window PH (parsimonious) ==\n")
for (w in WINDOWS) {
  d <- window_data("MACCE_pub", w$start, w$end)
  d <- droplevels(d[stats::complete.cases(d[, COV_PARS, drop = FALSE]), , drop = FALSE])
  m <- cox_model(d, COV_PARS)
  z <- survival::cox.zph(m)
  pf <- z$table["female", "p"]; pg <- z$table["GLOBAL", "p"]
  add(analysis = "R2-PH", window = w$lab, model = "parsimonious cox.zph",
      n = nrow(d), n_events = sum(d$.event), n_parameters = length(coef(m)),
      HR_fmt = "", p_fmt = "",
      note = sprintf("zph female p=%.3f; global p=%.3f", pf, pg))
  cat(sprintf("  %-26s zph female p=%.3f / global p=%.3f\n", w$lab, pf, pg))
}

## ---- R3: 창별 첫 사건 성분 분해 ---------------------------------------------
## 첫 MACCE 사건 시점(.time)과 성분별 사건시각을 대조. 동시(tie)는 위계
## 사망 > 뇌졸중 > MI > 재관류 로 1건 귀속 (관례적 위계, 각주 필요).
cat("== R3: first-event composition by window/sex ==\n")
comp_def <- list(death  = c("Death", "death_fudur"),
                 stroke = c("CVA_re", "CVA_fudur"),
                 mi     = c("MI", "MI_fudur"),
                 revasc = NULL)  # 재관류 = 나머지 (TVR/NLR/CABG 중 최솟값)
rv_t <- pmin(suppressWarnings(as.numeric(an$TVR_fudur2)),
             suppressWarnings(as.numeric(an$NLR_fudur2)),
             suppressWarnings(as.numeric(an$CABG_fudur)), na.rm = TRUE)
for (w in WINDOWS) {
  d <- window_data("MACCE_pub", w$start, w$end)
  ev <- d[d$.event == 1, , drop = FALSE]
  t_abs <- ev$.time + w$start                       # 원 시간축 사건 시각
  # day-0 규칙(0.5일)을 성분 시각과 대조하기 위해 반내림 허용 오차 1일
  cls <- function(row_i, t) {
    r <- ev[row_i, ]
    hit <- function(evv, tmv) !is.na(r[[evv]]) && r[[evv]] == 1 &&
      !is.na(r[[tmv]]) && abs(as.numeric(r[[tmv]]) - t) <= 1
    if (hit("Death", "death_fudur"))  return("death")
    if (hit("CVA_re", "CVA_fudur"))   return("stroke")
    if (hit("MI", "MI_fudur"))        return("mi")
    "revasc"
  }
  typ <- vapply(seq_len(nrow(ev)), function(i) cls(i, t_abs[i]), character(1))
  tb <- table(sex = ifelse(ev$female == 1, "Women", "Men"), component = typ)
  pr <- prop.table(tb, 1) * 100
  for (sx in rownames(tb)) {
    add(analysis = "R3-comp", window = w$lab, model = sx,
        n = sum(tb[sx, ]), n_events = sum(tb[sx, ]), n_parameters = NA,
        HR_fmt = "", p_fmt = "",
        note = paste(sprintf("%s %d (%.0f%%)", colnames(tb), tb[sx, ], pr[sx, ]), collapse = "; "))
    cat(sprintf("  %-26s %-6s %s\n", w$lab, sx, tail(out_rows,1)[[1]]$note))
  }
}

tab <- dplyr::bind_rows(out_rows)
tab$caveat <- "Exploratory v2 analysis (frozen v1 outputs untouched); not prospectively registered."
write.csv(tab, file.path(DIR_V2, "table23_landmark_remedies.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("saved: interim_report/v2_nace_primary/table23_landmark_remedies.csv\n")
