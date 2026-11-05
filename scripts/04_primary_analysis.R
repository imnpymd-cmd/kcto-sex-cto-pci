###############################################################################
# K-CTO Registry 2.0 — Sex Difference Analysis
# Phase 2 (primary analysis): 보정 다변량 · 경쟁위험 · 인과추론 보강 · 민감도
#
# 실행:  Rscript scripts/04_primary_analysis.R
# 입력:  data/processed/kcto_sex_analysis.rds  (via scripts/00_cohort_prep.R)
# 출력:  output/tables/table6*..table12*, tableS8*..tableS14*
#        interim_report/phase2_stats/reports/*
#
# 코호트: flag_cto_pci & index_yr 2014-2023 -> N = 6,777 (CLAUDE.md §1)
# Primary outcome: DMS = all-cause death + MI + stroke (재관류 제외) (CLAUDE.md §1-2)
# 노출: sex (female vs male). ⚠️ 무작위 배정이 아니며 조작 가능한 개입도 아니다.
#       모든 추정치는 **연관성(association)** 이며 인과 해석 금지.
#
# 핵심 질문: 여성의 초과 사망이 연령 보정 후에도 남는가?
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
  library(splines)
})

set.seed(42)

PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "00_cohort_prep.R"))

DIR_TAB  <- file.path(PROJ, "output", "tables")
DIR_REPX <- file.path(PROJ, "interim_report", "phase2_stats", "reports")
dir.create(DIR_TAB,  recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_REPX, recursive = TRUE, showWarnings = FALSE)

co <- build_cohort()
stopifnot(TAU5 == 1825)   # 5년 (일) — 00_cohort_prep.R 정의

###############################################################################
# 1. 공변량 세트 · outcome 목록 · 모델 helper
#
#  ⚠️ 전부 scripts/00_cohort_prep.R 에 단일 정의되어 있다 (04 / 05 공용).
#     COV_AGE / COV_CLIN / COV_LES / COV_ADMIN / COV_M2..COV_M5 / COV_M4L /
#     COV_PARS / COV_PS / COV_M4_NOCENTER / TAU5 / OUT / STRAT_VARS /
#     make_surv() / fit_cox() / tidy_cox() / cox_model() / logit_sex() / crude_or()
#     여기서 재정의하지 말 것 (정의가 갈라지면 두 스크립트 결과가 어긋난다).
#
#  희소 사건(events per parameter < 10) 대응: Model 4 는 파라미터 25개 내외이므로
#  MI(76) · stroke(127) · 대출혈(48) 등에서는 과적합·분리(separation) 위험이 있어
#  COV_PARS(절약형 세트)를 병기한다.
###############################################################################

###############################################################################
# 2. 분석 데이터셋 · 공통 complete-case 집합
#    단계적 보정 비교는 **동일 표본**에서 수행해야 "연령 보정으로 얼마나
#    소실되는가"를 표본 변화와 혼동하지 않는다.
###############################################################################

ALL_MODEL_VARS <- unique(c(COV_M4, "preEF", "cto_success", "CR_yes",
                           "dc_antiplt2_f", "dc_statin",
                           "female", "sex_f", "center_anon", "index_yr",
                           "age75", "JCTO_bi", "index_half_f"))

cc_m4 <- stats::complete.cases(co[, COV_M4, drop = FALSE])
an    <- co[cc_m4, , drop = FALSE]          # 주 분석 표본 (complete case, Model 4)
message(sprintf("[analysis set] complete case for Model 4: N = %d (%.1f%% of %d); F = %d",
                nrow(an), 100 * nrow(an) / nrow(co), nrow(co), sum(an$female)))

###############################################################################
# 4. Table 6 — PRIMARY: 단계적 보정 (DMS), 동일 표본
###############################################################################

message("== Table 6: staged adjustment for the primary outcome (DMS) ==")

d_dms <- make_surv(an, "DMS")

staged <- rbind(
  fit_cox(d_dms, character(0), "Model 1: unadjusted", OUT$DMS$lab,
          "sex only"),
  fit_cox(d_dms, COV_M2, "Model 2: + age", OUT$DMS$lab,
          "age (continuous)"),
  fit_cox(d_dms, COV_M3, "Model 3: + age, clinical covariates", OUT$DMS$lab,
          "age, BMI, smoking, HTN, DM, CKD, prior PCI/CABG/MI/stroke/PAD, AF, stable angina"),
  fit_cox(d_dms, COV_M4, "Model 4: fully adjusted (PRIMARY)", OUT$DMS$lab,
          "Model 3 + J-CTO, prox/mid CTO, ISR CTO, Rentrop 3, non-CTO PCI, vessel disease, index period, center volume"),
  fit_cox(d_dms, COV_M5, "Model 5: + procedural variables (exploratory, mediator-adjusted)", OUT$DMS$lab,
          "Model 4 + technical success, complete revascularization, discharge antiplatelet regimen, statin")
)

# 각 모델을 그 모델 자체의 complete case 로 적합한 결과 (표본 변화 영향 확인)
staged_own <- rbind(
  fit_cox(make_surv(co, "DMS"), character(0), "Model 1 (own complete case)", OUT$DMS$lab),
  fit_cox(make_surv(co[stats::complete.cases(co[, COV_M2, drop = FALSE]), ], "DMS"),
          COV_M2, "Model 2 (own complete case)", OUT$DMS$lab),
  fit_cox(make_surv(co[stats::complete.cases(co[, COV_M3, drop = FALSE]), ], "DMS"),
          COV_M3, "Model 3 (own complete case)", OUT$DMS$lab),
  fit_cox(make_surv(co[stats::complete.cases(co[, COV_M4, drop = FALSE]), ], "DMS"),
          COV_M4, "Model 4 (own complete case)", OUT$DMS$lab),
  fit_cox(make_surv(co[stats::complete.cases(co[, COV_M4L, drop = FALSE]), ], "DMS"),
          COV_M4L, "Model 4L: Model 4 + LVEF (own complete case)", OUT$DMS$lab)
)

# 연령 보정으로 소실되는 초과위험의 비율
hr1 <- staged$HR[staged$model == "Model 1: unadjusted"]
hr2 <- staged$HR[staged$model == "Model 2: + age"]
hr4 <- staged$HR[staged$model == "Model 4: fully adjusted (PRIMARY)"]
attenuation <- data.frame(
  metric = c("Crude HR", "Age-adjusted HR", "Fully adjusted HR",
             "% of excess risk explained by age alone",
             "% of excess risk explained by full adjustment"),
  value = c(sprintf("%.3f", hr1), sprintf("%.3f", hr2), sprintf("%.3f", hr4),
            sprintf("%.1f%%", 100 * (hr1 - hr2) / (hr1 - 1)),
            sprintf("%.1f%%", 100 * (hr1 - hr4) / (hr1 - 1))),
  stringsAsFactors = FALSE)

write.csv(rbind(staged, staged_own),
          file.path(DIR_TAB, "table6a_primary_staged_adjustment_DMS.csv"),
          row.names = FALSE, na = "")
write.csv(attenuation, file.path(DIR_TAB, "table6b_age_attenuation_DMS.csv"),
          row.names = FALSE, na = "")

# 완전보정 모델 전체 계수표
m4_dms <- cox_model(d_dms, COV_M4)
write.csv(tidy_cox(m4_dms, "Model 4: fully adjusted", OUT$DMS$lab),
          file.path(DIR_TAB, "table6c_full_model_coefficients_DMS.csv"),
          row.names = FALSE, na = "")

# --- 연령 유연 모형 (스플라인) --------------------------------------------
# 연령 선형성 검정이 기각되므로(§5-a), 연령을 자연 3차 스플라인으로 유연하게
# 모형화한 버전을 민감도 분석으로 병기한다. "연령 보정 부족으로 인한 잔여
# 교란" 이라는 리뷰어 지적을 차단하기 위한 핵심 분석이다.
# ⛔ 사전 계획을 시사하는 표현 금지 (CLAUDE.md §1-2) — 본 분석은 사전등록되지 않았다.
COV_M2S <- c("ns(Age, 3)")
COV_M4S <- c("ns(Age, 3)", setdiff(COV_M4, "Age"))
staged_spline <- rbind(
  fit_cox(d_dms, COV_M2S, "Model 2S: + age (natural cubic spline, 3 df)", OUT$DMS$lab,
          "flexible age adjustment"),
  fit_cox(d_dms, COV_M4S, "Model 4S: fully adjusted with spline age", OUT$DMS$lab,
          "flexible age adjustment"))
write.csv(staged_spline, file.path(DIR_TAB, "table6e_spline_age_models_DMS.csv"),
          row.names = FALSE, na = "")

# --- 연령층별 미보정 5년 누적발생률 (연령이 설명력의 대부분임을 시각화) ----
age_strat <- do.call(rbind, lapply(list(
  list(lab = "<65 years",    f = function(x) x$Age <  65),
  list(lab = "65-74 years",  f = function(x) x$Age >= 65 & x$Age < 75),
  list(lab = ">=75 years",   f = function(x) x$Age >= 75)), function(s) {
    dd <- d_dms[s$f(d_dms) %in% TRUE, , drop = FALSE]
    fit <- survival::survfit(Surv(.time, .event) ~ sex_f, data = dd, conf.type = "log-log")
    sm <- summary(fit, times = TAU5, extend = TRUE)
    cx <- survival::coxph(Surv(.time, .event) ~ female, data = dd, cluster = dd$center_anon,
                          robust = TRUE)
    cs <- summary(cx)
    do.call(rbind, lapply(c("Male", "Female"), function(g) {
      i <- which(as.character(sm$strata) == paste0("sex_f=", g))
      data.frame(age_stratum = s$lab, sex = g,
                 n = sum(dd$sex_f == g), events = sum(dd$.event[dd$sex_f == g] == 1),
                 median_age = round(stats::median(dd$Age[dd$sex_f == g], na.rm = TRUE), 1),
                 cum_inc_5yr_pct = round(100 * (1 - sm$surv[i]), 1),
                 ci_lo = round(100 * (1 - sm$upper[i]), 1),
                 ci_hi = round(100 * (1 - sm$lower[i]), 1),
                 within_stratum_HR = sprintf("%.2f (%.2f-%.2f)", cs$conf.int[1, 1],
                                             cs$conf.int[1, 3], cs$conf.int[1, 4]),
                 within_stratum_p = fmt_p(cs$coefficients[1, "Pr(>|z|)"]),
                 stringsAsFactors = FALSE)
    }))
  }))
write.csv(age_strat, file.path(DIR_TAB, "table6f_age_stratified_DMS.csv"),
          row.names = FALSE, na = "")

###############################################################################
# 5. 모델 진단 — 연령 선형성, VIF, PH 가정
###############################################################################

message("== Model diagnostics ==")

# --- (a) 연령 선형성 (자연 3차 스플라인 vs 선형, LR test; robust 아닌 모델로 검정)
m_lin <- survival::coxph(stats::as.formula(paste0(
  "Surv(.time, .event) ~ ", paste(c("female", COV_M4), collapse = " + "))),
  data = d_dms, ties = "efron")
m_spl <- survival::coxph(stats::as.formula(paste0(
  "Surv(.time, .event) ~ ", paste(c("female", "ns(Age, 3)", setdiff(COV_M4, "Age")),
                                  collapse = " + "))),
  data = d_dms, ties = "efron")
lr_age <- stats::anova(m_lin, m_spl)
p_age_nonlin <- lr_age$`Pr(>|Chi|)`[2]

# --- (b) VIF (다중공선성)
#     VIF 는 설계행렬 X 에만 의존하므로 절편이 있는 선형모형으로 계산한다
#     (coxph 는 절편이 없어 car::vif 가 경고를 낸다).
vifs <- tryCatch({
  d_vif <- d_dms
  d_vif$.dummy_y <- stats::rnorm(nrow(d_vif))
  lm_vif <- stats::lm(stats::as.formula(paste0(
    ".dummy_y ~ female + ", paste(COV_M4, collapse = " + "))), data = d_vif)
  v <- car::vif(lm_vif)
  if (is.matrix(v)) {
    out <- data.frame(term = rownames(v), VIF = (v[, "GVIF^(1/(2*Df))"])^2,
                      stringsAsFactors = FALSE, row.names = NULL)
  } else {
    out <- data.frame(term = names(v), VIF = as.numeric(v),
                      stringsAsFactors = FALSE, row.names = NULL)
  }
  out
}, error = function(e) data.frame(term = NA_character_, VIF = NA_real_))

# --- (c) PH 가정: 모든 outcome, Model 4 (희소 사건은 절약형 세트)
ph_rows <- list()
for (k in names(OUT)) {
  dk <- make_surv(an, k)
  cov_use <- if (sum(dk$.event) < 150) COV_PARS else COV_M4
  mk <- survival::coxph(stats::as.formula(paste0(
    "Surv(.time, .event) ~ ", paste(c("female", cov_use), collapse = " + "))),
    data = dk, ties = "efron")
  z <- tryCatch(survival::cox.zph(mk), error = function(e) NULL)
  if (is.null(z)) next
  tb <- z$table
  ph_rows[[k]] <- data.frame(
    outcome = OUT[[k]]$lab, outcome_key = k,
    covariate_set = ifelse(sum(dk$.event) < 150, "parsimonious (age, DM, CKD, AF)", "Model 4"),
    n_events = mk$nevent,
    zph_p_female = unname(tb["female", "p"]),
    zph_p_female_fmt = fmt_p(unname(tb["female", "p"])),
    zph_p_global = unname(tb["GLOBAL", "p"]),
    zph_p_global_fmt = fmt_p(unname(tb["GLOBAL", "p"])),
    ph_violation_female = unname(tb["female", "p"]) < 0.05,
    stringsAsFactors = FALSE)
}
ph_tab <- do.call(rbind, ph_rows)

diag_out <- rbind(
  data.frame(item = "Age linearity: LR test spline(3 df) vs linear, p",
             value = fmt_p(p_age_nonlin), stringsAsFactors = FALSE),
  data.frame(item = "Max VIF among Model 4 covariates",
             value = sprintf("%.2f (%s)", max(vifs$VIF, na.rm = TRUE),
                             vifs$term[which.max(vifs$VIF)]), stringsAsFactors = FALSE),
  data.frame(item = "N covariates with VIF >= 5",
             value = as.character(sum(vifs$VIF >= 5, na.rm = TRUE)), stringsAsFactors = FALSE),
  data.frame(item = "Number of contributing centers (clustering unit)",
             value = as.character(length(unique(an$center_anon))), stringsAsFactors = FALSE),
  data.frame(item = "Events per covariate parameter, Model 4 (DMS)",
             value = sprintf("%.1f", m4_dms$nevent / length(m4_dms$coefficients)),
             stringsAsFactors = FALSE))
write.csv(diag_out, file.path(DIR_TAB, "tableS8a_model_diagnostics.csv"), row.names = FALSE, na = "")
write.csv(vifs,    file.path(DIR_TAB, "tableS8b_vif_model4.csv"), row.names = FALSE, na = "")
write.csv(ph_tab,  file.path(DIR_TAB, "table8a_ph_assumption_by_outcome.csv"), row.names = FALSE, na = "")

###############################################################################
# 6. PH 위반 대응 — (i) 구간별 time-varying coefficient, (ii) RMST
###############################################################################

message("== PH violation handling: time-varying coefficient + RMST ==")

# (i) 구간 [0,1y), [1y,3y), [3y+) 별 female HR
tv_by_interval <- function(dat, key) {
  d <- make_surv(dat, key)
  # survSplit 은 관측시간이 모두 0보다 커야 한다. .time == 0 인 행은 전부
  # 시점 0 즉시 중도절단(무추적)이라 어떤 Cox 모형에도 정보를 주지 않으므로 제외한다.
  d <- d[d$.time > 0, , drop = FALSE]
  d$.id <- seq_len(nrow(d))
  sp <- survival::survSplit(Surv(.time, .event) ~ ., data = d,
                            cut = c(365, 1095), episode = "period", id = ".id2")
  sp$period_f <- factor(sp$period, levels = 1:3,
                        labels = c("0-1 year", "1-3 years", ">3 years"))
  f <- stats::as.formula(paste0(
    "Surv(tstart, .time, .event) ~ female:period_f + strata(period_f) + ",
    paste(COV_M4, collapse = " + ")))
  m <- survival::coxph(f, data = sp, cluster = sp$center_anon, robust = TRUE, ties = "efron")
  s <- summary(m)
  idx <- grep("^female:period_f", rownames(s$coefficients))
  data.frame(outcome = OUT[[key]]$lab, outcome_key = key,
             period = sub("^female:period_f", "", rownames(s$coefficients)[idx]),
             HR = s$conf.int[idx, 1], lo = s$conf.int[idx, 3], hi = s$conf.int[idx, 4],
             HR_fmt = sprintf("%.2f (%.2f-%.2f)", s$conf.int[idx, 1],
                              s$conf.int[idx, 3], s$conf.int[idx, 4]),
             p_value = s$coefficients[idx, "Pr(>|z|)"],
             p_fmt = fmt_p(s$coefficients[idx, "Pr(>|z|)"]),
             stringsAsFactors = FALSE, row.names = NULL)
}

tv_keys <- c("DMS", "Death", "NCVdeath", "MACCE_pub", "MACE_pub", "NACE_pub")
tv_tab <- do.call(rbind, lapply(tv_keys, function(k) tv_by_interval(an, k)))
write.csv(tv_tab, file.path(DIR_TAB, "table8b_time_varying_female_HR.csv"), row.names = FALSE, na = "")

# (ii) RMST (tau = 5년), 미보정 및 공변량 보정 (survRM2 ANCOVA 방식)
rmst_one <- function(dat, key, tau = TAU5) {
  d <- make_surv(dat, key)
  # 추적 0일(시점 0 즉시 중도절단) 행은 RMST 추정에 정보를 주지 않으며
  # survRM2 의 IPCW 회귀가 time > 0 을 요구하므로 제외한다.
  d <- d[d$.time > 0, , drop = FALSE]
  # tau 는 각 군 최대 관측시간의 최소값을 넘을 수 없음
  tmax <- min(tapply(d$.time, d$female, max))
  tau_use <- min(tau, tmax)
  X <- stats::model.matrix(stats::as.formula(paste0("~ ", paste(COV_M4, collapse = " + "))),
                           data = d)[, -1, drop = FALSE]
  un <- survRM2::rmst2(time = d$.time, status = d$.event, arm = d$female, tau = tau_use)
  ad <- tryCatch(survRM2::rmst2(time = d$.time, status = d$.event, arm = d$female,
                                tau = tau_use, covariates = X),
                 error = function(e) NULL)
  out <- data.frame(
    outcome = OUT[[key]]$lab, outcome_key = key, tau_days = tau_use,
    n_analyzed = nrow(d),
    rmst_male_days   = un$RMST.arm0$rmst["Est."],
    rmst_female_days = un$RMST.arm1$rmst["Est."],
    rmst_diff_days = un$unadjusted.result[1, 1],
    rmst_diff_lo = un$unadjusted.result[1, 2], rmst_diff_hi = un$unadjusted.result[1, 3],
    rmst_diff_p = un$unadjusted.result[1, 4],
    stringsAsFactors = FALSE, row.names = NULL)
  out$rmst_diff_fmt <- sprintf("%.1f (%.1f to %.1f)", out$rmst_diff_days,
                               out$rmst_diff_lo, out$rmst_diff_hi)
  if (!is.null(ad)) {
    r <- ad$RMST.difference.adjusted
    out$adj_rmst_diff_days <- r["arm", "coef"]
    out$adj_rmst_diff_lo   <- r["arm", "lower .95"]
    out$adj_rmst_diff_hi   <- r["arm", "upper .95"]
    out$adj_rmst_diff_p    <- r["arm", "p"]
    out$adj_rmst_diff_fmt  <- sprintf("%.1f (%.1f to %.1f)", r["arm", "coef"],
                                      r["arm", "lower .95"], r["arm", "upper .95"])
  }
  out
}
rmst_tab <- do.call(rbind, lapply(c("DMS", "Death", "NCVdeath", "CVdeath",
                                    "MACCE_pub", "MACE_pub", "NACE_pub", "AnyRevasc"),
                                  function(k) tryCatch(rmst_one(an, k), error = function(e) NULL)))
write.csv(rmst_tab, file.path(DIR_TAB, "table8c_rmst_5yr.csv"), row.names = FALSE, na = "")

###############################################################################
# 7. Table 7 — 모든 outcome 에 대한 보정 HR (Model 4), secondary outcome 전부 포함
#    ⛔ 사전 계획을 시사하는 표현 금지 (CLAUDE.md §1-2)
###############################################################################

message("== Table 7: adjusted HR for all reported outcomes ==")

adj_all <- do.call(rbind, lapply(names(OUT), function(k) {
  dk <- make_surv(an, k)
  nev <- sum(dk$.event)
  r <- rbind(
    fit_cox(dk, character(0), "Unadjusted",       OUT[[k]]$lab),
    fit_cox(dk, COV_M2,       "Age-adjusted",     OUT[[k]]$lab),
    fit_cox(dk, COV_PARS,     "Parsimoniously adjusted (age, DM, CKD, AF)", OUT[[k]]$lab),
    fit_cox(dk, COV_M4,       "Fully adjusted",   OUT[[k]]$lab))
  r$outcome_key <- k
  # events per parameter — 희소 사건에서 Model 4 의 신뢰성 경고
  r$events_per_parameter <- round(nev / c(1, 2, length(COV_PARS) + 1, 26), 1)
  r$sparse_warning <- ifelse(r$events_per_parameter < 10,
                             "events per parameter < 10; estimate may be unstable", "")
  r
}))
write.csv(adj_all, file.path(DIR_TAB, "table7_adjusted_HR_all_outcomes.csv"),
          row.names = FALSE, na = "")

# 원고용 wide 포맷 (outcome 정의 민감도 = 원고 방어 핵심)
adj_wide <- adj_all %>%
  dplyr::select(outcome_key, outcome, model, HR_fmt, p_fmt) %>%
  tidyr::pivot_wider(names_from = model, values_from = c(HR_fmt, p_fmt))
write.csv(adj_wide, file.path(DIR_TAB, "table7b_outcome_definition_sensitivity.csv"),
          row.names = FALSE, na = "")

###############################################################################
# 8. 경쟁위험 — Fine-Gray (CV death / non-CV death 상호 경쟁)
###############################################################################

message("== Competing risk: Fine-Gray subdistribution hazard ==")

finegray_sex <- function(dat, cause = 1, covs = COV_M4, meth_lab = "fully adjusted (Model 4)") {
  d <- dat
  t <- d$death_fudur
  ev <- ifelse(d$CVdeath == 1, 1, ifelse(d$NCVdeath == 1, 2, 0))
  t[ev > 0 & t == 0] <- 0.5
  ok <- !is.na(t) & !is.na(ev)
  d <- d[ok, , drop = FALSE]; t <- t[ok]; ev <- ev[ok]
  X <- stats::model.matrix(stats::as.formula(paste0("~ female + ", paste(covs, collapse = " + "))),
                           data = d)[, -1, drop = FALSE]
  fg <- cmprsk::crr(ftime = t, fstatus = ev, cov1 = X, failcode = cause, cencode = 0)
  s <- summary(fg)
  i <- which(rownames(s$coef) == "female")
  data.frame(
    event_type = ifelse(cause == 1, "Cardiovascular death", "Non-cardiovascular death"),
    method = paste0("Fine-Gray subdistribution hazard (", meth_lab, ")"),
    n = length(t), n_events = sum(ev == cause), n_competing = sum(ev == (3 - cause)),
    sHR = unname(s$conf.int[i, 1]),
    lo = unname(s$conf.int[i, 3]), hi = unname(s$conf.int[i, 4]),
    sHR_fmt = sprintf("%.2f (%.2f-%.2f)", s$conf.int[i, 1], s$conf.int[i, 3], s$conf.int[i, 4]),
    p_value = unname(s$coef[i, 5]), p_fmt = fmt_p(unname(s$coef[i, 5])),
    stringsAsFactors = FALSE)
}
fg_tab <- rbind(
  tryCatch(finegray_sex(an, 1, COV_PARS, "parsimoniously adjusted: age, DM, CKD, AF"), error = function(e) NULL),
  tryCatch(finegray_sex(an, 2, COV_PARS, "parsimoniously adjusted: age, DM, CKD, AF"), error = function(e) NULL),
  tryCatch(finegray_sex(an, 1), error = function(e) NULL),
  tryCatch(finegray_sex(an, 2), error = function(e) NULL))

# 미보정 sHR (참고)
finegray_crude <- function(dat, cause = 1) {
  t <- dat$death_fudur
  ev <- ifelse(dat$CVdeath == 1, 1, ifelse(dat$NCVdeath == 1, 2, 0))
  t[ev > 0 & t == 0] <- 0.5
  ok <- !is.na(t) & !is.na(ev)
  fg <- cmprsk::crr(ftime = t[ok], fstatus = ev[ok],
                    cov1 = matrix(dat$female[ok], ncol = 1, dimnames = list(NULL, "female")),
                    failcode = cause, cencode = 0)
  s <- summary(fg)
  data.frame(event_type = ifelse(cause == 1, "Cardiovascular death", "Non-cardiovascular death"),
             method = "Fine-Gray subdistribution hazard (unadjusted)",
             n = sum(ok), n_events = sum(ev[ok] == cause), n_competing = sum(ev[ok] == (3 - cause)),
             sHR = unname(s$conf.int[1, 1]), lo = unname(s$conf.int[1, 3]), hi = unname(s$conf.int[1, 4]),
             sHR_fmt = sprintf("%.2f (%.2f-%.2f)", s$conf.int[1, 1], s$conf.int[1, 3], s$conf.int[1, 4]),
             p_value = unname(s$coef[1, 5]), p_fmt = fmt_p(unname(s$coef[1, 5])),
             stringsAsFactors = FALSE)
}
fg_tab <- rbind(finegray_crude(an, 1), finegray_crude(an, 2), fg_tab)
write.csv(fg_tab, file.path(DIR_TAB, "table9a_finegray_competing_risk.csv"), row.names = FALSE, na = "")

# Aalen-Johansen CIF (보정 없음, 분석표본 기준) + Gray 검정
aj <- local({
  t <- an$death_fudur
  ev <- ifelse(an$CVdeath == 1, 1, ifelse(an$NCVdeath == 1, 2, 0))
  t[ev > 0 & t == 0] <- 0.5
  ci <- cmprsk::cuminc(ftime = t, fstatus = ev, group = as.character(an$sex_f), cencode = 0)
  est <- cmprsk::timepoints(ci, times = c(365, 1095, 1825))
  gp <- ci$Tests[, "pv"]
  out <- list()
  for (rn in rownames(est$est)) for (j in 1:3) {
    out[[length(out) + 1]] <- data.frame(
      sex = ifelse(grepl("^Male", rn), "Male", "Female"),
      event_type = ifelse(grepl(" 1$", rn), "CV death", "Non-CV death"),
      timepoint = c("1-year", "3-year", "5-year")[j],
      cif_pct = round(100 * est$est[rn, j], 2),
      se_pct = round(100 * sqrt(est$var[rn, j]), 2),
      gray_p = ifelse(grepl(" 1$", rn), gp["1"], gp["2"]),
      stringsAsFactors = FALSE)
  }
  r <- do.call(rbind, out); r$gray_p_fmt <- fmt_p(r$gray_p); r
})
write.csv(aj, file.path(DIR_TAB, "table9b_aalen_johansen_CIF.csv"), row.names = FALSE, na = "")

###############################################################################
# 9. 인과추론 보강 — IPTW(표준화 가중) · PSM · E-value
#
# ⚠️ 성별은 개입(intervention)이 아니다. 따라서 아래 가중은 엄밀한 의미의
#    propensity score weighting(치료 확률의 역수)이 아니라, 공변량 분포를
#    동일하게 맞춘 **표준화(standardization) 가중**으로 해석해야 한다.
#    "여성에게 남성의 공변량 분포를 부여했을 때"의 반사실은 정의되지 않으며,
#    결과는 "관측된 임상 프로파일 차이를 제거했을 때 남는 성별 관련 차이"로만
#    기술한다. 인과 용어 금지.
###############################################################################

message("== IPTW (standardization weights) and PSM ==")

ps_form <- stats::as.formula(paste("female ~", paste(COV_PS, collapse = " + ")))

# --- (a) IPTW ---------------------------------------------------------------
W <- WeightIt::weightit(ps_form, data = an, method = "glm", estimand = "ATE",
                        stabilize = TRUE)
an$w_iptw <- W$weights
an$ps     <- W$ps

bal <- cobalt::bal.tab(W, un = TRUE, stats = c("mean.diffs"), binary = "std",
                       s.d.denom = "pooled")
bal_df <- data.frame(variable = rownames(bal$Balance), bal$Balance,
                     stringsAsFactors = FALSE, row.names = NULL)

# 추가 진단: 체격 변수의 잔여 불균형 (positivity 위반 근거)
extra_bal <- cobalt::bal.tab(as.formula("female ~ Height + Weight + bsa + BMI_re"),
                             data = an, weights = an$w_iptw, un = TRUE,
                             s.d.denom = "pooled")
extra_df <- data.frame(variable = rownames(extra_bal$Balance), extra_bal$Balance,
                       stringsAsFactors = FALSE, row.names = NULL)
extra_df$note <- "Body-size variables: NOT included in the weighting model (positivity violation); shown to document non-overlap"
bal_df$note <- "Included in the weighting model"
write.csv(rbind(bal_df[, intersect(names(bal_df), c("variable", "Type", "Diff.Un", "Diff.Adj", "note"))],
                extra_df[, intersect(names(extra_df), c("variable", "Type", "Diff.Un", "Diff.Adj", "note"))]),
          file.path(DIR_TAB, "tableS9a_balance_iptw.csv"), row.names = FALSE, na = "")

ps_diag <- data.frame(
  item = c("PS range, male", "PS range, female",
           "Stabilized weight: min / median / max",
           "N with stabilized weight > 10",
           "Common support (PS overlap region)"),
  value = c(sprintf("%.4f - %.4f", min(an$ps[an$female == 0]), max(an$ps[an$female == 0])),
            sprintf("%.4f - %.4f", min(an$ps[an$female == 1]), max(an$ps[an$female == 1])),
            sprintf("%.3f / %.3f / %.3f", min(an$w_iptw), median(an$w_iptw), max(an$w_iptw)),
            as.character(sum(an$w_iptw > 10)),
            sprintf("%.4f - %.4f", max(min(an$ps[an$female == 0]), min(an$ps[an$female == 1])),
                    min(max(an$ps[an$female == 0]), max(an$ps[an$female == 1])))),
  stringsAsFactors = FALSE)
write.csv(ps_diag, file.path(DIR_TAB, "tableS9b_ps_diagnostics.csv"), row.names = FALSE, na = "")

# --- (a-2) 엔트로피 균형화 (평균 정확 균형; 로지스틱 PS 의 잔여 불균형 보완) --
Web <- tryCatch(WeightIt::weightit(ps_form, data = an, method = "ebal", estimand = "ATE"),
                error = function(e) NULL)
if (!is.null(Web)) {
  an$w_ebal <- Web$weights
  bal_e <- cobalt::bal.tab(Web, un = TRUE, s.d.denom = "pooled")
  bal_e_df <- data.frame(variable = rownames(bal_e$Balance), bal_e$Balance,
                         stringsAsFactors = FALSE, row.names = NULL)
  write.csv(bal_e_df, file.path(DIR_TAB, "tableS9e_balance_entropy.csv"),
            row.names = FALSE, na = "")
} else an$w_ebal <- NA_real_

iptw_fit <- function(key, trim = FALSE, wvar = "w_iptw") {
  d <- make_surv(an, key)
  w <- d[[wvar]]
  lab <- if (wvar == "w_ebal") "Entropy balancing weights (ATE)"
         else "IPTW (standardization weights, ATE)"
  if (trim) {
    q <- stats::quantile(w, c(0.01, 0.99))
    w <- pmin(pmax(w, q[1]), q[2]); lab <- "IPTW, weights trimmed at 1st/99th pct"
  }
  m <- survival::coxph(Surv(.time, .event) ~ female, data = d, weights = w,
                       cluster = d$center_anon, robust = TRUE)
  s <- summary(m)
  data.frame(outcome = OUT[[key]]$lab, outcome_key = key, method = lab,
             n = m$n, n_events = m$nevent,
             HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4],
             HR_fmt = sprintf("%.2f (%.2f-%.2f)", s$conf.int[1, 1], s$conf.int[1, 3], s$conf.int[1, 4]),
             p_value = s$coefficients[1, "Pr(>|z|)"],
             p_fmt = fmt_p(s$coefficients[1, "Pr(>|z|)"]), stringsAsFactors = FALSE)
}

# --- (b) PSM ----------------------------------------------------------------
M <- MatchIt::matchit(ps_form, data = an, method = "nearest", distance = "glm",
                      caliper = 0.2, std.caliper = TRUE, ratio = 1, replace = FALSE)
md <- MatchIt::match.data(M)
bal_m <- cobalt::bal.tab(M, un = TRUE, s.d.denom = "pooled")
bal_m_df <- data.frame(variable = rownames(bal_m$Balance), bal_m$Balance,
                       stringsAsFactors = FALSE, row.names = NULL)
bal_m_df$note <- "Included in the matching model"
extra_m <- cobalt::bal.tab(as.formula("female ~ Height + Weight + bsa + BMI_re"),
                           data = md, un = TRUE, s.d.denom = "pooled")
extra_m_df <- data.frame(variable = rownames(extra_m$Balance), extra_m$Balance,
                         stringsAsFactors = FALSE, row.names = NULL)
extra_m_df$note <- "Body-size variables in matched set: NOT matched on (non-overlap)"
write.csv(dplyr::bind_rows(bal_m_df, extra_m_df),
          file.path(DIR_TAB, "tableS9c_balance_psm.csv"), row.names = FALSE, na = "")

psm_fit <- function(key) {
  d <- make_surv(md, key)
  m <- survival::coxph(Surv(.time, .event) ~ female, data = d, cluster = d$subclass, robust = TRUE)
  s <- summary(m)
  data.frame(outcome = OUT[[key]]$lab, outcome_key = key,
             method = "1:1 propensity score matching (caliper 0.2 SD)",
             n = m$n, n_events = m$nevent,
             HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4],
             HR_fmt = sprintf("%.2f (%.2f-%.2f)", s$conf.int[1, 1], s$conf.int[1, 3], s$conf.int[1, 4]),
             p_value = s$coefficients[1, "Pr(>|z|)"],
             p_fmt = fmt_p(s$coefficients[1, "Pr(>|z|)"]), stringsAsFactors = FALSE)
}

KEYS_CAUSAL <- c("DMS", "Death", "NCVdeath", "CVdeath", "MACCE_pub", "MACE_pub",
                 "NACE_pub", "AnyRevasc")
causal_tab <- dplyr::bind_rows(
  lapply(KEYS_CAUSAL, function(k) {
    dk <- make_surv(an, k)
    r0 <- fit_cox(dk, COV_M4, "Multivariable Cox (Model 4)", OUT[[k]]$lab)
    r0$outcome_key <- k; r0$method <- "Multivariable Cox (Model 4)"
    r0 <- r0[, c("outcome", "outcome_key", "method", "n", "n_events",
                 "HR", "lo", "hi", "HR_fmt", "p_value", "p_fmt")]
    rows <- list(r0, iptw_fit(k), iptw_fit(k, trim = TRUE), psm_fit(k))
    if (!all(is.na(an$w_ebal))) rows[[length(rows) + 1]] <- iptw_fit(k, wvar = "w_ebal")
    dplyr::bind_rows(rows)
  }))
write.csv(causal_tab, file.path(DIR_TAB, "table10_causal_inference_comparison.csv"),
          row.names = FALSE, na = "")

psm_diag <- data.frame(
  item = c("Matched pairs", "Females matched", "Females unmatched",
           "Max |SMD| after matching (matched covariates)",
           "Max |SMD| after matching (body-size, unmatched-on)"),
  value = c(as.character(sum(!is.na(md$subclass)) / 2),
            as.character(sum(md$female == 1)),
            as.character(sum(an$female == 1) - sum(md$female == 1)),
            sprintf("%.3f", max(abs(bal_m_df$Diff.Adj), na.rm = TRUE)),
            sprintf("%.3f", max(abs(extra_m_df$Diff.Un), na.rm = TRUE))),
  stringsAsFactors = FALSE)
write.csv(psm_diag, file.path(DIR_TAB, "tableS9d_psm_diagnostics.csv"), row.names = FALSE, na = "")

# --- (c) E-value -------------------------------------------------------------
message("== E-values ==")
# ⚠️ F7 (Phase 2.5 QA): 미보정 모형에는 보정된 공변량이 없으므로 해석 문장에
#    "above and beyond the measured covariates" 를 붙이면 안 된다. adjusted 인자로
#    문장을 분기한다.
evalue_row <- function(hr, lo, hi, lab, method, adjusted = TRUE) {
  e <- EValue::evalues.HR(est = hr, lo = lo, hi = hi, rare = FALSE)
  ci_bound <- if (hr > 1) e["E-values", "lower"] else e["E-values", "upper"]
  qualifier <- if (adjusted) ", above and beyond the measured covariates," else ""
  data.frame(outcome = lab, method = method,
             HR = round(hr, 3), lo = round(lo, 3), hi = round(hi, 3),
             evalue_point = round(e["E-values", "point"], 2),
             evalue_ci = round(ci_bound, 2),
             interpretation = sprintf(
               paste0("An unmeasured confounder would need to be associated with both female sex and the outcome ",
                      "by a risk ratio of %.2f-fold each%s ",
                      "to fully explain the observed HR of %.2f; %.2f-fold to move the confidence interval to the null."),
               e["E-values", "point"], qualifier, hr, ci_bound),
             stringsAsFactors = FALSE)
}
ev_src <- causal_tab %>% filter(method == "Multivariable Cox (Model 4)")
ev_tab <- do.call(rbind, lapply(seq_len(nrow(ev_src)), function(i)
  evalue_row(ev_src$HR[i], ev_src$lo[i], ev_src$hi[i], ev_src$outcome[i],
             "Fully adjusted Cox (Model 4)")))
ev_tab <- rbind(ev_tab, do.call(rbind, lapply(which(causal_tab$method == "IPTW (standardization weights, ATE)"),
  function(i) evalue_row(causal_tab$HR[i], causal_tab$lo[i], causal_tab$hi[i],
                         causal_tab$outcome[i], "IPTW"))))

# 미보정 추정치의 E-value — "측정된 연령만으로 초과위험이 설명되었다"는 사실을
# 정량화하기 위한 대조군. 미보정 HR 을 설명하는 데 필요한 교란 강도와,
# 실제 관측된 연령-사망 연관 강도를 비교해 제시한다.
ev_crude <- do.call(rbind, lapply(c("DMS", "Death", "NCVdeath"), function(k) {
  r <- adj_all[adj_all$outcome_key == k & adj_all$model == "Unadjusted", ]
  evalue_row(r$HR, r$lo, r$hi, r$outcome, "Unadjusted (for comparison)", adjusted = FALSE)
}))
ev_tab <- rbind(ev_crude, ev_tab)

# 귀무 결과에 대한 해석 보조: 상한 CI 를 기준으로 "숨어 있을 수 있는 최대 유해효과"
ev_tab$null_result <- ev_tab$lo < 1 & ev_tab$hi > 1
ev_tab$note_for_null <- ifelse(
  ev_tab$null_result,
  paste0("Null result: the confidence interval already includes 1. ",
         "The E-value for the point estimate quantifies how much unmeasured confounding ",
         "would be needed to move the observed estimate to the null, which is small by ",
         "construction. The upper confidence limit (", sprintf("%.2f", ev_tab$hi),
         ") bounds the largest excess risk compatible with these data."),
  "Non-null result: E-value quantifies the unmeasured confounding needed to explain it away.")
write.csv(ev_tab, file.path(DIR_TAB, "table11_evalues.csv"), row.names = FALSE, na = "")

# 실제 관측된 연령의 연관 강도 (E-value 와의 대조용)
age_strength <- local({
  m <- cox_model(d_dms, COV_M4)
  s <- summary(m)
  hr10 <- exp(10 * coef(m)["Age"])
  data.frame(item = c("Adjusted HR per 10-year increase in age (DMS)",
                      "Median age, male vs female (years)",
                      "Implied crude HR contrast for the observed 9-year age gap"),
             value = c(sprintf("%.2f (%.2f-%.2f)", hr10,
                               exp(10 * (coef(m)["Age"] - 1.96 * s$coefficients["Age", "robust se"])),
                               exp(10 * (coef(m)["Age"] + 1.96 * s$coefficients["Age", "robust se"]))),
                       sprintf("%.0f vs %.0f", stats::median(an$Age[an$female == 0], na.rm = TRUE),
                               stats::median(an$Age[an$female == 1], na.rm = TRUE)),
                       sprintf("%.2f", exp(coef(m)["Age"] *
                         (stats::median(an$Age[an$female == 1], na.rm = TRUE) -
                          stats::median(an$Age[an$female == 0], na.rm = TRUE))))),
             stringsAsFactors = FALSE)
})
write.csv(age_strength, file.path(DIR_TAB, "table11b_age_association_strength.csv"),
          row.names = FALSE, na = "")

###############################################################################
# 10. Subgroup 분석 (primary outcome = DMS)
#     ⛔ 사전 계획을 시사하는 표현 금지 (CLAUDE.md §1-2). 본 분석은 사전등록되지
#        않았으므로 "본 분석에서 정의한 6개 subgroup" 으로만 기술한다.
###############################################################################

message("== Subgroup analyses (6 subgroups defined for this analysis) ==")

SUBGROUPS <- list(
  list(var = "age75",       lab = "Age",             levs = c("<75 years", ">=75 years"),
       drop = character(0)),
  list(var = "DM",          lab = "Diabetes",        levs = c("No", "Yes"), drop = "DM"),
  list(var = "CKD_2",       lab = "CKD (eGFR<60)",   levs = c("No", "Yes"), drop = "CKD_2"),
  list(var = "cto_success", lab = "Technical success", levs = c("Failure", "Success"),
       drop = character(0)),
  list(var = "JCTO_bi",     lab = "J-CTO score",     levs = c("<2", ">=2"), drop = "JCTO_re"),
  list(var = "index_half",  lab = "Index period",    levs = c("2014-2018", "2019-2023"),
       drop = "index_period_f")
)

subgroup_analysis <- function(dat, key) {
  d0 <- make_surv(dat, key)
  d0$index_half <- as.numeric(d0$index_half_f == "2019-2023")
  res <- list()
  for (sg in SUBGROUPS) {
    v <- sg$var
    covs <- setdiff(COV_M4, sg$drop)
    dd <- d0[!is.na(d0[[v]]), , drop = FALSE]
    # 층별 HR
    for (lv in 0:1) {
      ds <- dd[dd[[v]] == lv, , drop = FALSE]
      r <- tryCatch(fit_cox(ds, covs, sg$levs[lv + 1], OUT[[key]]$lab),
                    error = function(e) NULL)
      if (is.null(r)) next
      r$subgroup <- sg$lab; r$subgroup_var <- v; r$level <- sg$levs[lv + 1]
      res[[length(res) + 1]] <- r
    }
    # 상호작용 (Wald, cluster-robust)
    f <- stats::as.formula(paste0("Surv(.time, .event) ~ female * ", v, " + ",
                                  paste(covs, collapse = " + ")))
    mi <- tryCatch(survival::coxph(f, data = dd, cluster = dd$center_anon,
                                   robust = TRUE, ties = "efron"),
                   error = function(e) NULL)
    p_int <- NA_real_
    if (!is.null(mi)) {
      s <- summary(mi)
      ii <- grep(paste0("^female:", v, "$|^", v, ":female$"), rownames(s$coefficients))
      if (length(ii) == 1) p_int <- s$coefficients[ii, "Pr(>|z|)"]
    }
    for (j in seq_along(res)) if (res[[j]]$subgroup_var[1] == v) res[[j]]$p_interaction <- p_int
  }
  out <- do.call(rbind, res)
  out$p_interaction_fmt <- fmt_p(out$p_interaction)
  out[, c("outcome", "subgroup", "subgroup_var", "level", "n", "n_events",
          "events_female", "events_male", "HR", "lo", "hi", "HR_fmt",
          "p_value", "p_fmt", "p_interaction", "p_interaction_fmt")]
}

sg_dms   <- subgroup_analysis(an, "DMS")
sg_death <- subgroup_analysis(an, "Death")
sg_all <- rbind(sg_dms, sg_death)
# 다중비교: 12개 상호작용 검정 (6 subgroup x 2 outcome) -> Holm 보정 병기
int_only <- unique(sg_all[, c("outcome", "subgroup_var", "p_interaction")])
int_only$p_interaction_holm <- stats::p.adjust(int_only$p_interaction, method = "holm")
sg_all <- dplyr::left_join(sg_all, int_only, by = c("outcome", "subgroup_var", "p_interaction"))
sg_all$p_interaction_holm_fmt <- fmt_p(sg_all$p_interaction_holm)
sg_all$multiplicity_note <- sprintf(
  "%d interaction tests performed (6 subgroups defined for this analysis x 2 outcomes); Holm-adjusted p reported",
  nrow(int_only))
write.csv(sg_all, file.path(DIR_TAB, "table12_subgroup_analysis.csv"), row.names = FALSE, na = "")

###############################################################################
# 11. 민감도 분석 배터리
###############################################################################

message("== Sensitivity battery ==")

sens <- list()
addS <- function(x) sens[[length(sens) + 1]] <<- x

# (S1) 주 결과 재게시 (참조행)
r <- fit_cox(make_surv(an, "DMS"), COV_M4, "Primary analysis (Model 4, complete case)", OUT$DMS$lab)
r$sensitivity <- "Reference"; addS(r)

# (S2) day-0 event 제외 규칙
r <- fit_cox(make_surv(an, "DMS", day0 = "exclude"), COV_M4,
             "Day-0 events excluded (instead of 0.5-day imputation)", OUT$DMS$lab)
r$sensitivity <- "Day-0 event handling"; addS(r)

# (S2b) 연령 스플라인 (선형성 검정 기각 대응)
r <- fit_cox(make_surv(an, "DMS"), COV_M4S, "Model 4 with natural cubic spline for age", OUT$DMS$lab)
r$sensitivity <- "Flexible age adjustment"; addS(r)
r <- fit_cox(make_surv(an, "DMS"), COV_M2S, "Age only, natural cubic spline", OUT$DMS$lab)
r$sensitivity <- "Flexible age adjustment"; addS(r)

# (S3) LVEF 추가 (표본 감소)
d4l <- make_surv(co[stats::complete.cases(co[, COV_M4L, drop = FALSE]), ], "DMS")
r <- fit_cox(d4l, COV_M4L, "Model 4 + LVEF (complete case)", OUT$DMS$lab)
r$sensitivity <- "Additional covariate: LVEF"; addS(r)

# (S4) 등록기간 2014-2022 제한
d_win <- make_surv(an[an$index_yr <= 2022, ], "DMS")
r <- fit_cox(d_win, COV_M4, "Restricted to index year 2014-2022", OUT$DMS$lab)
r$sensitivity <- "Registration window"; addS(r)

# (S5) 5년 행정적 중도절단
d_adm <- make_surv(an, "DMS")
d_adm$.event <- ifelse(d_adm$.time > TAU5, 0, d_adm$.event)
d_adm$.time  <- pmin(d_adm$.time, TAU5)
r <- fit_cox(d_adm, COV_M4, "Administrative censoring at 5 years", OUT$DMS$lab)
r$sensitivity <- "Follow-up truncation"; addS(r)

# (S6) center shared frailty (clustered robust SE 대안)
d_fr <- make_surv(an, "DMS")
m_fr <- survival::coxph(stats::as.formula(paste0(
  "Surv(.time, .event) ~ female + ", paste(COV_M4, collapse = " + "),
  " + frailty(center_anon)")), data = d_fr, ties = "efron")
sfr <- summary(m_fr)
r <- data.frame(outcome = OUT$DMS$lab, model = "Shared gamma frailty by center",
                n = m_fr$n, n_events = m_fr$nevent,
                n_female = sum(d_fr$female == 1),
                events_female = sum(d_fr$.event[d_fr$female == 1] == 1),
                n_male = sum(d_fr$female == 0),
                events_male = sum(d_fr$.event[d_fr$female == 0] == 1),
                HR = exp(coef(m_fr)["female"]),
                lo = exp(coef(m_fr)["female"] - 1.96 * sqrt(m_fr$var[1, 1])),
                hi = exp(coef(m_fr)["female"] + 1.96 * sqrt(m_fr$var[1, 1])),
                HR_fmt = NA_character_, se_robust = sqrt(m_fr$var[1, 1]),
                p_value = sfr$coefficients["female", "p"],
                p_fmt = fmt_p(sfr$coefficients["female", "p"]),
                note = sprintf("frailty variance = %s",
                               tryCatch(sprintf("%.4f", m_fr$history[[1]]$theta),
                                        error = function(e) "not extractable")),
                stringsAsFactors = FALSE, row.names = NULL)
r$HR_fmt <- sprintf("%.2f (%.2f-%.2f)", r$HR, r$lo, r$hi)
r$sensitivity <- "Clustering approach"; addS(r)

# (S7) 표준 (non-robust) SE — clustered robust SE 의 영향 정량화
d_nr <- make_surv(an, "DMS")
m_nr <- survival::coxph(stats::as.formula(paste0(
  "Surv(.time, .event) ~ female + ", paste(COV_M4, collapse = " + "))),
  data = d_nr, ties = "efron")
snr <- summary(m_nr)
r <- data.frame(outcome = OUT$DMS$lab, model = "Model-based (non-robust) SE",
                n = m_nr$n, n_events = m_nr$nevent,
                n_female = sum(d_nr$female == 1),
                events_female = sum(d_nr$.event[d_nr$female == 1] == 1),
                n_male = sum(d_nr$female == 0),
                events_male = sum(d_nr$.event[d_nr$female == 0] == 1),
                HR = snr$conf.int["female", 1], lo = snr$conf.int["female", 3],
                hi = snr$conf.int["female", 4],
                HR_fmt = sprintf("%.2f (%.2f-%.2f)", snr$conf.int["female", 1],
                                 snr$conf.int["female", 3], snr$conf.int["female", 4]),
                se_robust = snr$coefficients["female", "se(coef)"],
                p_value = snr$coefficients["female", "Pr(>|z|)"],
                p_fmt = fmt_p(snr$coefficients["female", "Pr(>|z|)"]),
                note = "no clustering adjustment", stringsAsFactors = FALSE, row.names = NULL)
r$sensitivity <- "Clustering approach"; addS(r)

# (S8) 전체 코호트(결측 포함) 미보정 — complete-case 선택 편향 점검
r <- fit_cox(make_surv(co, "DMS"), character(0), "Unadjusted, all 6,777 (incl. incomplete cases)", OUT$DMS$lab)
r$sensitivity <- "Complete-case selection"; addS(r)
r <- fit_cox(make_surv(an, "DMS"), character(0), "Unadjusted, complete-case analysis set", OUT$DMS$lab)
r$sensitivity <- "Complete-case selection"; addS(r)

sens_tab <- dplyr::bind_rows(sens)
write.csv(sens_tab, file.path(DIR_TAB, "tableS10_sensitivity_battery_DMS.csv"),
          row.names = FALSE, na = "")

###############################################################################
# 12. 다중대치 (mice) 민감도 — 산발적 결측만 대치
#     구조적 결측(역행성 미시도, 스텐트 미삽입 등)은 대치 대상 아님
###############################################################################

message("== Multiple imputation (mice, m=20, seed 42) ==")

# center_anon 은 Center_volume_f 를 완전히 결정하므로 대치모형에서 제외한다
# (포함 시 완전 다중공선성 -> ridge penalty 경고 발생).
mi_vars <- unique(c(COV_M4, "preEF", "female"))
mi_df <- co[, mi_vars, drop = FALSE]
tmp <- apply_day0(co$DMS_fudur, co$DMS, "half_day")
mi_df$.time  <- tmp$time
mi_df$.event <- tmp$event
# van Buuren 권고: 생존 outcome 은 Nelson-Aalen 누적위험 + event indicator 로 대치모형에 포함
na_est <- survival::basehaz(survival::coxph(Surv(.time, .event) ~ 1, data = mi_df), centered = FALSE)
mi_df$.nelsonaalen <- stats::approx(na_est$time, na_est$hazard, xout = mi_df$.time,
                                    method = "constant", rule = 2)$y

pred <- mice::make.predictorMatrix(mi_df)
pred[, ".time"] <- 0                       # raw time 은 예측자에서 제외
pred[".time", ] <- 0; pred[".event", ] <- 0; pred[".nelsonaalen", ] <- 0
meth <- mice::make.method(mi_df)
meth[c(".time", ".event", ".nelsonaalen")] <- ""

imp <- mice::mice(mi_df, m = 20, maxit = 10, printFlag = FALSE, seed = 42,
                  predictorMatrix = pred, method = meth)

mi_fit <- with(imp, survival::coxph(
  Surv(.time, .event) ~ female + Age + BMI_re + Current_sm + HTN + DM + CKD_2 +
    prePCI + preCABG + preMI + preCVA + prePVD + AF + stable_angina +
    JCTO_re + CTO_location_pm + ISR_CTO + Rentrop3 + non_CTO.PCI +
    Vessel.Disease_f + index_period_f + Center_volume_f + preEF))
pool_obj <- mice::pool(mi_fit)
pl <- summary(pool_obj, conf.int = TRUE)
prow <- pl[pl$term == "female", ]
fmi_female <- tryCatch(pool_obj$pooled$fmi[pool_obj$pooled$term == "female"],
                       error = function(e) NA_real_)
mi_tab <- data.frame(
  outcome = OUT$DMS$lab,
  model = "Model 4 + LVEF, multiple imputation (m=20, mice, seed 42)",
  n_imputed_datasets = 20, n = nrow(mi_df),
  HR = exp(prow$estimate), lo = exp(prow$`2.5 %`), hi = exp(prow$`97.5 %`),
  HR_fmt = sprintf("%.2f (%.2f-%.2f)", exp(prow$estimate), exp(prow$`2.5 %`), exp(prow$`97.5 %`)),
  p_value = prow$p.value, p_fmt = fmt_p(prow$p.value),
  fmi = round(as.numeric(fmi_female), 3),
  note = "SE not clustered by center (mice pooling); compare with complete-case robust model",
  stringsAsFactors = FALSE)
write.csv(mi_tab, file.path(DIR_TAB, "tableS11_multiple_imputation_DMS.csv"),
          row.names = FALSE, na = "")

# 결측률 요약 (대치 대상 변수)
miss_sum <- data.frame(variable = mi_vars,
                       pct_missing = round(100 * sapply(mi_vars, function(v) mean(is.na(co[[v]]))), 2),
                       stringsAsFactors = FALSE, row.names = NULL)
miss_sum$missing_type <- ifelse(miss_sum$pct_missing == 0, "none", "sporadic (imputed)")
write.csv(miss_sum, file.path(DIR_TAB, "tableS11b_imputed_variable_missingness.csv"),
          row.names = FALSE, na = "")

###############################################################################
# 13. 원내 결과 · 시술 성공의 보정 분석 (로지스틱, cluster-robust)
###############################################################################

message("== Adjusted logistic models for procedural / in-hospital outcomes ==")

# logit_sex() / crude_or() 는 00_cohort_prep.R 에 정의 (04 / 05 공용)
# 시술 성공은 시술 변수를 공변량으로 쓰면 안 되므로 Model 4 공변량 사용
inhos_tab <- dplyr::bind_rows(
  logit_sex(an, "cto_success", "Technical success of CTO PCI"),
  logit_sex(an, "inhos_any",   "Any in-hospital adverse event"),
  logit_sex(an, "inhos_death", "In-hospital death"),
  logit_sex(an, "PeriMI",      "Periprocedural myocardial infarction"),
  logit_sex(an, "inhos_CVA",   "In-hospital stroke"),
  logit_sex(an, "CIN_2",       "Contrast-induced nephropathy"),
  logit_sex(an, "Final_TIMI3", "Final TIMI 3 flow"),
  logit_sex(an, "CR_yes",      "Complete revascularization"))
# 미보정 OR 병기
inhos_tab <- dplyr::bind_rows(
  dplyr::bind_rows(lapply(list(
    c("cto_success", "Technical success of CTO PCI"),
    c("inhos_any", "Any in-hospital adverse event"),
    c("inhos_death", "In-hospital death"),
    c("PeriMI", "Periprocedural myocardial infarction"),
    c("inhos_CVA", "In-hospital stroke"),
    c("CIN_2", "Contrast-induced nephropathy"),
    c("Final_TIMI3", "Final TIMI 3 flow"),
    c("CR_yes", "Complete revascularization")),
    function(z) crude_or(an, z[1], z[2]))),
  inhos_tab)

# --- 관리 전략(management strategy) 지표의 보정 분석 -------------------------
# 병변 복잡도·동반질환을 보정한 뒤에도 여성에게 덜 집약적인 전략이 적용되었는지
# STRAT_VARS(6개 지표 = 단일 검정 family) 와 retro_attempt_any 는
# 00_cohort_prep.R 에 정의되어 있다. 다중비교 보정은 05_qa_sensitivity.R (A2).
strat_tab <- dplyr::bind_rows(lapply(STRAT_VARS, function(z)
  dplyr::bind_rows(crude_or(an, z[1], z[2]), logit_sex(an, z[1], z[2]))))
strat_tab$note <- "Management-strategy indicator; adjusted for Model 4 covariates (comorbidity + lesion complexity)"
inhos_tab <- dplyr::bind_rows(inhos_tab, strat_tab)
write.csv(inhos_tab[order(inhos_tab$outcome, inhos_tab$model), ],
          file.path(DIR_TAB, "table13_adjusted_procedural_inhospital.csv"),
          row.names = FALSE, na = "")

###############################################################################
# 14. 보고서용 요약 · sessionInfo
###############################################################################

key_summary <- data.frame(
  item = c("Analysis cohort N", "Female n (%)",
           "Complete-case analysis set N (Model 4)",
           "DMS events (total / male / female)",
           "Crude HR, female vs male (DMS)",
           "Age-adjusted HR (DMS)",
           "Age-adjusted HR, spline age (DMS)",
           "Fully adjusted HR, Model 4 (DMS)",
           "Fully adjusted HR, spline age (DMS)",
           "Multiple-imputation HR (DMS)",
           "IPTW HR (DMS)", "PSM HR (DMS)",
           "E-value, unadjusted DMS (point / CI)",
           "E-value, fully adjusted DMS (point / CI)",
           "Fully adjusted HR, all-cause death",
           "Fully adjusted HR, non-CV death",
           "Fully adjusted HR, CV death",
           "Fully adjusted HR, MACCE (published def.)",
           "Fully adjusted HR, any revascularization"),
  value = c(
    as.character(nrow(co)),
    sprintf("%d (%.1f%%)", sum(co$female), 100 * mean(co$female)),
    sprintf("%d (%.1f%%)", nrow(an), 100 * nrow(an) / nrow(co)),
    sprintf("%d / %d / %d", sum(d_dms$.event), sum(d_dms$.event[d_dms$female == 0]),
            sum(d_dms$.event[d_dms$female == 1])),
    paste0(staged$HR_fmt[1], ", p ", staged$p_fmt[1]),
    paste0(staged$HR_fmt[2], ", p ", staged$p_fmt[2]),
    paste0(staged_spline$HR_fmt[1], ", p ", staged_spline$p_fmt[1]),
    paste0(staged$HR_fmt[4], ", p ", staged$p_fmt[4]),
    paste0(staged_spline$HR_fmt[2], ", p ", staged_spline$p_fmt[2]),
    paste0(mi_tab$HR_fmt, ", p ", mi_tab$p_fmt),
    paste0(causal_tab$HR_fmt[causal_tab$outcome_key == "DMS" &
                               causal_tab$method == "IPTW (standardization weights, ATE)"], ", p ",
           causal_tab$p_fmt[causal_tab$outcome_key == "DMS" &
                              causal_tab$method == "IPTW (standardization weights, ATE)"]),
    paste0(causal_tab$HR_fmt[causal_tab$outcome_key == "DMS" &
                               grepl("^1:1", causal_tab$method)], ", p ",
           causal_tab$p_fmt[causal_tab$outcome_key == "DMS" & grepl("^1:1", causal_tab$method)]),
    local({ r <- ev_tab[ev_tab$method == "Unadjusted (for comparison)" &
                          grepl("^DMS", ev_tab$outcome), ]
            sprintf("%.2f / %.2f", r$evalue_point[1], r$evalue_ci[1]) }),
    local({ r <- ev_tab[ev_tab$method == "Fully adjusted Cox (Model 4)" &
                          grepl("^DMS", ev_tab$outcome), ]
            sprintf("%.2f / %.2f", r$evalue_point[1], r$evalue_ci[1]) }),
    adj_all$HR_fmt[adj_all$outcome_key == "Death"     & adj_all$model == "Fully adjusted"],
    adj_all$HR_fmt[adj_all$outcome_key == "NCVdeath"  & adj_all$model == "Fully adjusted"],
    adj_all$HR_fmt[adj_all$outcome_key == "CVdeath"   & adj_all$model == "Fully adjusted"],
    adj_all$HR_fmt[adj_all$outcome_key == "MACCE_pub" & adj_all$model == "Fully adjusted"],
    adj_all$HR_fmt[adj_all$outcome_key == "AnyRevasc" & adj_all$model == "Fully adjusted"]),
  stringsAsFactors = FALSE)
write.csv(key_summary, file.path(DIR_TAB, "table6d_key_summary.csv"), row.names = FALSE, na = "")

# 가정 점검 보고서 (텍스트)
capture.output({
  cat("K-CTO sex difference analysis — assumption checks and E-value statements\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Cohort:", COHORT_LABEL, "\n\n")
  cat("--- 1. Age linearity (Cox, primary outcome DMS) ---\n")
  print(lr_age)
  cat("\nLR test p for non-linearity of age:", fmt_p(p_age_nonlin), "\n\n")
  cat("--- 2. Variance inflation factors (Model 4) ---\n"); print(vifs)
  cat("\n--- 3. Proportional hazards (cox.zph, Model 4) ---\n"); print(ph_tab)
  cat("\n--- 4. Time-varying female HR by interval ---\n"); print(tv_tab)
  cat("\n--- 4b. Staged adjustment (primary outcome) ---\n")
  print(staged[, c("model", "n", "n_events", "HR_fmt", "p_fmt")])
  print(staged_spline[, c("model", "n", "n_events", "HR_fmt", "p_fmt")])
  print(attenuation)
  cat("\n--- 4c. Age-stratified 5-year cumulative incidence (DMS) ---\n"); print(age_strat)
  cat("\n--- 4d. Strength of the age association (E-value benchmark) ---\n"); print(age_strength)
  cat("\n--- 5. RMST at 5 years ---\n"); print(as.data.frame(rmst_tab))
  cat("\n--- 6. Propensity/standardization weight diagnostics ---\n"); print(ps_diag)
  cat("\n--- 7. PSM diagnostics ---\n"); print(psm_diag)
  cat("\n--- 8. E-values ---\n")
  for (i in seq_len(nrow(ev_tab)))
    cat(sprintf("[%s | %s] E-value point = %.2f, CI limit = %.2f\n  %s\n",
                ev_tab$outcome[i], ev_tab$method[i], ev_tab$evalue_point[i],
                ev_tab$evalue_ci[i], ev_tab$interpretation[i]))
  cat("\n--- 9. Multiple imputation ---\n"); print(mi_tab)
}, file = file.path(DIR_REPX, "assumption_checks_04_primary.txt"))

capture.output({
  cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Script: scripts/04_primary_analysis.R\n")
  cat("Cohort:", COHORT_LABEL, "\n")
  cat("set.seed(42)\n\n")
  print(sessionInfo())
}, file = file.path(DIR_REPX, "sessionInfo_04_primary.txt"))

message("== DONE (04_primary_analysis.R) ==")
