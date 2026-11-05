###############################################################################
# K-CTO Registry 2.0 — Sex Difference Analysis
# Phase 2 (descriptive): Table 1-4 by sex + 미보정 생존분석
#
# 실행:  Rscript scripts/02_descriptive_by_sex.R
# 입력:  data/processed/kcto_sex_analysis.rds  (Phase 1 superset, 9,316 x 268)
# 출력:  output/tables/*.csv, output/tables/tables1to4_by_sex.docx
#
# 코호트: scripts/00_cohort_prep.R 의 build_cohort()
#         flag_cto_pci & index_yr 2014-2023 -> N = 6,777 (CLAUDE.md §1 확정)
# Endpoint: derive_endpoints() 파생 DMS / MACCE_pub / MACE_pub / NACE_pub
#           ⛔ 기존 MACE / NACE 컬럼 직접 사용 금지 (CLAUDE.md §1-1)
#
# 이 스크립트에서 하지 않는 것: 보정 다변량 모델 -> scripts/04_primary_analysis.R
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
  library(MASS)      # ginv (다범주 SMD)
})

set.seed(42)

PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "00_cohort_prep.R"))

DIR_TAB  <- file.path(PROJ, "output", "tables")
DIR_REPX <- file.path(PROJ, "interim_report", "phase2_stats", "reports")
dir.create(DIR_TAB,  recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_REPX, recursive = TRUE, showWarnings = FALSE)

co <- build_cohort()

###############################################################################
# 1. 통계 helper
###############################################################################

fmt_num <- function(x, d = 1) formatC(x, format = "f", digits = d)

# --- SMD ---------------------------------------------------------------------
smd_continuous <- function(x, g) {
  # 부호: Female - Male (g 의 2번째 수준 - 1번째 수준)
  lv <- levels(g)
  x1 <- x[g == lv[1] & !is.na(x)]; x2 <- x[g == lv[2] & !is.na(x)]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  s <- sqrt((stats::var(x1) + stats::var(x2)) / 2)
  if (!is.finite(s) || s == 0) return(NA_real_)
  (mean(x2) - mean(x1)) / s
}

smd_categorical <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- droplevels(factor(x[ok])); g <- droplevels(factor(g[ok]))
  lv <- levels(g)
  if (length(lv) != 2 || nlevels(x) < 2) return(NA_real_)
  p1 <- as.numeric(prop.table(table(x[g == lv[1]])))
  p2 <- as.numeric(prop.table(table(x[g == lv[2]])))
  k <- length(p1)
  if (k == 2) {
    # 이진: 부호 유지 (Female - Male)
    den <- sqrt((p1[2] * (1 - p1[2]) + p2[2] * (1 - p2[2])) / 2)
    if (!is.finite(den) || den == 0) return(NA_real_)
    return((p2[2] - p1[2]) / den)
  }
  T1 <- p1[-k]; T2 <- p2[-k]
  S1 <- diag(p1[-k], k - 1) - tcrossprod(p1[-k])
  S2 <- diag(p2[-k], k - 1) - tcrossprod(p2[-k])
  S  <- (S1 + S2) / 2
  dvec <- T1 - T2
  out <- tryCatch(as.numeric(sqrt(t(dvec) %*% MASS::ginv(S) %*% dvec)),
                  error = function(e) NA_real_)
  if (is.na(out) || !is.finite(out)) NA_real_ else out
}

# --- 정규성 판정: n>5000 이면 seed 고정 5,000 표본으로 Shapiro-Wilk ----------
shapiro_flex <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(NA_real_)
  if (length(unique(x)) < 3) return(0)
  if (length(x) > 5000) x <- sample(x, 5000)
  tryCatch(stats::shapiro.test(x)$p.value, error = function(e) NA_real_)
}

###############################################################################
# 2. 요약 함수 (long format 행 생성)
###############################################################################

EMPTY_ROW <- function(n = 1) {
  data.frame(
    table = NA_character_, section = NA_character_, variable = NA_character_,
    label = NA_character_, level = NA_character_, type = NA_character_,
    denominator = NA_character_,
    n_overall = NA_integer_, n_male = NA_integer_, n_female = NA_integer_,
    n_missing = NA_integer_, pct_missing = NA_real_, missing_type = NA_character_,
    overall = NA_character_, male = NA_character_, female = NA_character_,
    p_value = NA_real_, p_fmt = NA_character_, p_method = NA_character_,
    smd = NA_real_, smd_abs = NA_real_,
    overall_mean_sd = NA_character_, male_mean_sd = NA_character_,
    female_mean_sd = NA_character_,
    overall_median_iqr = NA_character_, male_median_iqr = NA_character_,
    female_median_iqr = NA_character_,
    p_ttest = NA_real_, p_wilcox = NA_real_, normality = NA_character_,
    stringsAsFactors = FALSE)[rep(1, n), , drop = FALSE]
}

sum_continuous <- function(dat, var, label, tbl, section,
                           denom = "Full analysis cohort", d = 1,
                           missing_type = "sporadic") {
  x <- dat[[var]]; g <- dat$sex_f
  ok <- !is.na(x)
  xm <- x[ok & g == "Male"]; xf <- x[ok & g == "Female"]

  ms <- function(v) if (length(v) == 0) NA_character_ else
    paste0(fmt_num(mean(v), d), " +/- ", fmt_num(stats::sd(v), d))
  mi <- function(v) if (length(v) == 0) NA_character_ else
    paste0(fmt_num(stats::median(v), d), " [",
           fmt_num(stats::quantile(v, .25), d), ", ",
           fmt_num(stats::quantile(v, .75), d), "]")

  p_sw_m <- shapiro_flex(xm); p_sw_f <- shapiro_flex(xf)
  normal <- (!is.na(p_sw_m) && p_sw_m > 0.05) && (!is.na(p_sw_f) && p_sw_f > 0.05)

  p_t <- tryCatch(stats::t.test(xm, xf)$p.value, error = function(e) NA_real_)
  p_w <- tryCatch(suppressWarnings(
    stats::wilcox.test(xm, xf, exact = FALSE)$p.value), error = function(e) NA_real_)

  r <- EMPTY_ROW(1)
  r$table <- tbl; r$section <- section; r$variable <- var
  r$label <- label; r$level <- ""; r$type <- "continuous"
  r$denominator <- denom
  r$n_overall <- sum(ok); r$n_male <- length(xm); r$n_female <- length(xf)
  r$n_missing <- sum(!ok); r$pct_missing <- round(100 * mean(!ok), 1)
  r$missing_type <- missing_type
  r$overall_mean_sd <- ms(x[ok]); r$male_mean_sd <- ms(xm); r$female_mean_sd <- ms(xf)
  r$overall_median_iqr <- mi(x[ok]); r$male_median_iqr <- mi(xm); r$female_median_iqr <- mi(xf)
  r$overall <- if (normal) r$overall_mean_sd else r$overall_median_iqr
  r$male    <- if (normal) r$male_mean_sd    else r$male_median_iqr
  r$female  <- if (normal) r$female_mean_sd  else r$female_median_iqr
  r$p_value <- if (normal) p_t else p_w
  r$p_method <- if (normal) "Welch t-test" else "Wilcoxon rank-sum"
  r$p_fmt <- fmt_p(r$p_value)
  r$p_ttest <- p_t; r$p_wilcox <- p_w
  r$normality <- if (normal) "normal (Shapiro-Wilk p>0.05 both groups)" else "non-normal"
  r$smd <- smd_continuous(x, g); r$smd_abs <- abs(r$smd)
  r
}

sum_binary <- function(dat, var, label, tbl, section,
                       denom = "Full analysis cohort", missing_type = "sporadic",
                       yes_label = NULL) {
  x <- dat[[var]]; g <- dat$sex_f
  f <- factor(ifelse(is.na(x), NA, ifelse(x == 1, "Yes", "No")), levels = c("No", "Yes"))
  ok <- !is.na(f)
  tab <- table(g[ok], f[ok])
  n_m <- sum(ok & g == "Male"); n_f <- sum(ok & g == "Female")
  yes_m <- sum(ok & g == "Male" & f == "Yes"); yes_f <- sum(ok & g == "Female" & f == "Yes")
  yes_o <- sum(ok & f == "Yes")

  p <- NA_real_; meth <- NA_character_
  if (all(dim(tab) == c(2, 2)) && sum(tab) > 0) {
    ex <- tryCatch(suppressWarnings(stats::chisq.test(tab)$expected), error = function(e) NULL)
    if (!is.null(ex) && min(ex) >= 5) {
      p <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$p.value)
      meth <- "Pearson chi-square"
    } else {
      p <- tryCatch(stats::fisher.test(tab)$p.value, error = function(e) NA_real_)
      meth <- "Fisher's exact"
    }
  }

  r <- EMPTY_ROW(1)
  r$table <- tbl; r$section <- section; r$variable <- var
  r$label <- if (is.null(yes_label)) label else paste0(label, ", ", yes_label)
  r$level <- ""; r$type <- "binary"; r$denominator <- denom
  r$n_overall <- sum(ok); r$n_male <- n_m; r$n_female <- n_f
  r$n_missing <- sum(!ok); r$pct_missing <- round(100 * mean(!ok), 1)
  r$missing_type <- missing_type
  r$overall <- sprintf("%d (%.1f)", yes_o, 100 * yes_o / max(sum(ok), 1))
  r$male    <- sprintf("%d (%.1f)", yes_m, 100 * yes_m / max(n_m, 1))
  r$female  <- sprintf("%d (%.1f)", yes_f, 100 * yes_f / max(n_f, 1))
  r$p_value <- p; r$p_fmt <- fmt_p(p); r$p_method <- meth
  r$smd <- smd_categorical(f, g); r$smd_abs <- abs(r$smd)
  r
}

sum_categorical <- function(dat, var, label, tbl, section,
                            denom = "Full analysis cohort", missing_type = "sporadic") {
  f <- dat[[var]]; g <- dat$sex_f
  if (!is.factor(f)) f <- factor(f)
  ok <- !is.na(f)
  tab <- table(f[ok], g[ok])
  n_m <- sum(ok & g == "Male"); n_f <- sum(ok & g == "Female")

  p <- NA_real_; meth <- NA_character_
  if (nrow(tab) >= 2 && ncol(tab) == 2) {
    ex <- tryCatch(suppressWarnings(stats::chisq.test(tab)$expected), error = function(e) NULL)
    if (!is.null(ex) && min(ex) >= 5) {
      p <- suppressWarnings(stats::chisq.test(tab)$p.value); meth <- "Pearson chi-square"
    } else {
      p <- tryCatch(stats::fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value,
                    error = function(e) NA_real_)
      meth <- "Fisher's exact (Monte Carlo, B=1e5)"
    }
  }

  hdr <- EMPTY_ROW(1)
  hdr$table <- tbl; hdr$section <- section; hdr$variable <- var
  hdr$label <- label; hdr$level <- ""; hdr$type <- "categorical"
  hdr$denominator <- denom
  hdr$n_overall <- sum(ok); hdr$n_male <- n_m; hdr$n_female <- n_f
  hdr$n_missing <- sum(!ok); hdr$pct_missing <- round(100 * mean(!ok), 1)
  hdr$missing_type <- missing_type
  hdr$p_value <- p; hdr$p_fmt <- fmt_p(p); hdr$p_method <- meth
  hdr$smd <- smd_categorical(f, g); hdr$smd_abs <- abs(hdr$smd)

  lvs <- levels(droplevels(f[ok]))
  body <- EMPTY_ROW(length(lvs))
  for (i in seq_along(lvs)) {
    L <- lvs[i]
    n_o <- sum(ok & f == L); nm <- sum(ok & f == L & g == "Male")
    nf <- sum(ok & f == L & g == "Female")
    body$table[i] <- tbl; body$section[i] <- section; body$variable[i] <- var
    body$label[i] <- paste0("  ", L); body$level[i] <- L
    body$type[i] <- "categorical_level"; body$denominator[i] <- denom
    body$overall[i] <- sprintf("%d (%.1f)", n_o, 100 * n_o / max(sum(ok), 1))
    body$male[i]    <- sprintf("%d (%.1f)", nm, 100 * nm / max(n_m, 1))
    body$female[i]  <- sprintf("%d (%.1f)", nf, 100 * nf / max(n_f, 1))
  }
  rbind(hdr, body)
}

###############################################################################
# 3. Table 1 — Baseline characteristics by sex
###############################################################################

build_table1 <- function(dat, tbl = "Table 1") {
  out <- list(); A <- function(x) out[[length(out) + 1]] <<- x

  sec <- "Demographics"
  A(sum_continuous(dat, "Age", "Age, years", tbl, sec, d = 1))
  A(sum_binary(dat, "age75", "Age >=75 years", tbl, sec, yes_label = "n (%)"))
  A(sum_continuous(dat, "Height", "Height, cm", tbl, sec, d = 1))
  A(sum_continuous(dat, "Weight", "Weight, kg", tbl, sec, d = 1))
  A(sum_continuous(dat, "BMI_re", "Body mass index, kg/m2", tbl, sec, d = 1))
  A(sum_continuous(dat, "bsa", "Body surface area, m2", tbl, sec, d = 2))

  sec <- "Risk factors and comorbidities"
  A(sum_binary(dat, "Current_sm", "Current smoker", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "HTN", "Hypertension", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "Dyslipidemia", "Dyslipidemia", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "DM", "Diabetes mellitus", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "DM_Insulin", "Insulin-treated diabetes", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "AF", "Atrial fibrillation", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "CKD_2", "Chronic kidney disease (eGFR <60)", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "esrd_eq", "ESRD-equivalent (eGFR <15)", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "prePCI", "Previous PCI", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "preCABG", "Previous CABG", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "preMI", "Previous myocardial infarction", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "preCVA", "Previous stroke", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "prePVD", "Previous peripheral artery disease", tbl, sec, yes_label = "n (%)"))

  sec <- "Clinical presentation"
  A(sum_categorical(dat, "Clindx_f", "Clinical diagnosis, n (%)", tbl, sec))
  A(sum_binary(dat, "stable_angina", "Stable angina", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "acute_MI", "Myocardial infarction at presentation", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "acs", "Acute coronary syndrome at presentation", tbl, sec, yes_label = "n (%)"))

  sec <- "Laboratory and imaging"
  A(sum_continuous(dat, "Cr", "Serum creatinine, mg/dL", tbl, sec, d = 2))
  A(sum_continuous(dat, "egfr", "eGFR, mL/min/1.73m2", tbl, sec, d = 1))
  A(sum_continuous(dat, "preEF", "LV ejection fraction, %", tbl, sec, d = 1))
  A(sum_binary(dat, "lvef40", "LV ejection fraction <40%", tbl, sec, yes_label = "n (%)"))
  A(sum_continuous(dat, "CKMB_pre", "Pre-procedural CK-MB, ng/mL", tbl, sec, d = 1))

  sec <- "Index procedure setting"
  A(sum_continuous(dat, "index_yr", "Index PCI year", tbl, sec, d = 0))
  A(sum_categorical(dat, "index_period_f", "Index PCI period, n (%)", tbl, sec))
  A(sum_categorical(dat, "Center_volume_f", "Center volume (registry contribution), n (%)", tbl, sec))
  A(sum_categorical(dat, "hos_location_f", "Center location, n (%)", tbl, sec))

  sec <- "Follow-up"
  A(sum_continuous(dat, "fudur2", "Follow-up duration, days", tbl, sec, d = 0))
  A(sum_binary(dat, "fu_zero", "No follow-up beyond index date (fudur = 0)", tbl, sec, yes_label = "n (%)"))

  do.call(rbind, out)
}

###############################################################################
# 4. Table 2 — Coronary lesion / angiographic characteristics by sex
###############################################################################

build_table2 <- function(dat, tbl = "Table 2") {
  out <- list(); A <- function(x) out[[length(out) + 1]] <<- x

  sec <- "Extent of coronary disease"
  A(sum_categorical(dat, "Vessel.Disease_f", "Number of diseased vessels, n (%)", tbl, sec))
  A(sum_binary(dat, "multivessel", "Multivessel disease (>=2 vessels)", tbl, sec, yes_label = "n (%)"))
  A(sum_categorical(dat, "CTO_number_f", "Number of CTO vessels, n (%)", tbl, sec))

  sec <- "CTO lesion location"
  A(sum_categorical(dat, "CTO_Vessel_f", "CTO vessel, n (%)", tbl, sec))
  A(sum_categorical(dat, "CTO_location_f", "CTO segment, n (%)", tbl, sec))
  A(sum_binary(dat, "CTO_location_pm", "Proximal or mid CTO segment", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "LM.or.LAD.OS", "Ostial LAD or LM CTO", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "ISR_CTO", "In-stent restenotic CTO", tbl, sec, yes_label = "n (%)"))

  sec <- "Collateral circulation"
  A(sum_categorical(dat, "number_donor_vessel_f", "Number of collateral donor vessels, n (%)", tbl, sec))
  A(sum_binary(dat, "lad_donor", "LAD as donor vessel", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "lcx_donor", "LCX as donor vessel", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "rca_donor", "RCA as donor vessel", tbl, sec, yes_label = "n (%)"))
  A(sum_categorical(dat, "Rentrop_f", "Rentrop collateral grade, n (%)", tbl, sec))
  A(sum_binary(dat, "Rentrop3", "Rentrop grade 3", tbl, sec, yes_label = "n (%)"))

  sec <- "Lesion complexity (J-CTO)"
  A(sum_continuous(dat, "JCTO_re", "J-CTO score", tbl, sec, d = 1))
  A(sum_categorical(dat, "JCTO_g4_f", "J-CTO score, n (%)", tbl, sec))
  A(sum_binary(dat, "JCTO_bi", "J-CTO score >=2", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "Blunt", "Blunt stump", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "Calcific", "Heavy calcification", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "Tortuocity", "Severe tortuosity", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "Length20", "Occlusion length >=20 mm", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "Retry", "Retry (previously attempted) lesion", tbl, sec, yes_label = "n (%)"))

  sec <- "Non-CTO disease"
  A(sum_binary(dat, "non_CTO.PCI", "PCI of non-CTO lesion(s) at index", tbl, sec, yes_label = "n (%)"))

  do.call(rbind, out)
}

###############################################################################
# 5. Table 3 — Procedural characteristics and discharge medication by sex
#    구조적 결측 변수는 적용 대상 하위집단으로 분모를 제한하고 denominator 명시
###############################################################################

build_table3 <- function(dat, tbl = "Table 3") {
  out <- list(); A <- function(x) out[[length(out) + 1]] <<- x

  sec <- "Procedural outcome"
  A(sum_binary(dat, "cto_success", "Technical success of CTO PCI", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "Final_TIMI3", "Final TIMI 3 flow", tbl, sec, yes_label = "n (%)"))
  A(sum_categorical(dat, "Final_TIMI_f", "Final TIMI flow, n (%)", tbl, sec))
  A(sum_binary(dat, "CR_yes", "Complete revascularization (recoded: codebook Yes=1/No=2)",
               tbl, sec, yes_label = "n (%)"))

  sec <- "Procedural strategy"
  A(sum_binary(dat, "Bilati_inj", "Bilateral injection", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "IVUS", "IVUS use", tbl, sec, yes_label = "n (%)"))
  A(sum_categorical(dat, "approach_direction_f", "Approach direction, n (%)", tbl, sec))
  A(sum_binary(dat, "ante_approach_only", "Antegrade approach only", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "AWE", "Antegrade wire escalation", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "ADR", "Antegrade dissection re-entry", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "RWE", "Retrograde wire escalation", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "rCART", "Reverse CART", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "CART", "CART", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "Fail_tech", "Failure of initial wiring technique", tbl, sec, yes_label = "n (%)"))
  A(sum_categorical(dat, "Colla_approach_f", "Collateral channel approach, n (%)", tbl, sec))
  A(sum_categorical(dat, "Retro_attempt_f", "Retrograde attempt, n (%)", tbl, sec))
  A(sum_binary(dat, "Retro_first", "Retrograde-first strategy", tbl, sec, yes_label = "n (%)"))

  # 구조적 결측: Retro_success 는 retrograde 시도 환자에게만 적용
  sub_retro <- dat[dat$Retro_attempt_re %in% c(1, 2), , drop = FALSE]
  A(sum_binary(sub_retro, "Retro_success", "Retrograde success", tbl, sec,
               denom = "Patients with a retrograde attempt (Retro_attempt_re = 1 or 2)",
               missing_type = "structural (restricted denominator)", yes_label = "n (%)"))

  sec <- "Lesion treatment"
  A(sum_binary(dat, "TV_DEB", "Drug-coated balloon", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "TV_POBA", "Balloon angioplasty only (POBA)", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "TV_STENT", "Stent implantation", tbl, sec, yes_label = "n (%)"))

  # 구조적 결측: stent 상세는 stent 를 넣은 환자에게만 적용
  sub_stent <- dat[dat$TV_STENT %in% 1, , drop = FALSE]
  dn <- "Patients receiving >=1 stent (TV_STENT = 1)"
  mt <- "structural (restricted denominator)"
  A(sum_categorical(sub_stent, "n_stent_f", "Number of stents, n (%)", tbl, sec, denom = dn, missing_type = mt))
  A(sum_continuous(sub_stent, "stent_dia_mean", "Mean stent diameter, mm", tbl, sec, denom = dn, d = 2, missing_type = mt))
  A(sum_continuous(sub_stent, "stent_dia_min", "Minimum stent diameter, mm", tbl, sec, denom = dn, d = 2, missing_type = mt))
  A(sum_continuous(sub_stent, "stent_dia_max", "Maximum stent diameter, mm", tbl, sec, denom = dn, d = 2, missing_type = mt))
  A(sum_continuous(sub_stent, "stent_len_total", "Total stent length, mm", tbl, sec, denom = dn, d = 1, missing_type = mt))

  sec <- "Procedural resource use"
  A(sum_categorical(dat, "n_wire_f", "Number of guidewires, n (%)", tbl, sec,
                    missing_type = "sporadic (not recorded)"))
  A(sum_categorical(dat, "n_micro_f", "Number of microcatheters, n (%)", tbl, sec,
                    missing_type = "sporadic (not recorded)"))
  A(sum_continuous(dat, "fluorotime", "Fluoroscopy time, min", tbl, sec, d = 1,
                   missing_type = "sporadic (not recorded)"))
  A(sum_continuous(dat, "dye", "Contrast volume, mL", tbl, sec, d = 0,
                   missing_type = "sporadic (not recorded)"))

  sec <- "Discharge medication"
  A(sum_binary(dat, "dc_asp", "Aspirin", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "dc_clopid", "Clopidogrel", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "dc_ticagrelor", "Ticagrelor", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "dc_prasugrel", "Prasugrel", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "dc_cilostazol", "Cilostazol", tbl, sec, yes_label = "n (%)"))
  A(sum_categorical(dat, "dc_antiplt2_f", "Antiplatelet regimen at discharge, n (%)", tbl, sec))
  sub_dapt <- dat[dat$dc_antiplt_2 %in% c("ac", "ap", "at"), , drop = FALSE]
  A(sum_binary(sub_dapt, "dapt_potent", "Potent P2Y12 inhibitor (prasugrel or ticagrelor)", tbl, sec,
               denom = "Patients discharged on dual antiplatelet therapy (ac/ap/at)",
               missing_type = "structural (restricted denominator)", yes_label = "n (%)"))
  A(sum_binary(dat, "dc_statin", "Statin", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "dc_BB", "Beta-blocker", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "dc_ACE_ARB", "ACE inhibitor or ARB", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "dc_entresto", "ARNI", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "dc_CCB", "Calcium channel blocker", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "dc_nitro", "Nitrate", tbl, sec, yes_label = "n (%)"))

  do.call(rbind, out)
}

###############################################################################
# 6. Table 4 — In-hospital outcomes by sex
###############################################################################

build_table4 <- function(dat, tbl = "Table 4") {
  out <- list(); A <- function(x) out[[length(out) + 1]] <<- x

  sec <- "In-hospital clinical events"
  A(sum_binary(dat, "inhos_any", "Any in-hospital adverse event (composite)", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "inhos_death", "In-hospital death", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "PeriMI", "Periprocedural myocardial infarction", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "inhos_CVA", "In-hospital stroke", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "urgent_PCI", "Urgent repeat PCI", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "urgent_CABG", "Urgent CABG", tbl, sec, yes_label = "n (%)"))
  A(sum_binary(dat, "urgent_tampon_intervention", "Urgent tamponade intervention", tbl, sec, yes_label = "n (%)"))

  sec <- "Contrast-associated renal injury"
  A(sum_binary(dat, "CIN_2", "Contrast-induced nephropathy", tbl, sec,
               denom = "Patients with both baseline and day-2 creatinine",
               missing_type = "sporadic (day-2 creatinine not measured)", yes_label = "n (%)"))
  A(sum_continuous(dat, "Cr_2d", "Day-2 serum creatinine, mg/dL", tbl, sec, d = 2,
                   denom = "Patients with day-2 creatinine",
                   missing_type = "sporadic (not measured)"))
  A(sum_continuous(dat, "cr_increase", "Absolute creatinine increase, mg/dL", tbl, sec, d = 2,
                   denom = "Patients with baseline and day-2 creatinine",
                   missing_type = "sporadic (not measured)"))

  do.call(rbind, out)
}

###############################################################################
# 7. 미보정 생존분석
#    ⚠️ endpoint 는 전부 derive_endpoints() 파생본 (CLAUDE.md §1-1)
#    day-0 event 처리: event==1 & fudur==0 -> 0.5일 대체 (주 규칙)
###############################################################################

OUTCOMES <- list(
  list(key = "DMS",       ev = "DMS",            tm = "DMS_fudur",
       lab = "PRIMARY: DMS (all-cause death, MI, or stroke)"),
  list(key = "MACCE_pub", ev = "MACCE_pub",      tm = "MACCE_pub_fudur",
       lab = "MACCE, published definition (death, MI, stroke, any revascularization)"),
  list(key = "MACE_pub",  ev = "MACE_pub",       tm = "MACE_pub_fudur",
       lab = "MACE, published definition (death, MI, any revascularization)"),
  list(key = "NACE_pub",  ev = "NACE_pub",       tm = "NACE_pub_fudur",
       lab = "NACE, published definition (MACCE + TIMI major bleeding)"),
  list(key = "Death",     ev = "Death",          tm = "death_fudur", lab = "All-cause death"),
  list(key = "CVdeath",   ev = "CVdeath",        tm = "death_fudur", lab = "Cardiovascular death"),
  list(key = "NCVdeath",  ev = "NCVdeath",       tm = "death_fudur", lab = "Non-cardiovascular death"),
  list(key = "AnyRevasc", ev = "AnyRevasc",      tm = "AnyRevasc_fudur",
       lab = "Any coronary revascularization (TVR, NLR, or CABG)"),
  list(key = "MI",        ev = "MI",             tm = "MI_fudur",    lab = "Myocardial infarction"),
  list(key = "Stroke",    ev = "CVA_re",         tm = "CVA_fudur",   lab = "Stroke (CVA)"),
  list(key = "TVR",       ev = "TVR_re2",        tm = "TVR_fudur2",  lab = "Target-vessel revascularization"),
  list(key = "NLR",       ev = "NLR_re2",        tm = "NLR_fudur2",  lab = "Non-target-lesion revascularization"),
  list(key = "CABG",      ev = "CABG_re",        tm = "CABG_fudur",  lab = "Coronary artery bypass grafting"),
  list(key = "MajBld",    ev = "Bleeding_major", tm = "Bleeding_major_fudur", lab = "TIMI major bleeding")
)

TIMES <- c(365, 1095, 1825)
TIME_LAB <- c("1-year", "3-year", "5-year")

km_one <- function(dat, ev, tm, lab, day0 = c("half_day", "exclude"), cohort_label) {
  day0 <- match.arg(day0)
  a <- apply_day0(dat[[tm]], dat[[ev]], rule = day0)
  keep <- a$keep & !is.na(dat$sex_f)
  n_day0 <- sum(!is.na(dat[[tm]]) & !is.na(dat[[ev]]) & dat[[ev]] == 1 & dat[[tm]] == 0)
  d <- data.frame(t = a$time[keep], e = a$event[keep], g = droplevels(dat$sex_f[keep]))

  fit  <- survival::survfit(survival::Surv(t, e) ~ g, data = d, conf.type = "log-log")
  sd_  <- survival::survdiff(survival::Surv(t, e) ~ g, data = d)
  p_lr <- stats::pchisq(sd_$chisq, df = length(sd_$n) - 1, lower.tail = FALSE)
  cx   <- survival::coxph(survival::Surv(t, e) ~ g, data = d)
  cs   <- summary(cx)
  ph   <- tryCatch(survival::cox.zph(cx)$table["g", "p"], error = function(x) NA_real_)

  s <- summary(fit, times = TIMES, extend = TRUE)
  strata_names <- as.character(s$strata)
  rows <- list()
  for (i in seq_along(TIMES)) {
    for (gg in c("Male", "Female")) {
      sel <- which(strata_names == paste0("g=", gg) & s$time == TIMES[i])
      if (length(sel) == 0) next
      ci <- 100 * (1 - s$surv[sel]); lo <- 100 * (1 - s$upper[sel]); hi <- 100 * (1 - s$lower[sel])
      rows[[length(rows) + 1]] <- data.frame(
        cohort = cohort_label, day0_rule = day0, outcome = lab, outcome_var = ev,
        timepoint = TIME_LAB[i], time_days = TIMES[i], sex = gg,
        n_at_risk = s$n.risk[sel], n_events_since_prev_timepoint = s$n.event[sel],
        cum_incidence_pct = round(ci, 2),
        ci_lower_pct = round(lo, 2), ci_upper_pct = round(hi, 2),
        estimate_fmt = sprintf("%.1f (%.1f-%.1f)", ci, lo, hi),
        stringsAsFactors = FALSE)
    }
  }
  ti <- do.call(rbind, rows)

  smry <- d %>% group_by(g) %>%
    summarise(n = n(), events = sum(e == 1), py = sum(t) / 365.25, .groups = "drop") %>%
    mutate(rate_per_100py = round(100 * events / py, 2))

  hd <- data.frame(
    cohort = cohort_label, day0_rule = day0, outcome = lab, outcome_var = ev,
    n_total = nrow(d),
    n_male = smry$n[smry$g == "Male"], n_female = smry$n[smry$g == "Female"],
    events_male = smry$events[smry$g == "Male"], events_female = smry$events[smry$g == "Female"],
    py_male = round(smry$py[smry$g == "Male"], 0), py_female = round(smry$py[smry$g == "Female"], 0),
    rate_male_per100py = smry$rate_per_100py[smry$g == "Male"],
    rate_female_per100py = smry$rate_per_100py[smry$g == "Female"],
    n_day0_events = n_day0,
    logrank_p = p_lr, logrank_p_fmt = fmt_p(p_lr),
    crude_HR_female_vs_male = round(cs$conf.int[1, 1], 3),
    HR_lower = round(cs$conf.int[1, 3], 3), HR_upper = round(cs$conf.int[1, 4], 3),
    HR_fmt = sprintf("%.2f (%.2f-%.2f)", cs$conf.int[1, 1], cs$conf.int[1, 3], cs$conf.int[1, 4]),
    cox_p = cs$coefficients[1, 5], cox_p_fmt = fmt_p(cs$coefficients[1, 5]),
    ph_assumption_p = round(ph, 4),
    stringsAsFactors = FALSE)

  list(timepoints = ti, summary = hd)
}

run_km_battery <- function(dat, cohort_label, day0 = "half_day", which_out = NULL) {
  sel <- if (is.null(which_out)) OUTCOMES else Filter(function(o) o$key %in% which_out, OUTCOMES)
  ti <- list(); sm <- list()
  for (o in sel) {
    r <- km_one(dat, o$ev, o$tm, o$lab, day0 = day0, cohort_label = cohort_label)
    ti[[length(ti) + 1]] <- r$timepoints
    sm[[length(sm) + 1]] <- r$summary
  }
  list(timepoints = do.call(rbind, ti), summary = do.call(rbind, sm))
}

# --- KM 곡선 좌표 (illustrator 전달용) ---------------------------------------
km_coordinates <- function(dat, ev, tm, lab, max_days = 1825) {
  a <- apply_day0(dat[[tm]], dat[[ev]], rule = "half_day")
  keep <- a$keep & !is.na(dat$sex_f)
  d <- data.frame(t = a$time[keep], e = a$event[keep], g = droplevels(dat$sex_f[keep]))
  fit <- survival::survfit(survival::Surv(t, e) ~ g, data = d, conf.type = "log-log")
  s <- summary(fit)
  out <- data.frame(
    outcome = lab, outcome_var = ev,
    sex = sub("^g=", "", as.character(s$strata)),
    time_days = s$time, n_risk = s$n.risk, n_event = s$n.event,
    cum_incidence_pct = round(100 * (1 - s$surv), 4),
    ci_lower_pct = round(100 * (1 - s$upper), 4),
    ci_upper_pct = round(100 * (1 - s$lower), 4),
    stringsAsFactors = FALSE)
  out[out$time_days <= max_days, , drop = FALSE]
}

# --- 잠재 추적기간(reverse KM) -----------------------------------------------
reverse_km <- function(dat, cohort_label) {
  d <- data.frame(t = dat$death_fudur, e = 1 - dat$Death, g = dat$sex_f)
  f <- survival::survfit(survival::Surv(t, e) ~ g, data = d)
  q <- stats::quantile(f, probs = c(0.25, 0.5, 0.75))
  data.frame(cohort = cohort_label, sex = c("Male", "Female"),
             median_potential_fu_days = as.numeric(q$quantile[, "50"]),
             q25 = as.numeric(q$quantile[, "25"]), q75 = as.numeric(q$quantile[, "75"]),
             stringsAsFactors = FALSE)
}

# --- 경쟁위험 (CV death vs non-CV death) --------------------------------------
competing_risk_cvdeath <- function(dat, cohort_label) {
  t  <- dat$death_fudur
  ev <- ifelse(dat$CVdeath == 1, 1, ifelse(dat$NCVdeath == 1, 2, 0))
  t[ev > 0 & t == 0] <- 0.5
  g <- dat$sex_f
  ok <- !is.na(t) & !is.na(ev) & !is.na(g)
  t <- t[ok]; ev <- ev[ok]; g <- droplevels(g[ok])

  ci <- cmprsk::cuminc(ftime = t, fstatus = ev, group = as.character(g), cencode = 0)
  gray_p <- ci$Tests[, "pv"]
  est <- cmprsk::timepoints(ci, times = TIMES)
  out <- list()
  for (rn in rownames(est$est)) for (j in seq_along(TIMES)) {
    out[[length(out) + 1]] <- data.frame(
      cohort = cohort_label, method = "Aalen-Johansen CIF (competing risk)",
      strata = rn, timepoint = TIME_LAB[j], time_days = TIMES[j],
      cif_pct = round(100 * est$est[rn, j], 2),
      se_pct = round(100 * sqrt(est$var[rn, j]), 2), stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, out)
  res$gray_test_p <- NA_real_
  res$gray_test_p[grepl(" 1$", res$strata)] <- gray_p["1"]
  res$gray_test_p[grepl(" 2$", res$strata)] <- gray_p["2"]
  res$gray_test_p_fmt <- fmt_p(res$gray_test_p)
  res$event_type <- ifelse(grepl(" 1$", res$strata), "CV death", "Non-CV death")
  res$sex <- ifelse(grepl("^Male", res$strata), "Male", "Female")
  res$strata <- NULL
  res[, c("cohort", "method", "sex", "event_type", "timepoint", "time_days",
          "cif_pct", "se_pct", "gray_test_p", "gray_test_p_fmt")]
}

###############################################################################
# 8. 실행 — 주 코호트
###############################################################################

message("== Tables 1-4 (main cohort, N = ", nrow(co), ") ==")
t1 <- build_table1(co); t2 <- build_table2(co)
t3 <- build_table3(co); t4 <- build_table4(co)

write.csv(t1, file.path(DIR_TAB, "table1_baseline_by_sex.csv"), row.names = FALSE, na = "")
write.csv(t2, file.path(DIR_TAB, "table2_lesion_by_sex.csv"), row.names = FALSE, na = "")
write.csv(t3, file.path(DIR_TAB, "table3_procedural_by_sex.csv"), row.names = FALSE, na = "")
write.csv(t4, file.path(DIR_TAB, "table4_inhospital_outcomes_by_sex.csv"), row.names = FALSE, na = "")

message("== Unadjusted survival analysis (main cohort) ==")
km_main <- run_km_battery(co, COHORT_LABEL, day0 = "half_day")
km_excl <- run_km_battery(co, COHORT_LABEL, day0 = "exclude")

write.csv(km_main$timepoints, file.path(DIR_TAB, "table5a_km_timepoints_by_sex.csv"),
          row.names = FALSE, na = "")
write.csv(km_main$summary, file.path(DIR_TAB, "table5b_km_summary_by_sex.csv"),
          row.names = FALSE, na = "")

# KM 곡선 좌표 (Phase 4 illustrator 전달용)
kmc <- do.call(rbind, lapply(OUTCOMES, function(o) km_coordinates(co, o$ev, o$tm, o$lab)))
write.csv(kmc, file.path(DIR_TAB, "table5e_km_curve_coordinates.csv"), row.names = FALSE, na = "")

rkm <- reverse_km(co, COHORT_LABEL)
write.csv(rkm, file.path(DIR_TAB, "table5c_potential_followup_reverseKM.csv"),
          row.names = FALSE, na = "")

cr_cv <- if (requireNamespace("cmprsk", quietly = TRUE)) competing_risk_cvdeath(co, COHORT_LABEL) else NULL
if (!is.null(cr_cv))
  write.csv(cr_cv, file.path(DIR_TAB, "table5d_competing_risk_cvdeath.csv"), row.names = FALSE, na = "")

###############################################################################
# 9. 민감도 분석 — 등록기간 2014-2022 제한 (parent Methods 기술 정합)
#    ⛔ 사전 계획을 시사하는 표현 금지 (CLAUDE.md §1-2) — 본 분석은 사전등록되지 않았다.
###############################################################################

co_win <- co[co$index_yr <= 2022, , drop = FALSE]
message("== Sensitivity cohort (2014-2022): N = ", nrow(co_win), " ==")

s1 <- build_table1(co_win, tbl = "Table S1 (2014-2022)")
s2 <- build_table2(co_win, tbl = "Table S1 (2014-2022)")
s3 <- build_table3(co_win, tbl = "Table S1 (2014-2022)")
s4 <- build_table4(co_win, tbl = "Table S1 (2014-2022)")
write.csv(rbind(s1, s2, s3, s4),
          file.path(DIR_TAB, "tableS1_sensitivity_2014_2022_tables1to4.csv"),
          row.names = FALSE, na = "")

km_win <- run_km_battery(co_win, "Sensitivity cohort 2014-2022", day0 = "half_day")
write.csv(km_win$timepoints, file.path(DIR_TAB, "tableS2_km_timepoints_2014_2022.csv"),
          row.names = FALSE, na = "")
write.csv(km_win$summary, file.path(DIR_TAB, "tableS3_km_summary_2014_2022.csv"),
          row.names = FALSE, na = "")

cmp_key <- rbind(t1, t2, t3, t4) %>%
  filter(type %in% c("continuous", "binary", "categorical")) %>%
  dplyr::select(variable, label, main_male = male, main_female = female,
                main_p = p_fmt, main_smd = smd)
cmp_sen <- rbind(s1, s2, s3, s4) %>%
  filter(type %in% c("continuous", "binary", "categorical")) %>%
  dplyr::select(variable, sens_male = male, sens_female = female,
                sens_p = p_fmt, sens_smd = smd)
cmp <- dplyr::left_join(cmp_key, cmp_sen, by = "variable") %>%
  mutate(smd_delta = round(sens_smd - main_smd, 3),
         direction_consistent = ifelse(is.na(main_smd) | is.na(sens_smd), NA,
                                       sign(main_smd) == sign(sens_smd)))
write.csv(cmp, file.path(DIR_TAB, "tableS4_main_vs_2014_2022_comparison.csv"),
          row.names = FALSE, na = "")

km_cmp <- dplyr::left_join(
  km_main$summary %>% dplyr::select(outcome_var, outcome, main_HR = HR_fmt, main_lr_p = logrank_p_fmt),
  km_win$summary  %>% dplyr::select(outcome_var, sens_HR = HR_fmt, sens_lr_p = logrank_p_fmt),
  by = "outcome_var")
add5 <- function(tp, pref) {
  tp %>% filter(time_days == 1825) %>%
    dplyr::select(outcome_var, sex, estimate_fmt) %>%
    tidyr::pivot_wider(names_from = sex, values_from = estimate_fmt,
                       names_prefix = paste0(pref, "_5yr_"))
}
km_cmp <- km_cmp %>%
  dplyr::left_join(add5(km_main$timepoints, "main"), by = "outcome_var") %>%
  dplyr::left_join(add5(km_win$timepoints, "sens"), by = "outcome_var")
write.csv(km_cmp, file.path(DIR_TAB, "tableS4b_km_main_vs_2014_2022.csv"), row.names = FALSE, na = "")

# --- day-0 처리 민감도 --------------------------------------------------------
d0 <- dplyr::left_join(
  km_main$summary %>% dplyr::select(outcome_var, outcome, half_HR = HR_fmt,
                                    half_lr_p = logrank_p_fmt, n_day0_events),
  km_excl$summary %>% dplyr::select(outcome_var, excl_HR = HR_fmt, excl_lr_p = logrank_p_fmt),
  by = "outcome_var")
d0 <- dplyr::left_join(d0,
  km_main$timepoints %>% filter(time_days == 1825) %>%
    dplyr::select(outcome_var, sex, half_5yr = estimate_fmt) %>%
    tidyr::pivot_wider(names_from = sex, values_from = half_5yr, names_prefix = "half_5yr_"),
  by = "outcome_var")
d0 <- dplyr::left_join(d0,
  km_excl$timepoints %>% filter(time_days == 1825) %>%
    dplyr::select(outcome_var, sex, excl_5yr = estimate_fmt) %>%
    tidyr::pivot_wider(names_from = sex, values_from = excl_5yr, names_prefix = "excl_5yr_"),
  by = "outcome_var")
write.csv(d0, file.path(DIR_TAB, "tableS5_day0_event_handling_sensitivity.csv"),
          row.names = FALSE, na = "")

###############################################################################
# 10. 결측 프로파일 (구조적 vs 산발적) + 성별 차등결측 검정
###############################################################################

allrows <- rbind(t1, t2, t3, t4)
miss_vars <- allrows %>% filter(type %in% c("continuous", "binary", "categorical")) %>%
  dplyr::select(table, section, variable, label, denominator, missing_type,
                n_overall, n_male, n_female, n_missing, pct_missing) %>%
  distinct()

diff_miss <- lapply(seq_len(nrow(miss_vars)), function(i) {
  v <- miss_vars$variable[i]
  base <- if (grepl("restricted denominator", miss_vars$missing_type[i])) {
    if (v == "Retro_success") co[co$Retro_attempt_re %in% c(1, 2), ]
    else if (v %in% c("n_stent_f", "stent_dia_mean", "stent_dia_min",
                      "stent_dia_max", "stent_len_total")) co[co$TV_STENT %in% 1, ]
    else if (v == "dapt_potent") co[co$dc_antiplt_2 %in% c("ac", "ap", "at"), ]
    else co
  } else co
  if (!v %in% names(base)) return(data.frame(variable = v, p_diff_missing = NA_real_))
  m <- is.na(base[[v]])
  if (length(unique(m)) < 2) return(data.frame(variable = v, p_diff_missing = NA_real_))
  p <- tryCatch(suppressWarnings(stats::chisq.test(table(base$sex_f, m))$p.value),
                error = function(e) NA_real_)
  data.frame(variable = v, p_diff_missing = p)
}) %>% bind_rows()

miss_out <- dplyr::left_join(miss_vars, diff_miss, by = "variable") %>%
  mutate(p_diff_missing_fmt = fmt_p(p_diff_missing))
write.csv(miss_out, file.path(DIR_TAB, "tableS6_missingness_profile_by_sex.csv"),
          row.names = FALSE, na = "")

###############################################################################
# 11. 코호트 진단 로그 (STROBE flow 용 수치 포함)
###############################################################################

full <- build_cohort(filter_cohort = FALSE, verbose = FALSE)
diag_rows <- data.frame(
  item = c("Registry total",
           "  CTO PCI performed (flag_cto_pci)",
           "  CTO PCI & index year 2014-2023 (ANALYSIS COHORT)",
           "    Male", "    Female",
           "  Excluded: CTO PCI but index year <2014 or >2023",
           "Contributing centers (anonymized)",
           "Sensitivity cohort: index year 2014-2022",
           "Zero follow-up & no event (death_fudur=0 & Death=0)",
           "  Male", "  Female",
           "Day-0 events: DMS (primary)", "Day-0 events: all-cause death",
           "Day-0 events: MACCE_pub", "Day-0 events: NACE_pub",
           "Day-0 events: CV death",
           "flag_logic_invalid_last_FU_text in cohort",
           "flag_logic_fudur_exceeds_total in cohort"),
  value = c(nrow(full), sum(full$flag_cto_pci %in% TRUE), nrow(co),
            sum(co$sex_f == "Male"), sum(co$sex_f == "Female"),
            sum(full$flag_cto_pci %in% TRUE) - nrow(co),
            length(unique(co$center_anon)), nrow(co_win),
            sum(co$Death == 0 & co$death_fudur == 0),
            sum(co$Death == 0 & co$death_fudur == 0 & co$sex_f == "Male"),
            sum(co$Death == 0 & co$death_fudur == 0 & co$sex_f == "Female"),
            sum(co$DMS == 1 & co$DMS_fudur == 0),
            sum(co$Death == 1 & co$death_fudur == 0),
            sum(co$MACCE_pub == 1 & co$MACCE_pub_fudur == 0),
            sum(co$NACE_pub == 1 & co$NACE_pub_fudur == 0),
            sum(co$CVdeath == 1 & co$death_fudur == 0),
            sum(co$flag_logic_invalid_last_FU_text %in% TRUE),
            sum(co$flag_logic_fudur_exceeds_total %in% TRUE)),
  stringsAsFactors = FALSE)
p_fu0 <- suppressWarnings(stats::chisq.test(
  table(co$sex_f, co$Death == 0 & co$death_fudur == 0))$p.value)
diag_rows <- rbind(diag_rows,
  data.frame(item = "Chi-square p, zero follow-up by sex (Limitations only; no sensitivity analysis per investigator)",
             value = round(p_fu0, 5)))
write.csv(diag_rows, file.path(DIR_TAB, "tableS7_analysis_diagnostics.csv"), row.names = FALSE, na = "")

###############################################################################
# 12. docx 출력
###############################################################################

try({
  if (requireNamespace("flextable", quietly = TRUE) && requireNamespace("officer", quietly = TRUE)) {
    mk <- function(df, title) {
      d <- df %>% dplyr::transmute(
        Variable = label,
        `N available` = ifelse(is.na(n_overall), "", as.character(n_overall)),
        `Overall` = ifelse(is.na(overall), "", overall),
        `Male` = ifelse(is.na(male), "", male),
        `Female` = ifelse(is.na(female), "", female),
        `p value` = ifelse(is.na(p_fmt), "", p_fmt),
        `SMD` = ifelse(is.na(smd), "", sprintf("%.3f", smd)))
      flextable::autofit(flextable::set_caption(flextable::flextable(d), title))
    }
    doc <- officer::read_docx()
    doc <- flextable::body_add_flextable(doc, mk(t1, "Table 1. Baseline characteristics by sex (N=6,777)"))
    doc <- officer::body_add_break(doc)
    doc <- flextable::body_add_flextable(doc, mk(t2, "Table 2. Coronary lesion characteristics by sex"))
    doc <- officer::body_add_break(doc)
    doc <- flextable::body_add_flextable(doc, mk(t3, "Table 3. Procedural characteristics and discharge medication by sex"))
    doc <- officer::body_add_break(doc)
    doc <- flextable::body_add_flextable(doc, mk(t4, "Table 4. In-hospital outcomes by sex"))
    print(doc, target = file.path(DIR_TAB, "tables1to4_by_sex.docx"))
  }
}, silent = TRUE)

###############################################################################
# 13. sessionInfo
###############################################################################

capture.output({
  cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Script: scripts/02_descriptive_by_sex.R\n")
  cat("Cohort:", COHORT_LABEL, "\n")
  cat("set.seed(42)\n\n")
  print(sessionInfo())
}, file = file.path(DIR_REPX, "sessionInfo_02_descriptive.txt"))

message("== DONE (02_descriptive_by_sex.R) ==")
