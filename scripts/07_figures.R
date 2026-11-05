###############################################################################
# K-CTO Registry 2.0 — Sex Difference Analysis
# Phase 4a (illustrator): Figure 1-4 (main manuscript figures)
#
# 실행:  Rscript scripts/07_figures.R
# 입력:  output/tables/*.csv  (Phase 2 statistician 산출물 — 하드코딩 금지,
#        모든 수치는 아래에서 CSV로부터 직접 읽는다)
# 출력:  output/figures/Figure{1..3}_*.tiff  (1200 dpi; PNG 중단 — 연구자 확정 2026-08-12)
#        output/figures/captions.md
#
# Figure list
#   Figure 1  STROBE flow diagram (registry -> analysis cohort -> analytic n)
#   Figure 2  Age: 5-year incidence by age band + hazard ratio within age band
#   Figure 3  Procedural intensity: score distribution + bilateral injection
#   (구 Figure 4 및 구 2A/3B 는 Central Illustration 으로 이동 — 아래 주석 참조)
#
# R 전용, ggplot2 + patchwork. Colorblind-safe (Okabe-Ito) + shape/linetype
# encoding for sex. All figures exported at line-art dpi (1200) satisfying
# both JACC:CI and Circulation:CVI thresholds (journal_config.json).
#
# 금지 문구 준수 (CLAUDE.md §1): "동일한 병변 복잡도", "장기 예후 동등",
# "치료 과정 격차 5개", "unexplained disparity", "pre-specified", 인과 언어.
# PHI: id/patient_key/id_code/Center_name 미사용 (output/tables/*.csv만 사용,
# 이미 익명화된 center_anon 집계치만 포함됨).
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(stringr)
  library(forcats)
})

PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "_theme.R"))

rt <- function(f) read.csv(file.path(TABLES_DIR, f), stringsAsFactors = FALSE, check.names = FALSE)

###############################################################################
## FIGURE 1 — STROBE flow diagram
##   Source: output/tables/tableS7_analysis_diagnostics.csv (registry -> 6,777
##   cohort counts, zero-follow-up); tableS15b_person_time_contributing_n.csv
##   (model-specific analytic samples: complete-case 6,259; person-time
##   contributing for DMS 5,597).
###############################################################################

s7  <- rt("tableS7_analysis_diagnostics.csv")
s15b <- rt("tableS15b_person_time_contributing_n.csv")

gv <- function(item_pattern) {
  hit <- s7$value[grepl(item_pattern, s7$item)]
  if (length(hit) == 0) stop("tableS7 item not found: ", item_pattern)
  hit[1]
}

n_registry   <- gv("^Registry total$")
n_ctopci     <- gv("CTO PCI performed")
n_cohort     <- gv("ANALYSIS COHORT")
n_excluded   <- gv("Excluded: CTO PCI but index year")
n_zerofu     <- gv("Zero follow-up")
n_centers    <- gv("Contributing centers")

# tableS7 item column: cohort Male/Female rows are indented 4 spaces
# ("    Male"), zero-follow-up Male/Female rows are indented 2 spaces
# ("  Male") — indentation level, not row order, distinguishes the two blocks.
n_cohort_male    <- s7$value[s7$item == "    Male"]
n_cohort_female  <- s7$value[s7$item == "    Female"]
n_zerofu_male    <- s7$value[s7$item == "  Male"]
n_zerofu_female  <- s7$value[s7$item == "  Female"]
p_zerofu <- as.numeric(gv("Chi-square p, zero follow-up"))
p_zerofu_fmt <- ifelse(p_zerofu < 0.001, "P<0.001", sprintf("P=%.3f", p_zerofu))

n_male_pct   <- round(100 * n_cohort_male   / n_cohort, 1)
n_female_pct <- round(100 * n_cohort_female / n_cohort, 1)
zerofu_male_pct   <- round(100 * n_zerofu_male   / n_cohort_male, 1)
zerofu_female_pct <- round(100 * n_zerofu_female / n_cohort_female, 1)

# ①(2026-08-27) 첫 제외 단계 인원. 그림에 이미 인쇄된 두 값의 차이이며 새 수치가 아니다.
n_nopci <- n_registry - n_ctopci

n_completecase <- s15b$n_model_reported[s15b$outcome_key == "MACCE_pub"]        # v2: primary = MACCE
n_persontime   <- s15b$n_contributing_person_time[s15b$outcome_key == "MACCE_pub"]  # 5,650

# ---- SUBMISSION-VERSION NOTE (illustrator, 2026-08-10, QA BLOCK1) ----------
# In-figure title ("Figure 1. STROBE flow diagram") and the bottom
# source/provenance caption (internal CSV paths, CLAUDE.md section refs,
# raw variable names such as death_fudur/flag_cto_pci) are intentionally NOT
# rendered into the submission image. The journal typesets the figure number
# and title; the caption text below is submitted as manuscript body text
# (see interim_report/phase4_manuscript/01_manuscript.md, FIGURE LEGENDS,
# and interim_report/phase4_manuscript/figure_inventory.md). Provenance is
# preserved here in this script comment and in figure_inventory.md, not lost.
#   Source: output/tables/tableS7_analysis_diagnostics.csv;
#   tableS15b_person_time_contributing_n.csv.
#   Cohort definition: flag_cto_pci == TRUE & index year in [2014, 2023]
#   (CLAUDE.md section 1, "확정 사항").
#   "No follow-up beyond index procedure" box = death_fudur == 0 and no
#   in-hospital event recorded (raw variable name intentionally not printed
#   inside the box in the submission image; see figure_inventory.md).

# ---- Build flow diagram with geom_rect + geom_segment (manual coordinates)
box <- function(xc, yc, w, h, label, fill = "white", border = "grey20", size = 2.5, fontface = "plain") {
  list(
    geom_rect(aes(xmin = xc - w/2, xmax = xc + w/2, ymin = yc - h/2, ymax = yc + h/2),
              fill = fill, color = border, linewidth = 0.45),
    annotate("text", x = xc, y = yc, label = label, size = size, lineheight = 0.92, fontface = fontface)
  )
}
arrow_down <- function(x, y1, y2) {
  geom_segment(aes(x = x, xend = x, y = y1, yend = y2),
               arrow = arrow(length = unit(2, "mm"), type = "closed"), linewidth = 0.45, color = "grey20")
}
arrow_right <- function(x1, x2, y) {
  geom_segment(aes(x = x1, xend = x2, y = y, yend = y),
               arrow = arrow(length = unit(2, "mm"), type = "closed"), linewidth = 0.45, color = "grey30")
}

# NOTE(Lead, 2026-08-10): 레이아웃 재구성.
#  (1) 제외 화살표의 화살촉이 상자 내부(x=8.85 < 상자 좌변 7.7)에 그려져 보이지 않던 문제 수정.
#  (2) 모델표본 주석 상자로 가는 점선이 y=6.55에서 끝나 상자(y<=5.25)에 닿지 않고 공중에 뜨던
#      문제를 elbow 연결자로 수정.
#  (3) STROBE 원칙에 맞춰 **주 흐름(실선 화살표)은 분석 코호트에서 끝낸다.** 추적 0일 인원과
#      모델별 표본은 제외 단계가 아니라 같은 코호트의 부가 정보이므로 **점선 주석**으로 분리.
dash_elbow <- function(x1, y1, x2, y2, xmid) {
  list(geom_segment(aes(x = x1, xend = xmid, y = y1, yend = y1),
                    linewidth = 0.4, color = "grey45", linetype = "42"),
       geom_segment(aes(x = xmid, xend = xmid, y = y1, yend = y2),
                    linewidth = 0.4, color = "grey45", linetype = "42"),
       geom_segment(aes(x = xmid, xend = x2, y = y2, yend = y2),
                    linewidth = 0.4, color = "grey45", linetype = "42"))
}

p1 <- ggplot() +
  coord_cartesian(xlim = c(0, 12.3), ylim = c(3.3, 10.7), clip = "off") +
  # ---- 주 흐름 (실선) : 등록 -> CTO PCI -> 분석 코호트에서 종료 ----
  box(3.4, 10.05, 6.6, 1.0,
      sprintf("K-CTO Registry 2.0, all patients\nN = %s", format(n_registry, big.mark = ","))) +
  arrow_down(3.4, 9.55, 8.87) +
  # ① 첫 단계 제외 분기 (STROBE 13b: 각 단계의 제외 인원·사유)
  arrow_right(3.4, 7.61, 9.24) +
  box(9.9, 9.24, 4.5, 0.85,
      sprintf("Excluded (n = %s): CTO PCI\nnot performed", format(n_nopci, big.mark = ",")),
      size = 2.3) +
  box(3.4, 8.35, 6.6, 1.0,
      sprintf("CTO PCI performed\nN = %s", format(n_ctopci, big.mark = ","))) +
  arrow_down(3.4, 7.85, 7.20) +
  box(3.4, 6.55, 6.6, 1.25,
      sprintf("Analysis cohort: CTO PCI & index year 2014-2023\nN = %s  (Men %s [%.1f%%] / Women %s [%.1f%%])\n%d contributing centers",
              format(n_cohort, big.mark = ","), format(n_cohort_male, big.mark = ","), n_male_pct,
              format(n_cohort_female, big.mark = ","), n_female_pct, as.integer(n_centers)),
      fill = "grey90", size = 2.5) +   # ② 음영은 분석 코호트에만
  # ---- 제외 분기 (실선 화살표) : 세로 화살표에서 갈라져 상자 좌변에서 정확히 종료 ----
  arrow_right(3.4, 7.61, 7.52) +
  box(9.9, 7.52, 4.5, 1.15,
      sprintf("Excluded (n = %s): CTO PCI with\nindex procedure in 2003-2013\n(outside the 2014-2023 window)",
              format(n_excluded, big.mark = ",")), size = 2.3) +
  # ---- 부가 정보 (점선) : 제외가 아니라 같은 코호트에 대한 주석 ----
  dash_elbow(3.4, 5.925, 3.4, 5.18, 3.4) +
  box(3.4, 4.42, 6.6, 1.45,   # ④ 우측 상자와 높이 통일
      sprintf(paste0(
        "Retained in the analysis cohort, but with no follow-up beyond\nthe index procedure (and no in-hospital event): n = %s (%.1f%%)\n",
        "Men %s (%.1f%% of men) / Women %s (%.1f%% of women), %s"),
        format(n_zerofu, big.mark = ","), 100*n_zerofu/n_cohort,
        format(n_zerofu_male, big.mark = ","), zerofu_male_pct,
        format(n_zerofu_female, big.mark = ","), zerofu_female_pct,
        p_zerofu_fmt),
      border = "grey45", size = 2.3) +
  dash_elbow(6.70, 6.55, 7.65, 4.42, 7.20) +
  box(9.9, 4.42, 4.5, 1.45,   # ④ 좌측 상자와 높이 통일
      sprintf(paste0(
        # ③ 6,259 / 5,650 이 별개 코호트로 읽히지 않도록 "subsets" 명시
        "Model-specific analytic samples (primary outcome),\nall subsets of the same %s-patient cohort\n",
        "Complete-case (Model 4 covariates): n = %s (%.1f%%)\n",
        "Contributing follow-up time: n = %s"),
        format(n_cohort, big.mark = ","),
        format(n_completecase, big.mark = ","), 100*n_completecase/n_cohort,
        format(n_persontime, big.mark = ",")),
      border = "grey45", size = 2.3) +
  # 점선 = 부가 정보(제외 아님) 임을 명시
  theme_flow(base_size = 8)

save_fig(p1, file.path(FIGURES_DIR, "Figure1_STROBE_flow"), width_mm = 183, height_mm = 108)
cat("Figure 1 saved.\n")

###############################################################################
## FIGURE 2 — Staged adjustment (age attenuation) forest plot +
##             age-stratified 5-year cumulative incidence (composition effect)
##   Source: output/tables/table7_adjusted_HR_all_outcomes.csv (staged HRs,
##   4 outcomes x 3 stages); table6b_age_attenuation_DMS.csv (% excess risk
##   explained by age, DMS); table6f_age_stratified_DMS.csv (5-yr cumulative
##   incidence and within-age-stratum HR, male vs female).
###############################################################################

t7  <- rt("table7_adjusted_HR_all_outcomes.csv")
t6b <- rt("table6b_age_attenuation_DMS.csv")
t6f <- rt("table6f_age_stratified_DMS.csv")

outcome_lab <- c(
  "DMS (death, MI, or stroke) - PRIMARY" = "DMS (death, MI, or\nstroke) — Primary",
  "All-cause death"                      = "All-cause death",
  "MACCE (published definition)"         = "MACCE",
  "Any coronary revascularization"       = "Any revascularization"
)
model_lab <- c("Unadjusted" = "Unadjusted", "Age-adjusted" = "Age-adjusted", "Fully adjusted" = "Fully adjusted")

f2a_dat <- t7 %>%
  filter(outcome %in% names(outcome_lab), model %in% names(model_lab)) %>%
  mutate(
    outcome_f = factor(outcome_lab[outcome], levels = rev(unname(outcome_lab))),
    model_f   = factor(model_lab[model], levels = unname(model_lab)),
    lab       = sprintf("%s", HR_fmt),
    label_x   = hi * 1.14   # multiplicative offset (log scale) so labels never overlap the CI bar
  )

# ---- QA BLOCK3 fix (2026-08-10): the Panel A subtitle previously reused the
# already-rounded string values in table6b_age_attenuation_DMS.csv (which
# stores "1.005" -> round()-to-2dp gives 1.00 because of floating-point
# representation) instead of the unrounded numeric HR column in table7 (which
# gives 1.01, matching the forest plot drawn from the same table7 in this
# panel, the main text, and CLAUDE.md). To prevent this recurring, the staged
# HR values quoted in the subtitle are now read directly from the same t7
# rows used to draw the forest plot above (single source of truth); only the
# "% of excess risk explained by age" figure (not present in table7) is still
# taken from table6b.
dms_stage <- f2a_dat %>% filter(outcome == "DMS (death, MI, or stroke) - PRIMARY")
hr_unadj_num <- dms_stage$HR[dms_stage$model == "Unadjusted"]
hr_age_num   <- dms_stage$HR[dms_stage$model == "Age-adjusted"]
hr_full_num  <- dms_stage$HR[dms_stage$model == "Fully adjusted"]

pA <- ggplot(f2a_dat, aes(x = HR, y = model_f, color = model_f, shape = model_f)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.35) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.22, linewidth = 0.5, orientation = "y") +
  geom_point(size = 2.0) +
  geom_text(aes(x = label_x, label = lab), hjust = 0, size = 2.05, color = "grey15", show.legend = FALSE) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2), limits = c(0.55, 4.6)) +
  scale_color_manual(values = setNames(stage_palette[c(1,3,5)], unname(model_lab)), name = "Adjustment stage") +
  scale_shape_manual(values = c(15, 16, 17), name = "Adjustment stage") +
  facet_grid(rows = vars(outcome_f), switch = "y", scales = "free_y", space = "free_y") +
  labs(x = "Hazard ratio, female vs male (log scale)", y = NULL,
       title = "A. Staged covariate adjustment attenuates the female excess hazard",) +
  theme_kcto(base_size = 8) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold", size = rel(0.85)),
        legend.position = "top", panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = rel(0.85)))

# ---- Panel B: age-stratified 5-year cumulative incidence (composition effect)
# v2(2026-08-14): DMS(table6f) → MACCE(table26) 로 교체 (연구자 결정)
t26 <- rt("table26_age_stratified_MACCE.csv")
t6f2 <- t26 %>%
  mutate(
    age_stratum_f = factor(age_stratum, levels = c("<65 years", "65-74 years", ">=75 years"),
                            labels = c("<65 y", "65-74 y", "≥75 y")),
    sex = sex_factor(sex)
  )

t6f2 <- t6f2 %>% mutate(label_y = ci_hi + 0.9)
hr_lab <- t6f2 %>% distinct(age_stratum_f, within_stratum_HR, within_stratum_p) %>%
  mutate(lab = sprintf("HR %s, P=%s", within_stratum_HR, within_stratum_p))

pB <- ggplot(t6f2, aes(x = age_stratum_f, y = cum_inc_5yr_pct, fill = sex)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65, color = "grey15", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), position = position_dodge(width = 0.72),
                width = 0.18, linewidth = 0.4, color = "grey15") +
  geom_text(aes(y = label_y, label = sprintf("%.1f", cum_inc_5yr_pct)), position = position_dodge(width = 0.72),
            vjust = 0, size = 2.15) +
  scale_fill_manual(values = sex_fills, name = NULL) +   # 범례 제목 제거 (연구자 지시 2026-08-12)
  scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Age stratum", y = "5-year cumulative incidence of MACCE, %",) +
  # 그림 내부 제목 삭제(2026-08-28, 검토 메모 id 50) — 제목은 FIGURE LEGENDS 가 담당.
  # (구 제목 "Five-year cumulative incidence by age band and sex"; QA BLOCK2 문구 이력은 git/legacy 참조)
  theme_kcto(base_size = 8) +
  theme(legend.position = "top")

# ---- NOTE(Lead, 2026-08-10): Figure 2 는 pB4(연령대내 HR)가 정의된 뒤에
# 조립한다. 스크립트 맨 아래 "FIGURE 2 (재구성)" 절 참조.
# 구 Panel A(staged adjustment forest)는 Central Illustration 으로 이동했다.

###############################################################################
## FIGURE 3 — Process-of-care gap: intensity-score distribution + gap that
##            present at every complexity stratum. NOTE(2026-08-10): 상호작용 검정
##            결과(A18) 반영 — "widens/gradient" 주장 철회, 순서만 기술.
##   Source: output/tables/table15b_procedural_intensity_distribution.csv
##   (panel A, score distribution); table15c_intensity_score_robustness.csv
##   (panel B, complexity-stratified score gap, 4- and 3-component scores);
##   table14_within_center_management.csv, table18_escalation_conditional_
##   on_indication.csv (panel C, bilateral injection, the single most robust
##   individual indicator).
###############################################################################

t15b <- rt("table15b_procedural_intensity_distribution.csv")
t15c <- rt("table15c_intensity_score_robustness.csv")
t14  <- rt("table14_within_center_management.csv")
t18  <- rt("table18_escalation_conditional_on_indication.csv")

## ---- Panel A: distribution of the 4-component intensity-of-care score ----
distA <- t15b %>%
  filter(score_definition == "Complete case (all 4 components observed)") %>%
  select(sex, pct_0, pct_1, pct_2, pct_3, pct_4, n, mean_sd) %>%
  pivot_longer(cols = starts_with("pct_"), names_to = "score", values_to = "pct") %>%
  mutate(score = factor(sub("pct_", "", score), levels = as.character(0:4)),
         sex = sex_factor(sex))

score_pal <- colorRampPalette(c("grey90", unname(okabe_ito["blue"])))(5)

n_distA_total <- sum(distA$n[!duplicated(distA$sex)])

# NOTE(Lead, 2026-08-10): coord_flip 기본 적층은 마지막 수준(4)이 왼쪽에 오므로
# 막대의 시각 순서(4->0)와 범례 순서(0->4)가 반대였다. reverse=TRUE 로 적층을
# 0->4 오름차순으로 바꿔 범례와 일치시킨다.
pA3 <- ggplot(distA, aes(x = sex, y = pct, fill = score)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.3,
           position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = ifelse(pct >= 4, sprintf("%.0f", pct), "")),
            position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.1, color = "grey15") +
  scale_fill_manual(values = setNames(score_pal, as.character(0:4)),
                     name = "Intensity-of-care\nscore (0-4)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_flip() +
  labs(x = NULL, y = "Patients, %",) +
  # 그림 내부 제목 삭제(2026-08-28, 검토 메모 id 51) — 제목은 FIGURE LEGENDS 가 담당.
  theme_kcto(base_size = 8) +
  theme(legend.position = "right", axis.text.y = element_text(size = rel(1.05)))

## ---- Panel B: complexity-stratified score gap (4- and 3-component) -------
## Ordered by increasing baseline lesion complexity; the ordering is descriptive only
## directly top-to-bottom: Overall -> J-CTO <2 (low) -> J-CTO >=2 (high) ->
## occlusion length >=20 mm (high, independent complexity measure).
strat_lab <- c(
  "Full analysis cohort"                                        = "Overall",
  "Baseline complexity stratum: J-CTO score <2 (contrast)"       = "J-CTO score <2",
  "Baseline complexity stratum: J-CTO score >=2"                 = "J-CTO score ≥2",
  "Baseline complexity stratum: occlusion length >=20 mm"        = "Occlusion length ≥20 mm"
)
strat_order <- c("Overall", "J-CTO score <2", "J-CTO score ≥2", "Occlusion length ≥20 mm")

get_score_rows <- function(outcome_str, n_target) {
  t15c %>%
    filter(model == "Linear mixed model, centre random intercept",
           outcome == outcome_str,
           covariate_set %in% c("Model 4 covariates minus centre volume (absorbed by centre)", NA) | is.na(covariate_set),
           stratum %in% names(strat_lab),
           n %in% n_target) %>%
    distinct(stratum, n, n_male, n_female, estimate, lo, hi, p_fmt)
}
# Overall rows for A16-0/A16-1 lack a covariate_set value (NA) — filtered via n instead
b4 <- t15c %>% filter(analysis == "A16-0", model == "Linear mixed model, centre random intercept", n == 5097) %>%
  distinct(stratum, n, n_male, n_female, estimate, lo, hi, p_fmt) %>% mutate(stratum = "Full analysis cohort")
b4strat <- t15c %>% filter(analysis == "A16-4", model == "Linear mixed model, centre random intercept",
                            outcome == "Intensity score 0-4 (4 components)",
                            covariate_set == "Model 4 covariates minus centre volume (absorbed by centre)",
                            stratum %in% names(strat_lab)) %>%
  distinct(stratum, n, n_male, n_female, estimate, lo, hi, p_fmt)
score4 <- bind_rows(b4, b4strat) %>% mutate(score_type = "4-component (bilateral, retrograde,\nIVUS, potent P2Y12)")

b3 <- t15c %>% filter(analysis == "A16-1", model == "Linear mixed model, centre random intercept", n == 5879,
                       outcome == "Intensity score 0-3, excluding potent P2Y12 at discharge") %>%
  distinct(stratum, n, n_male, n_female, estimate, lo, hi, p_fmt) %>% mutate(stratum = "Full analysis cohort")
b3strat <- t15c %>% filter(analysis == "A16-4", model == "Linear mixed model, centre random intercept",
                            outcome == "Intensity score 0-3, excluding potent P2Y12 at discharge",
                            covariate_set == "Model 4 covariates minus centre volume (absorbed by centre)",
                            stratum %in% names(strat_lab)) %>%
  distinct(stratum, n, n_male, n_female, estimate, lo, hi, p_fmt)
score3 <- bind_rows(b3, b3strat) %>% mutate(score_type = "3-component (primary: bilateral,\nretrograde, IVUS only)")

f3b_dat <- bind_rows(score4, score3) %>%
  mutate(stratum_f = factor(strat_lab[stratum], levels = rev(strat_order)),
         score_type = factor(score_type, levels = c("3-component (primary: bilateral,\nretrograde, IVUS only)",
                                                      "4-component (bilateral, retrograde,\nIVUS, potent P2Y12)")),
         lab = sprintf("%.3f (%.3f to %.3f)", estimate, lo, hi))

pB3 <- ggplot(f3b_dat, aes(x = estimate, y = stratum_f, color = score_type, shape = score_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.35) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.18, linewidth = 0.55,
                position = position_dodge(width = 0.5), orientation = "y") +
  geom_point(size = 2.1, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c(unname(okabe_ito["vermillion"]), unname(okabe_ito["blue"])), name = NULL) +
  scale_shape_manual(values = c(17, 16), name = NULL) +
  scale_x_continuous(limits = c(-0.32, 0.06)) +
  labs(x = "Adjusted mean score difference,\nfemale vs male (score points)", y = NULL,
       # NOTE(Lead, 2026-08-10): "at every level"은 전칭 명제로, J-CTO <2 층의 경계값
# (P=0.049)에 의존한다. 원고에서 동일 표현을 철거했으므로 figure도 기술적 표현으로 통일.
       title = str_wrap("B. Score gap by baseline lesion complexity", width = 42),) +
  theme_kcto(base_size = 8) +
  theme(legend.position = "top", legend.text = element_text(size = rel(0.72)),
        plot.title = element_text(size = rel(1.0)), plot.subtitle = element_text(size = rel(0.85)))

## ---- Panel C: bilateral injection — the single most robust indicator -----
bilat_overall <- t14 %>% filter(outcome == "Bilateral injection",
                                 model == "(2) Within-centre: random intercept for centre (full Laplace)") %>%
  transmute(stratum_f = "Overall", n, estimate, lo, hi, p_fmt)
bilat_strat <- t18 %>% filter(outcome == "Bilateral injection", grepl("centre random intercept", model),
                               grepl("J-CTO score >=2|occlusion length >=20 mm", model)) %>%
  mutate(stratum_f = case_when(
    grepl("J-CTO score >=2", model) ~ "J-CTO score ≥2",
    grepl("occlusion length >=20 mm", model) ~ "Occlusion length\n≥20 mm"
  )) %>%
  transmute(stratum_f, n, estimate, lo, hi, p_fmt)
f3c_dat <- bind_rows(bilat_overall, bilat_strat) %>%
  mutate(stratum_f = factor(stratum_f, levels = rev(c("Overall", "J-CTO score ≥2", "Occlusion length\n≥20 mm"))),
         lab = sprintf("%.2f (%.2f-%.2f)\nn=%s", estimate, lo, hi, format(n, big.mark=",")),
         label_x = 1.12)

pC3 <- ggplot(f3c_dat, aes(x = estimate, y = stratum_f)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.35) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.16, linewidth = 0.55, color = unname(okabe_ito["blue"])) +
  geom_point(size = 2.3, shape = 18, color = unname(okabe_ito["blue"])) +
  geom_text(aes(x = label_x, label = lab), hjust = 0, size = 2.0, color = "grey15", lineheight = 0.85) +
  # NOTE(Lead): 패널 폭 대비 눈금이 많아 "0.751.00" 처럼 라벨이 겹쳤다. 눈금 축소.
  scale_x_log10(limits = c(0.45, 2.45), breaks = c(0.5, 0.75, 1.0, 1.5, 2.0),
                labels = c("0.50", "0.75", "1.00", "1.50", "2.00")) +
  labs(x = "Adjusted odds ratio, female vs male (log scale)", y = NULL,
       # ①(2026-08-27) 패널 제목은 기술만 — "most robust" 단정 제거
       title = "B. Bilateral contrast injection, overall and by lesion complexity",) +
  theme_kcto(base_size = 8) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(lineheight = 0.9),
        plot.subtitle = element_text(size = rel(0.85)))

# ---- NOTE(Lead, 2026-08-10): 구 Panel B(복잡도 층별 강도 점수 평균차)는
# Central Illustration 으로 이동했다.
## 연구자 지시(2026-08-12): ① A-B 사이 1줄 간격(pC3 상단 여백) ② 183→140mm 축소
## 연구자 지시(2026-08-28, 검토 메모 id 30): 구 Panel B(bilateral forest, pC3) 삭제 —
## CI 우하단 층별 forest 와 시각 구조가 중복. bilateral 층별 수치는 본문·Table 5·S9 에
## 인쇄되어 손실 없음. pC3 는 출처 기록으로 정의만 남긴다(저장 안 함).
p3 <- pA3 + theme(plot.margin = margin(3, 8, 3, 3))
save_fig(p3, file.path(FIGURES_DIR, "Figure3_process_of_care_gap"),
         width_mm = 140, height_mm = 60)
cat("Figure 3 saved.\n")

###############################################################################
## FIGURE 4 — Kaplan-Meier curves for DMS by sex, with number-at-risk table,
##            and an age-stratified panel showing the gap is a composition
##            effect (women are older, not at higher risk within age band).
##   Source: output/tables/table5e_km_curve_coordinates.csv (KM step curve);
##   table5b_km_summary_by_sex.csv (N, log-rank P, crude HR); table6f_age_
##   stratified_DMS.csv (age-stratified 5-yr cumulative incidence).
###############################################################################

t5e <- rt("table5e_km_curve_coordinates.csv")
t5b <- rt("table5b_km_summary_by_sex.csv")

km_dms <- t5e %>% filter(outcome_var == "DMS") %>%
  mutate(sex = sex_factor(sex), time_yr = time_days / 365.25)
b_dms  <- t5b %>% filter(outcome_var == "DMS")
b_dms_p <- b_dms$logrank_p_fmt[1]

t5a <- rt("table5a_km_timepoints_by_sex.csv")
end5 <- t5a %>%
  filter(outcome_var == "DMS", time_days >= 1825, time_days <= 1827) %>%
  mutate(sex = sex_factor(sex),
         lab = sprintf("%.1f%% (%.1f-%.1f)", cum_incidence_pct, ci_lower_pct, ci_upper_pct))
lr_lab <- sprintf("Log-rank P%s", ifelse(b_dms_p == "<0.001", "<0.001", paste0("=", b_dms_p)))


# ---- QA "그 밖에" fix (2026-08-10, NC-06/M7): table5b's crude HR (1.48,
# 1.25-1.74) is fitted on the full n=6,777 cohort, while Figure 2 / the main
# text / CLAUDE.md report the unadjusted DMS HR from the same complete-case
# sample used throughout the staged-adjustment chain (n=6,259; 1.48,
# 1.23-1.76, table7_adjusted_HR_all_outcomes.csv). Two different unadjusted
# CIs for the same-looking "1.48" printed on two figures is a within-package
# inconsistency; unify on the table7 value (same source as Figure 2) so a
# single unadjusted DMS estimate is used everywhere in the figure set.
hr_unadj_fig4 <- t7$HR_fmt[t7$outcome == "DMS (death, MI, or stroke) - PRIMARY" & t7$model == "Unadjusted"]

## ---- Panel A: KM curves (cumulative incidence, 1 - survival) -------------
## NOTE (연구자 지시, 2026-08-09): KM 곡선에 신뢰구간 밴드를 그리지 않는다.
## 구간 추정의 불확실성은 number-at-risk 표와 본문 수치(1/3/5년 누적발생률 + 95% CI,
## table5a)로 전달한다. 리뷰어가 CI 표시를 요구하면 아래 geom_ribbon 한 줄을 되살리면 된다.
## (ci_lower_pct / ci_upper_pct 컬럼은 table5e_km_curve_coordinates.csv 에 그대로 있음)
pA4 <- ggplot(km_dms, aes(x = time_yr, y = cum_incidence_pct, color = sex, fill = sex, linetype = sex)) +
  # geom_ribbon(aes(ymin = ci_lower_pct, ymax = ci_upper_pct), alpha = 0.14, color = NA) +
  geom_step(linewidth = 0.62) +
  # 연구자 요청(2026-08-10): 곡선 끝에 5년 누적발생률(95% CI)과 log-rank P 표시
  geom_text(data = end5, aes(x = 5.08, y = cum_incidence_pct, label = lab, color = sex),
            inherit.aes = FALSE, hjust = 0, size = 2.2, show.legend = FALSE) +
  annotate("text", x = 3.0, y = 2.2, hjust = 0, size = 2.4, color = "grey25",
           label = lr_lab) +
  scale_color_manual(values = sex_colors, name = "Sex") +
  scale_fill_manual(values = sex_fills, name = "Sex") +
  scale_linetype_manual(values = sex_linetypes, name = "Sex") +
  scale_x_continuous(breaks = 0:5, limits = c(-0.15, 6.6), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(limits = c(0, 26), expand = expansion(mult = c(0, 0.02))) +
  # 연구자 지시(2026-08-09): x축 눈금·제목은 KM 곡선 바로 아래에 두고,
  # number-at-risk 표는 그 아래로 내린다.
  labs(x = "Years since index CTO PCI", y = "Cumulative incidence of DMS, %",
       title = "A. Death, myocardial infarction, or stroke (DMS), by sex",) +
  theme_kcto(base_size = 8) +
  # NOTE(Lead, 2026-08-10, 연구자 요청): KM 곡선의 세로 보조선 제거,
  # y축선을 x축과 동일한 굵기·색으로, 축 타이틀 여백을 x축 기준에 맞춰 축소.
  theme(legend.position = c(0.14, 0.86),
        legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.2),
        panel.grid.major.x = element_blank(),
        axis.line.y  = element_line(color = "grey20", linewidth = 0.35),
        axis.ticks.y = element_line(color = "grey20", linewidth = 0.3),
        axis.title.y = element_text(margin = margin(r = 0)),
        axis.title.x = element_text(margin = margin(t = 1.5, b = 0)),
        plot.margin = margin(4, 8, 0, 4))

## ---- Number-at-risk table, read off the same KM step coordinates ---------
at_risk_at <- function(df, days_target) {
  purrr_map <- lapply(days_target, function(dtar) {
    df2 <- df[df$time_days <= max(dtar, min(df$time_days)), ]
    if (nrow(df2) == 0) df2 <- df[1, ]
    df2[which.max(df2$time_days), ]
  })
  do.call(rbind, purrr_map)
}
landmarks <- c(0, 365.25, 730.5, 1095.75, 1461, 1826.25)
risk_tab <- bind_rows(lapply(sex_display, function(s) {
  sub <- km_dms %>% filter(sex == s) %>% arrange(time_days)
  rt2 <- at_risk_at(sub, landmarks)
  rt2$landmark_yr <- 0:5
  rt2
})) %>% select(sex, landmark_yr, n_risk)

## 연구자 지시(2026-08-09): "Number at risk" 를 떠 있는 제목이 아니라
## Male / Female 과 같은 행 라벨 열에 정렬한다. 빈 헤더 행을 만들어 y축 라벨로 처리.
risk_rows <- c("Women", "Men", "Number at risk")   # 아래 → 위 순서
risk_tab2 <- risk_tab %>%
  mutate(row = factor(as.character(sex), levels = risk_rows))

pRisk <- ggplot(risk_tab2, aes(x = landmark_yr, y = row)) +
  geom_text(aes(label = format(n_risk, big.mark=",")), size = 2.35) +
  scale_x_continuous(breaks = 0:5, limits = c(-0.15, 6.6), expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(limits = risk_rows, drop = FALSE) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_kcto(base_size = 8) +
  theme(panel.grid = element_blank(), panel.grid.major.x = element_blank(),
        # limits 순서(Female, Male, Number at risk)에 맞춰 face 벡터를 지정
        # 연구자 요청: "Number at risk" 제목을 숫자 크기에 맞추고 세 행 라벨 크기 통일
        axis.text.y = element_text(face = c("bold", "bold", "plain"),
                                   size = rel(1.0),
                                   colour = c("grey15", "grey15", "grey25")),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(-2, 8, 2, 3))

## ---- Panel B: within-age-band hazard ratio (forest) ----------------------
## NOTE(Lead, 2026-08-10, 연구자 요청): x축을 5년 누적발생률(%)에서 **위험비(HR)** 로 전환.
## 오른쪽에 HR 텍스트를 따로 붙이던 구성 대신, HR 자체를 x축에 두어 forest 로 그린다.
## (연령대별 5년 누적발생률 자체는 Figure 2B 가 이미 막대로 보여주므로 중복도 해소된다.)
t6f_hr <- t26 %>%
  distinct(age_stratum, within_stratum_HR, within_stratum_p) %>%
  mutate(age_stratum_f = factor(age_stratum,
                                levels = rev(c("<65 years", "65-74 years", ">=75 years")),
                                labels = rev(c("<65 y", "65-74 y", "≥75 y"))),
         hr = as.numeric(sub("^([0-9.]+).*$", "\\1", within_stratum_HR)),
         lo = as.numeric(sub("^.*\\(([0-9.]+)-.*$", "\\1", within_stratum_HR)),
         hi = as.numeric(sub("^.*-([0-9.]+)\\).*$", "\\1", within_stratum_HR)),
         lab = sprintf("%s\nP=%s", within_stratum_HR, within_stratum_p),
         # ②(2026-08-27) 라벨 열: 기준선·whisker 어느 쪽과도 겹치지 않도록
         # 가장 넓은 CI 상한(1.25)보다 오른쪽에서 시작한다.
         label_x = 1.32)

pB4 <- ggplot(t6f_hr, aes(x = hr, y = age_stratum_f)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.16, linewidth = 0.55,
                color = "grey25", orientation = "y") +
  geom_point(size = 2.4, shape = 22, fill = "grey25", color = "grey25") +
  geom_text(aes(x = label_x, label = lab), hjust = 0, size = 2.0, color = "grey15", lineheight = 0.85) +
  scale_x_log10(limits = c(0.58, 2.90), breaks = c(0.75, 1.0, 1.5, 2.0),
                labels = c("0.75", "1.00", "1.50", "2.00")) +
  labs(x = "Hazard ratio for MACCE, female vs male (log scale)", y = "Age at index procedure",
       title = "B. Hazard ratio for MACCE within each age band") +
  theme_kcto(base_size = 8) +
  theme(panel.grid.major.x = element_blank(),
        axis.line.y  = element_line(color = "grey20", linewidth = 0.35),
        axis.ticks.y = element_line(color = "grey20", linewidth = 0.3),
        axis.title.y = element_text(margin = margin(r = 0)),
        axis.title.x = element_text(margin = margin(t = 1.5, b = 0)),
        plot.margin = margin(4, 8, 0, 4))

###############################################################################
## FIGURE 2 (재구성) — 연령: 연령대별 5년 누적발생률 + 연령대내 HR
##
## Lead 결정(2026-08-10, 연구자 승인): JACC:CI 저자지침은 Central Illustration 이
## 본문 figure 를 복제하지 못하게 한다. 제목이 약속하는 두 축(clinical outcomes /
## procedural intensity)이 곧 구 Figure 2A·3B·4A 의 내용이므로, 복제를 피하려면
## 그림을 "옮겨야" 한다. 지침은 핵심 figure 를 CI 로 승격하는 것을 허용한다.
##   - 구 Figure 2A (staged adjustment forest) -> Central Illustration
##   - 구 Figure 3B (복잡도 층별 강도 점수)     -> Central Illustration
##   - 구 Figure 4A (KM 곡선 + number at risk)  -> Central Illustration
## 남은 연령 관련 두 패널을 하나의 Figure 2 로 묶는다.
## pA / pB3 / pA4 / pRisk 는 CI 로 이동한 패널의 출처로 정의만 남겨둔다(저장 안 함).
###############################################################################
## 연구자 지시(2026-08-12): ① A-B 사이 1줄 간격(pB4 상단 여백) ② x축 단축 → 183→140mm
## 연구자 지시(2026-08-28, 검토 메모 id 30): Panel B(연령대내 HR forest, pB4) 삭제 →
## Figure 2 는 연령대별 5년 누적발생률 단일 패널. pB4 는 출처 기록으로 정의만 남긴다
## (전기간 단일 HR 표시가 PH 위반 원칙과 긴장 + 메시지는 Panel A·본문이 이미 전달).
p2 <- pB + theme(plot.margin = margin(3, 8, 3, 3))
save_fig(p2, file.path(FIGURES_DIR, "Figure2_age_attenuation"),
         width_mm = 140, height_mm = 85)
cat("Figure 2 saved.\n")
