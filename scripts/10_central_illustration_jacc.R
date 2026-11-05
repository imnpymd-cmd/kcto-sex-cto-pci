###############################################################################
# K-CTO Registry 2.0 — Central Illustration, JACC:CI house style
#
# 실행:  Rscript scripts/10_central_illustration_jacc.R
# 출력:  output/figures/CentralIllustration.tiff  (1200 dpi, 채택본; PNG 중단 — 연구자 확정 2026-08-12)
#        BUILD_VARIANTS=TRUE 이면 검토 시안 v1/v2/v4 도 함께 생성
#
# 연구자가 제공한 실제 게재 Central Illustration 9편(JACC Cardiovasc Interv
# 2024-2026)에서 관찰된 house style 을 그대로 따른다:
#   1) 전체를 얇은 빨간 테두리로 감싼다
#   2) 최상단 연회색 머리띠: "CENTRAL ILLUSTRATION:" (빨강 볼드) + 제목(남색 볼드)
#   3) 진남색 배너에 흰 글씨로 연구 제목 + N
#   4) 섹션 구분은 중간 파랑 바 + 흰 볼드 가운데 정렬
#   5) 내용 상자는 연파랑 배경 + 진한 글씨
#   6) 플롯 패널 배경은 아주 연한 파랑
#   7) 계열색 2종(파랑/주황-빨강) — Okabe-Ito 유지(색맹 안전)
#   8) 맨 아래 진회색 상자에 흰 불릿 take-home
#   ※ 게재본에 있는 하단 인용줄(Author, et al. JACC...)은 투고본에는 넣지 않는다.
#
# 모든 수치는 output/tables/*.csv 에서 직접 읽는다.
# 금지 문구(CLAUDE.md §1) 준수: "widened"/"every level"/"equal"/인과 언어 등.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(stringr); library(forcats)
})

PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
source(file.path(PROJ, "scripts", "_theme.R"))
rt <- function(f) read.csv(file.path(TABLES_DIR, f), stringsAsFactors = FALSE, check.names = FALSE)

## 채택본만 생성한다(연구자 승인 2026-08-10). 검토 과정에서 만든 v1/v2/v4 는
## legacy/central_illustration_variants_2026-08-10/ 로 옮겼다. TRUE 로 두면 함께 재생성된다.
BUILD_VARIANTS <- FALSE

## ---- JACC house palette -----------------------------------------------------
JRED    <- "#C8102E"   # "CENTRAL ILLUSTRATION:" 및 외곽선
JNAVY   <- "#1F4E79"   # 연구 배너 / 제목 글씨
JBAR    <- "#5B9BD5"   # 섹션 바
JBOX    <- "#DCE9F7"   # 연파랑 내용 상자
JPANEL  <- "#F3F8FD"   # 플롯 패널 배경
JSTRIP  <- "#EAF1F6"   # 최상단 머리띠
JDARK   <- "#4D4D4D"   # take-home 상자
SER_B   <- unname(okabe_ito["blue"])        # 남성
SER_R   <- unname(okabe_ito["vermillion"])  # 여성

theme_jacc <- function(base = 8) {
  theme_minimal(base_size = base, base_family = "Helvetica") +
    theme(panel.background  = element_rect(fill = JPANEL, color = NA),
          panel.grid.minor  = element_blank(),
          panel.grid.major  = element_blank(),
          axis.line         = element_line(color = "grey25", linewidth = 0.35),
          axis.ticks        = element_line(color = "grey25", linewidth = 0.3),
          axis.text         = element_text(size = rel(0.95), color = "grey15"),
          axis.title        = element_text(size = rel(0.95), color = "grey10"),
          plot.title        = element_blank(),
          plot.margin       = margin(3, 5, 3, 5),
          legend.position   = "none")
}

## ---- building blocks --------------------------------------------------------
strip_header <- function(title, size = 3.35) {
  ggplot() + xlim(0, 1) + ylim(0, 1) +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = JSTRIP, color = NA) +
    annotate("text", x = 0.008, y = 0.5, label = "CENTRAL ILLUSTRATION:", hjust = 0, vjust = 0.5,
             color = JRED, fontface = "bold", size = size) +
    annotate("text", x = 0.248, y = 0.5, label = title, hjust = 0, vjust = 0.5,
             color = JNAVY, fontface = "bold", size = size) +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
}
bar <- function(txt, fill = JBAR, col = "white", size = 3.0, hj = 0.5, face = "bold") {
  ggplot() + xlim(0, 1) + ylim(0, 1) +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0.06, ymax = 0.94, fill = fill, color = NA) +
    annotate("text", x = hj, y = 0.5, label = txt, hjust = hj, vjust = 0.5,
             color = col, fontface = face, size = size, lineheight = 1.12) +
    theme_void() + theme(plot.margin = margin(1, 2, 1, 2))
}
pale_box <- function(txt, size = 2.7, face = "plain") {
  ggplot() + xlim(0, 1) + ylim(0, 1) +
    annotate("rect", xmin = 0.01, xmax = 0.99, ymin = 0.05, ymax = 0.95, fill = JBOX, color = NA) +
    annotate("text", x = 0.5, y = 0.5, label = txt, hjust = 0.5, vjust = 0.5,
             color = "grey10", fontface = face, size = size, lineheight = 1.3) +
    theme_void() + theme(plot.margin = margin(2, 3, 2, 3))
}
takehome <- function(bullets, size = 2.5, lineheight = 1.45) {
  # ③(2026-08-27) 줄바꿈된 연속 줄에 hanging indent.
  #   전에는 왼쪽 끝에 붙어 별도 항목처럼 읽혔다.
  bl  <- vapply(bullets, function(b) paste0("• ", gsub("\n", "\n   ", b, fixed = TRUE)),
                character(1), USE.NAMES = FALSE)
  txt <- paste(bl, collapse = "\n")
  ggplot() + xlim(0, 1) + ylim(0, 1) +
    annotate("rect", xmin = 0.004, xmax = 0.996, ymin = 0.04, ymax = 0.96, fill = JDARK, color = NA) +
    annotate("text", x = 0.018, y = 0.5, label = txt, hjust = 0, vjust = 0.5,
             color = "white", size = size, lineheight = lineheight) +
    theme_void() + theme(plot.margin = margin(2, 3, 2, 3))
}

## ---- data -------------------------------------------------------------------
t1 <- rt("table1_baseline_by_sex.csv"); t3 <- rt("table3_procedural_by_sex.csv")
t7 <- rt("table7_adjusted_HR_all_outcomes.csv"); t13 <- rt("table13_adjusted_procedural_inhospital.csv")
t14 <- rt("table14_within_center_management.csv"); t15c <- rt("table15c_intensity_score_robustness.csv")
s7 <- rt("tableS7_analysis_diagnostics.csv"); t5b <- rt("table5b_km_summary_by_sex.csv")
t5e <- rt("table5e_km_curve_coordinates.csv"); t6f <- rt("table6f_age_stratified_DMS.csv")

gv <- function(p) s7$value[grepl(p, s7$item)][1]
N_TOT <- gv("ANALYSIS COHORT"); N_M <- s7$value[s7$item == "    Male"]; N_F <- s7$value[s7$item == "    Female"]
N_CTR <- as.integer(gv("Contributing centers")); pct_f <- round(100 * N_F / N_TOT, 1)

stage <- t7 %>% filter(outcome == "DMS (death, MI, or stroke) - PRIMARY",
                       model %in% c("Unadjusted", "Age-adjusted", "Fully adjusted")) %>%
  mutate(st = factor(model, levels = c("Unadjusted", "Age-adjusted", "Fully adjusted"),
                     labels = c("Unadjusted", "+ Age", "+ Full"))) %>% arrange(st)
HRc <- stage$HR_fmt[1]; HRa <- stage$HR_fmt[2]; HRf <- stage$HR_fmt[3]

a75 <- t1 %>% filter(variable == "age75")
p75 <- c(as.numeric(sub("^.*\\(([0-9.]+)\\)$", "\\1", a75$male)),
         as.numeric(sub("^.*\\(([0-9.]+)\\)$", "\\1", a75$female)))
agev <- t1 %>% filter(variable == "Age")

sc <- t15c %>% filter(outcome == "Intensity score 0-3, excluding potent P2Y12 at discharge",
                      model == "Linear mixed model, centre random intercept")
sc_all <- sc %>% filter(stratum == "Full analysis cohort", n == 5879) %>% slice(1)
gapd <- bind_rows(sc_all %>% mutate(lab = "Overall"),
                  sc %>% filter(grepl("complexity stratum", stratum),
                                covariate_set == "Model 4 covariates minus centre volume (absorbed by centre)") %>%
                    mutate(lab = recode(stratum,
                      "Baseline complexity stratum: J-CTO score <2 (contrast)" = "J-CTO <2",
                      "Baseline complexity stratum: J-CTO score >=2"            = "J-CTO ≥2",
                      "Baseline complexity stratum: occlusion length >=20 mm"   = "Length\n≥20 mm"))) %>%
  mutate(lab = factor(lab, levels = rev(c("Overall", "J-CTO <2", "J-CTO ≥2", "Length\n≥20 mm"))))

bilat <- t14 %>% filter(outcome == "Bilateral injection",
                        model == "(2) Within-centre: random intercept for centre (full Laplace)") %>% slice(1)
tech  <- t13 %>% filter(outcome == "Technical success of CTO PCI", grepl("Adjusted", model)) %>% slice(1)
lrp   <- t5b$logrank_p_fmt[t5b$outcome_var == "MACCE_pub"][1]   # v2: 0.330

## 5년 시점 누적발생률 (연구자 지시 2026-08-12: KM 곡선 위·아래에 직접 표기)
## 근거: table5a — 남 15.65 → "15.7%", 여 22.4 → "22.4%" (원고 Results 인쇄값과 동일)
t5a  <- rt("table5a_km_timepoints_by_sex.csv")
km5  <- t5a %>% filter(outcome_var == "MACCE_pub", timepoint == "5-year") %>% distinct(sex, .keep_all = TRUE)
inc5_m <- km5$cum_incidence_pct[km5$sex == "Male"]
inc5_f <- km5$cum_incidence_pct[km5$sex == "Female"]
stopifnot(length(inc5_m) == 1, length(inc5_f) == 1)

## v2 landmark 창별 값 + RMST (table22 / table8c)
lm22 <- rt("table22_landmark_partition.csv")
.lmget <- function(mdl, w) {
  r <- lm22[lm22$outcome_key == "MACCE_pub" & startsWith(lm22$window, w) & lm22$model == mdl, ][1, ]
  r
}
lm_rows <- rbind(.lmget("Model 4 (re-fitted within window)", "0-1"),
                 .lmget("Model 4 (re-fitted within window)", "1-3"),
                 .lmget("parsimonious (age, DM, CKD, AF)",  ">3"))
lm_rows$wlab <- factor(c("0-1 year", "1-3 years", ">3 years"),
                       levels = rev(c("0-1 year", "1-3 years", ">3 years")))
lm_crude <- .lmget("unadjusted", ">3")$HR_fmt
lm_age   <- .lmget("age only",  ">3")$HR_fmt
t8c <- rt("table8c_rmst_5yr.csv")
rmst_fmt <- t8c$adj_rmst_diff_fmt[t8c$outcome_key == "MACCE_pub"]

grab <- function(p) { r <- t3[grepl(p, paste0(t3$label, t3$variable)), ][1, ]
  c(m = as.numeric(sub("^.*\\(([0-9.]+)\\)$", "\\1", r$male)),
    f = as.numeric(sub("^.*\\(([0-9.]+)\\)$", "\\1", r$female))) }
usg <- bind_rows(data.frame(item = "Bilateral injection", t(grab("Bilateral injection"))),
                 data.frame(item = "Intravascular ultrasound", t(grab("IVUS use"))),
                 data.frame(item = "Technical success", t(grab("Technical success of CTO PCI")))) %>%
  pivot_longer(c(m, f), names_to = "sex", values_to = "pct") %>%
  mutate(sex = sex_factor(ifelse(sex == "m", "Male", "Female")))

km <- t5e %>% filter(outcome_var == "MACCE_pub") %>%   # v2: primary = MACCE
  mutate(sex = sex_factor(sex), yr = time_days / 365.25)

## number-at-risk — Figure 4 가 CI 로 흡수되면서 위험집합 표도 함께 옮긴다.
## 좌표는 KM 과 동일한 table5e 스텝에서 읽는다(§2-3: 단순 뺄셈으로 검산하지 말 것).
risk_ci <- bind_rows(lapply(sex_display, function(sx) {
  sub <- km %>% filter(sex == sx) %>% arrange(time_days)
  do.call(rbind, lapply(c(0, 365.25, 730.5, 1095.75, 1461, 1826.25), function(d) {
    s2 <- sub[sub$time_days <= max(d, min(sub$time_days)), ]
    if (nrow(s2) == 0) s2 <- sub[1, ]
    s2[which.max(s2$time_days), ]
  })) %>% mutate(yr_lab = 0:5)
})) %>%
  transmute(row = factor(as.character(sex), levels = c("Women", "Men", "No. at risk")),
            yr_lab, n_risk)

## 행 라벨은 축 텍스트가 아니라 패널 왼쪽에 매다는 geom_text 로 그린다.
## 축 텍스트로 두면 patchwork 가 KM 의 y축 텍스트 열을 이 라벨 폭만큼 넓혀
## 타이틀-축 간격이 벌어진다 (연구자 지적 2026-08-12).
risk_lab <- data.frame(row = factor(c("Women", "Men", "No. at risk"),
                                    levels = c("Women", "Men", "No. at risk")))
p_risk <- ggplot(risk_ci, aes(yr_lab, row)) +
  geom_text(aes(label = format(n_risk, big.mark = ",")), size = 2.1, colour = "grey15") +
  geom_text(data = risk_lab[1:2, , drop = FALSE], aes(x = -0.26, y = row, label = row),
            hjust = 1, size = 2.1, colour = "grey15", fontface = "bold") +
  geom_text(data = risk_lab[3, , drop = FALSE], aes(x = -0.26, y = row, label = row),
            hjust = 1, size = 2.1, colour = "grey30") +
  scale_x_continuous(breaks = 0:5, limits = c(-0.26, 5.30), expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(limits = c("Women", "Men", "No. at risk"), drop = FALSE,
                   expand = expansion(add = c(0.30, 0.30))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +   # x 타이틀은 p_km 으로 이동 (연구자 지시 2026-08-12)
  theme_jacc() +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(0, 5, 1, 10))

## ---- panels -----------------------------------------------------------------
sex_leg <- c(Men   = sprintf("Men (n = %s)",   format(N_M, big.mark = ",")),
             Women = sprintf("Women (n = %s)", format(N_F, big.mark = ",")))
p_km <- ggplot(km, aes(yr, cum_incidence_pct, color = sex, linetype = sex)) +
  geom_step(linewidth = 0.6) +
  # NOTE(Lead): 미보정 HR 은 옆 패널(동일표본 n=6,259)과 표본이 달라 CI 가 어긋난다.
  # 한 그림 안에 두 값이 공존하지 않도록 여기서는 log-rank P 만 표시한다.
  # (연구자 지시 2026-08-12: log-rank 는 하단 배치, 5년 누적발생률은 곡선 위·아래 직접 표기)
  annotate("text", x = 3, y = 8.0, hjust = 0.5, size = 2.5, color = "grey20",
           label = paste0("Log-rank P", ifelse(grepl("^<", lrp), lrp, paste0("=", lrp)))) +
  annotate("text", x = 5.0, y = inc5_f + 2.2, hjust = 1, size = 2.4, fontface = "bold",
           color = SER_R, label = sprintf("%.1f%%", inc5_f)) +
  annotate("text", x = 5.0, y = inc5_m - 3.4, hjust = 1, size = 2.4, fontface = "bold",
           color = SER_B, label = sprintf("%.1f%%", inc5_m)) +
  scale_color_manual(values = c(Men = SER_B, Women = SER_R), labels = sex_leg) +
  scale_linetype_manual(values = c(Men = "solid", Women = "22"), labels = sex_leg) +
  scale_x_continuous(breaks = 0:5, limits = c(-0.26, 5.30), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 38), expand = expansion(mult = c(0, 0))) +
  labs(x = "Years since index CTO PCI", y = "Cumulative incidence of\nMACCE (%)",
       title = "Unadjusted cumulative incidence") +
  theme_jacc() + theme(plot.title = element_text(size = 7.6, face = "bold", colour = JNAVY,
                                                 hjust = 0.5, margin = margin(b = 4)),
                       axis.title.y = element_text(margin = margin(r = 1)),
                       axis.title.x = element_text(margin = margin(t = 2)),
                       legend.position = "inside",
                       legend.position.inside = c(0.03, 0.98),
                       legend.justification.inside = c(0, 1),
                       legend.title = element_blank(),
                       legend.text = element_text(size = 6.3, face = "bold", color = "grey15"),
                       legend.key.height = unit(8, "pt"), legend.key.width = unit(14, "pt"),
                       legend.background = element_blank(),
                       legend.spacing.y = unit(1, "pt"),
                       plot.margin = margin(3, 5, 1, 29))

## v2: staged 감쇠 패널 → landmark 창별 forest (연구자 확정 2026-08-14)
p_stage <- ggplot(lm_rows, aes(HR, wlab)) +
  # ③(2026-08-27) 파선을 추정치 행 구간으로 제한. 전에는 패널 전체를 가로질러
  #   아래 RMST 주석의 글자를 관통했다.
  annotate("segment", x = 1, xend = 1, y = 0.62, yend = 3.62,
           linetype = "dashed", color = "grey45", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.14, linewidth = 0.6,
                color = SER_B, orientation = "y") +
  geom_point(size = 3.2, color = SER_B) +
  geom_text(aes(label = HR_fmt), vjust = -1.5, size = 2.25, fontface = "bold", color = "grey10") +
  annotate("text", x = 1, y = 0.38, hjust = 0.5, size = 2.35, color = "grey25",
           label = sprintf("Adjusted 5-year restricted-mean difference +%s days", rmst_fmt)) +
  scale_x_log10(limits = c(0.60, 1.75), breaks = c(0.75, 1.0, 1.5),
                labels = c("0.75", "1.00", "1.50")) +
  scale_y_discrete(expand = expansion(add = c(1.0, 0.55))) +
  coord_cartesian(clip = "off") +
  labs(x = "Adjusted hazard ratio, female vs male",
       y = NULL, title = "Landmark-partitioned windows") +
  theme_jacc() + theme(axis.text.y = element_text(face = "bold"),
                       plot.title = element_text(size = 7.6, face = "bold", colour = JNAVY,
                                                 hjust = 0.5, margin = margin(b = 4)),
                       plot.margin = margin(3, 33, 3, 5))

p_ageband <- (function() {
  d <- t6f %>% mutate(band = factor(recode(age_stratum, "<65 years" = "<65 y",
                                           "65-74 years" = "65-74 y", ">=75 years" = "≥75 y"),
                                    levels = c("<65 y", "65-74 y", "≥75 y")),
                      sex = sex_factor(sex))
  ggplot(d, aes(band, cum_inc_5yr_pct, fill = sex)) +
    geom_col(position = position_dodge(width = 0.68), width = 0.6) +
    geom_text(aes(label = sprintf("%.1f", cum_inc_5yr_pct)), position = position_dodge(width = 0.68),
              vjust = -0.4, size = 2.1) +
    scale_fill_manual(values = c(Men = SER_B, Women = SER_R)) +
    scale_y_continuous(limits = c(0, 40), expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Age at index procedure", y = "5-year incidence (%)") +
    theme_jacc()
})()

# 연구자 지시(2026-08-28, 검토 메모 id 30): Figure 3B 삭제 후 이 패널이 층별 forest 의
# 유일한 그림이 되므로 95% CI 와 층별 n 을 함께 인쇄한다 (평균차 패널이라 OR 축 눈금은 해당 없음).
gap_lev <- levels(gapd$lab)
gap_lev_n <- setNames(sprintf("%s\nn=%s", sub("\n", " ", gap_lev),
                              formatC(as.integer(gapd$n[match(gap_lev, as.character(gapd$lab))]),
                                      big.mark = ",", format = "d")), gap_lev)
gapd$lab_n <- factor(gap_lev_n[as.character(gapd$lab)], levels = unname(gap_lev_n))

p_gap <- ggplot(gapd, aes(estimate, lab_n)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.15, linewidth = 0.6, color = SER_R) +
  geom_point(size = 2.9, color = SER_R, shape = 18) +
  geom_text(aes(label = sprintf("%.3f (%.3f to %.3f)", estimate, lo, hi)),
            vjust = -1.35, size = 1.95, fontface = "bold", color = "grey10") +
  scale_x_continuous(limits = c(-0.235, 0.045)) +
  labs(x = "Adjusted mean score difference, female vs male", y = NULL) +
  theme_jacc() + theme(axis.text.y = element_text(face = "bold", lineheight = 0.9),
                       plot.margin = margin(3, 33, 3, 5))

usg$item2 <- factor(recode(usg$item,
                           "Bilateral injection"      = "Bilateral\ninjection",
                           "Intravascular ultrasound" = "Intravascular\nultrasound",
                           "Technical success"        = "Technical\nsuccess"),
                    levels = c("Bilateral\ninjection", "Intravascular\nultrasound",
                               "Technical\nsuccess"))
p_usage <- ggplot(usg, aes(pct, fct_rev(item2),
                           fill = fct_rev(sex))) +
  geom_col(position = position_dodge(width = 0.66), width = 0.58) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.66),
            hjust = -0.15, size = 2.1) +
  scale_fill_manual(values = c(Men = SER_B, Women = SER_R)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Patients (%)", y = NULL) +
  theme_jacc() + theme(axis.text.y = element_text(face = "bold"),
                       legend.position = "inside",
                       legend.position.inside = c(0.97, 0.98),
                       legend.justification.inside = c(1, 1),
                       legend.title = element_blank(),
                       legend.text = element_text(size = 6.3, face = "bold", color = "grey15"),
                       legend.key.size = unit(8, "pt"),
                       legend.background = element_blank(),
                       plot.margin = margin(3, 15, 3, 28))

legend_row <- ggplot() + xlim(0, 1) + ylim(0, 1) +
  annotate("rect", xmin = 0.335, xmax = 0.359, ymin = 0.34, ymax = 0.66, fill = SER_B, color = NA) +
  annotate("text", x = 0.369, y = 0.5, hjust = 0, size = 2.75, color = "grey15", fontface = "bold",
           label = sprintf("Men (n = %s)", format(N_M, big.mark = ","))) +
  annotate("rect", xmin = 0.545, xmax = 0.569, ymin = 0.34, ymax = 0.66, fill = SER_R, color = NA) +
  annotate("text", x = 0.579, y = 0.5, hjust = 0, size = 2.75, color = "grey15", fontface = "bold",
           label = sprintf("Women (n = %s)", format(N_F, big.mark = ","))) +
  theme_void() + theme(plot.margin = margin(0, 2, 0, 2))

TITLE_SHORT <- "Sex Differences in Outcomes and Procedural Intensity After CTO PCI"
BANNER <- sprintf("Sex Differences in Clinical Outcomes and Procedural Intensity After Chronic Total Occlusion PCI:\nK-CTO Nationwide Retrospective Registry, %d centers, 2014-2023, N = %s",
                  N_CTR, format(N_TOT, big.mark = ","))

# 연구자 지시(2026-08-28): 핵심 메시지 3문장으로 압축, 글자 크기 = 패널 x축 타이틀(7.6pt).
# CI 붙은 전체값 대신 점추정치만 인쇄(전체값은 본문 Table 3·4 / Table 5). hedge "a late
# difference is not excluded" 는 QA 확정 문안이라 보존.
BULLETS <- c(
  sprintf("Women were %s of the cohort and %.0f years older; after adjustment, no MACCE difference was detected in any follow-up window, and the crude late-excess hazard ratio fell from %s to %s with age adjustment alone (a late difference is not excluded).",
          paste0(pct_f, "%"),
          as.numeric(sub("^([0-9.]+).*$", "\\1", agev$female)) - as.numeric(sub("^([0-9.]+).*$", "\\1", agev$male)),
          sub("^([0-9.]+).*$", "\\1", lm_crude), sub("^([0-9.]+).*$", "\\1", lm_age)),
  sprintf("Under the same adjustment, women received a lower procedural intensity score, and bilateral contrast injection was the only individual management strategy that remained lower (adjusted OR %s).",
          sub("^([0-9.]+).*$", "\\1", bilat$est_fmt)),
  "Technical success did not differ detectably by sex.")

W <- function(g) patchwork::wrap_elements(full = g)

frame <- function(p) p + plot_annotation(theme = theme(
  plot.background = element_rect(fill = "white", color = JRED, linewidth = 1.1),
  plot.margin = margin(4, 5, 4, 5)))

###############################################################################
## v1 — Outcomes / Procedural intensity, two stacked sections (NOTION-3 형식)
###############################################################################
if (BUILD_VARIANTS) {
v1 <- W(strip_header(TITLE_SHORT)) / W(bar(BANNER, fill = JNAVY, size = 3.0)) /
  W(legend_row) /
  W(bar("Clinical Outcomes")) / (p_km | p_stage) /
  W(bar("Procedural Intensity and Technical Success")) / (p_gap | p_usage) /
  W(takehome(str_wrap(BULLETS, 112))) +
  plot_layout(heights = c(0.40, 0.58, 0.20, 0.24, 2.15, 0.24, 2.15, 1.70))
save_fig(frame(v1), file.path(FIGURES_DIR, "CentralIllustration_jacc_v1"),
         width_mm = 210, height_mm = 208)
cat("JACC-style v1 saved.\n")
}

###############################################################################
## v2 — Study Population | Results 두 열 (SGLT2 형식)
###############################################################################
pop_box <- pale_box(sprintf(
  "K-CTO Registry 2.0\n%d Korean centers, 2014-2023\n\n%s patients undergoing CTO PCI\nMen %s   |   Women %s (%.1f%%)\n\nMedian age %s vs %s years (male vs female)",
  N_CTR, format(N_TOT, big.mark = ","), format(N_M, big.mark = ","),
  format(N_F, big.mark = ","), pct_f,
  sub("^([0-9.]+).*$", "\\1", agev$male), sub("^([0-9.]+).*$", "\\1", agev$female)), size = 2.6)

if (BUILD_VARIANTS) {
v2 <- W(strip_header(TITLE_SHORT)) / W(bar(BANNER, fill = JNAVY, size = 3.0)) /
  W(legend_row) /
  ((bar("Study Population") / pop_box / bar("Age Composition") / p_ageband +
      plot_layout(heights = c(0.3, 1.5, 0.3, 1.9))) |
     (bar("Outcomes: Staged Adjustment") / p_stage / bar("Procedural Intensity by Complexity") / p_gap +
        plot_layout(heights = c(0.3, 1.5, 0.3, 1.9)))) /
  W(takehome(str_wrap(BULLETS, 112))) +
  plot_layout(heights = c(0.40, 0.58, 0.20, 4.4, 1.70))
save_fig(frame(v2), file.path(FIGURES_DIR, "CentralIllustration_jacc_v2"),
         width_mm = 210, height_mm = 186)
cat("JACC-style v2 saved.\n")
}

###############################################################################
## v3 — KM 중심 (Sakuma / NOTION-3 하이브리드): 곡선을 크게, 아래 3패널
###############################################################################
## 연구자 지시(2026-08-12): 독립 범례 행 제거 — 범례는 p_km·p_usage 패널 내부로 이동.
## KM x 타이틀이 패널로 들어가며 위험표 행은 축소. 세로 배분 재조정.
v3 <- W(strip_header(TITLE_SHORT)) + W(bar(BANNER, fill = JNAVY, size = 3.0)) +
  W(bar("Primary Outcome: MACCE — Death, MI, Stroke, or Any Revascularization")) +
  p_km + p_stage + p_risk +
  W(bar("Process of Care and Technical Success")) +
  W(p_usage) + W(p_gap) +   # wrap: 상단 패널과의 가로 폭 정렬 해제 (연구자 지시 2026-08-12)
  # 연구자 지시(2026-08-28): 글자 크기 = x축 타이틀 7.6pt, 3문장 압축.
  # 2차 지시: 오른쪽 여백 제거 — 줄바꿈 110 → 140자(박스 폭 채움), 5줄이 되어 높이 1.55 → 1.35.
  W(takehome(str_wrap(BULLETS, 140), size = 7.6 / ggplot2::.pt, lineheight = 1.32)) +
  plot_layout(design = "AA\nBB\nCC\nDE\nFE\nGG\nHI\nJJ",
              heights = c(0.40, 0.58, 0.26, 2.42, 0.34, 0.26, 1.90, 1.35))
## 채택본(연구자 승인 2026-08-10). 구 Figure 2A / 3B / 4A 를 흡수한 정본이다.
save_fig(frame(v3), file.path(FIGURES_DIR, "CentralIllustration"),
         width_mm = 210, height_mm = 210)
cat("Central Illustration (채택본) saved.\n")

###############################################################################
## v4 — 본문 figure 와 중복되지 않는 구성
##
## JACC:CI 저자지침: Central Illustration "must not duplicate content from other
## figures in the manuscript". v1~v3 은 KM 곡선(Figure 4A), staged HR 플롯
## (Figure 2A), 복잡도별 평균차 forest(Figure 3B)를 그대로 재사용하므로 저촉된다.
## v4 는 어느 본문 figure 에도 그림으로 등장하지 않는 재료만 쓴다:
##   - 성별 연령 구성           (Table 1 계열 — figure 없음)
##   - staged HR 수치 도식       (플롯이 아닌 숫자 상자)
##   - 미보정 시행률            (Table 3 — figure 없음)
##   - 개별 지표 센터내 보정 OR (Table 13b — figure 없음)
###############################################################################

if (BUILD_VARIANTS) {
t13b <- rt("table13b_management_multiplicity.csv")

## 연령 구성 — 전체 분석코호트 N=6,777 (Table 1 의 >=75세 18.2 / 43.6 과 일치)
.d <- readRDS(file.path(PROJ, "data", "processed", "kcto_sex_analysis.rds"))
.d <- subset(.d, flag_cto_pci & index_yr >= 2014 & index_yr <= 2023)
stopifnot(nrow(.d) == 6777)
ageb <- .d %>%
  mutate(band = cut(Age, c(-Inf, 64.999, 74.999, Inf),
                    labels = c("<65 y", "65-74 y", "≥75 y")),
         sex = sex_factor(ifelse(female_bi == 1, "Female", "Male"))) %>%
  count(sex, band) %>% group_by(sex) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
stopifnot(all(abs(ageb$pct[ageb$band == "≥75 y"] - c(18.2, 43.6)) < 0.1))
rm(.d)

## 미보정 시행률 — Table 3 (전체 분석코호트). retrograde attempt 는 수준값의 여집합.
pctv <- function(s) as.numeric(sub("^.*\\(([0-9.]+)\\)$", "\\1", s))
row3 <- function(p) t3[grepl(p, t3$label, fixed = TRUE), ][1, ]
na_m <- pctv(row3("  Not attempted")$male); na_f <- pctv(row3("  Not attempted")$female)
IND <- c("Bilateral\ninjection", "Intravascular\nultrasound",
         "Retrograde\nattempt", "Potent P2Y12\nat discharge")
use4 <- data.frame(
  item = factor(IND, levels = IND),
  Male   = c(pctv(row3("Bilateral injection, n (%)")$male),
             pctv(row3("IVUS use, n (%)")$male),
             round(100 - na_m, 1),
             pctv(row3("Potent P2Y12 inhibitor (prasugrel or ticagrelor), n (%)")$male)),
  Female = c(pctv(row3("Bilateral injection, n (%)")$female),
             pctv(row3("IVUS use, n (%)")$female),
             round(100 - na_f, 1),
             pctv(row3("Potent P2Y12 inhibitor (prasugrel or ticagrelor), n (%)")$female))) %>%
  pivot_longer(c(Male, Female), names_to = "sex", values_to = "pct") %>%
  mutate(sex = sex_factor(sex))

## 개별 지표 센터내(랜덤절편) 보정 OR — 네 지표 모두 동일 모형
or4 <- t13b %>%
  filter(outcome %in% c("Bilateral injection", "Intravascular ultrasound use",
                        "Retrograde approach attempted",
                        "Potent P2Y12 inhibitor at discharge (among DAPT patients)")) %>%
  mutate(item = factor(recode(outcome,
           "Bilateral injection" = "Bilateral\ninjection",
           "Intravascular ultrasound use" = "Intravascular\nultrasound",
           "Retrograde approach attempted" = "Retrograde\nattempt",
           "Potent P2Y12 inhibitor at discharge (among DAPT patients)" =
             "Potent P2Y12\nat discharge"), levels = IND),
         lab = sprintf("%.2f (%.2f-%.2f)", ri_OR, ri_lo, ri_hi)) %>%
  arrange(item)

## ---- v4 panels --------------------------------------------------------------
ttl <- function(p, s) p + labs(title = s) +
  theme(plot.title = element_text(size = 7.6, face = "bold", colour = JNAVY,
                                  hjust = 0.5, margin = margin(b = 4)))

p_age <- ttl(ggplot(ageb, aes(band, pct, fill = sex)) +
  geom_col(position = position_dodge(width = 0.70), width = 0.62) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.70),
            vjust = -0.45, size = 2.2) +
  scale_fill_manual(values = c(Men = SER_B, Women = SER_R)) +
  scale_y_continuous(limits = c(0, 62), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Age at index procedure", y = "Patients within each sex (%)") +
  theme_jacc(), "Age composition")

## staged HR — 플롯이 아니라 숫자 도식 (본문 Figure 2A 와 그림이 겹치지 않도록)
chain <- (function() {
  yb <- c(0.815, 0.500, 0.185); h <- 0.115
  g <- ggplot() +
    annotate("text", x = 0.5, y = 0.975, hjust = 0.5, vjust = 1, size = 2.65,
             fontface = "bold", colour = JNAVY,
             label = "Hazard ratio, female vs male")
  for (i in 1:3) {
    g <- g +
      annotate("rect", xmin = 0, xmax = 1, ymin = yb[i] - h, ymax = yb[i] + h,
               fill = JBOX, color = NA) +
      annotate("text", x = 0.07, y = yb[i], hjust = 0, size = 2.5, colour = "grey25",
               label = c("Unadjusted", "+ Age", "+ Full adjustment")[i]) +
      annotate("text", x = 0.93, y = yb[i], hjust = 1, size = 3.1, fontface = "bold",
               colour = JNAVY, label = stage$HR_fmt[i])
  }
  g +
    annotate("segment", x = 0.5, xend = 0.5, y = yb[1] - h - 0.012, yend = yb[2] + h + 0.012,
             colour = JNAVY, linewidth = 0.6,
             arrow = grid::arrow(length = unit(1.5, "mm"), type = "closed")) +
    annotate("segment", x = 0.5, xend = 0.5, y = yb[2] - h - 0.012, yend = yb[3] + h + 0.012,
             colour = JNAVY, linewidth = 0.6,
             arrow = grid::arrow(length = unit(1.5, "mm"), type = "closed")) +
    annotate("text", x = 0.53, y = (yb[1] + yb[2]) / 2, hjust = 0, size = 2.3,
             colour = "grey30", label = "98.9% of the excess\nhazard removed") +
    annotate("text", x = 0.5, y = 0.035, hjust = 0.5, size = 2.25, colour = "grey35",
             label = "Death, myocardial infarction, or stroke; n = 6,259") +
    scale_x_continuous(limits = c(0, 1), expand = expansion(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(0, 0)) +
    theme_void() + theme(plot.margin = margin(3, 5, 3, 5))
})()

p_use4 <- ttl(ggplot(use4, aes(pct, fct_rev(item), fill = fct_rev(sex))) +
  geom_col(position = position_dodge(width = 0.68), width = 0.60) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.68),
            hjust = -0.16, size = 2.1) +
  scale_fill_manual(values = c(Men = SER_B, Women = SER_R)) +
  scale_x_continuous(limits = c(0, 52), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Patients (%)", y = NULL) +
  theme_jacc() + theme(axis.text.y = element_text(face = "bold")),
  "Unadjusted use")

p_or4 <- ttl(ggplot(or4, aes(ri_OR, fct_rev(item))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ri_lo, xmax = ri_hi), width = 0.16, linewidth = 0.6, colour = JNAVY) +
  geom_point(size = 2.9, colour = JNAVY) +
  geom_text(aes(label = lab), vjust = -1.25, size = 2.1, fontface = "bold", colour = "grey10") +
  scale_x_continuous(trans = "log", breaks = c(0.6, 0.8, 1.0, 1.2),
                     limits = c(0.55, 1.30)) +
  labs(x = "Adjusted odds ratio, female vs male", y = NULL) +
  theme_jacc() + theme(axis.text.y = element_blank()),
  "Adjusted, within centre")

v4 <- W(strip_header(TITLE_SHORT)) + W(bar(BANNER, fill = JNAVY, size = 3.0)) +
  W(legend_row) +
  W(bar("Patient Age and the Primary Outcome")) +
  p_age + chain +
  W(bar("Process of Care: Unadjusted Use and Adjusted Comparison")) +
  p_use4 + p_or4 +
  W(takehome(str_wrap(BULLETS, 112))) +
  plot_layout(design = "AA\nBB\nCC\nDD\nEF\nGG\nHI\nJJ",
              heights = c(0.40, 0.58, 0.20, 0.24, 2.4, 0.24, 1.95, 1.70))
save_fig(frame(v4), file.path(FIGURES_DIR, "CentralIllustration_jacc_v4"),
         width_mm = 210, height_mm = 210)
cat("JACC-style v4 saved.\n")
}
