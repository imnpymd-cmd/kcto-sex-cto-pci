# ============================================================================
# _theme.R — Shared ggplot2 theme, colorblind-safe palette, and I/O helpers
#             for K-CTO sex-difference manuscript figures (Phase 4a, illustrator)
#
# Sourced by: scripts/07_figures.R, scripts/10_central_illustration_jacc.R
#
# Design notes:
#  - Colour palette: Okabe-Ito (colorblind-safe, printed and screen tested).
#    Male = blue (#0072B2), Female = vermillion (#D55E00) — high-contrast,
#    distinguishable under deuteranopia/protanopia/tritanopia simulation.
#  - Per instructions, sex is additionally encoded by SHAPE and LINETYPE
#    (not colour alone): Male = circle / solid line, Female = triangle / dashed
#    line, so figures remain interpretable in grayscale print.
#  - Resolution: journal_config.json critical_downstream_alerts.figure_dpi
#    requires line_art >=1200 dpi and halftone/combination >=600 dpi to
#    satisfy JACC:CI + Circulation:CVI simultaneously. All ggplot output here
#    is vector/line-art in nature (points, lines, text) -> exported at
#    1200 dpi TIFF (LZW), which exceeds both thresholds. PNG 산출은 중단
#    (연구자 확정 2026-08-12: 제출용 TIFF + 편집용 PPTX 만 유지).
#  - Fonts: sans-serif (Helvetica/Arial equivalent), base size tuned for
#    print at ~89-183 mm width (single/double column).
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ---- Resolution / format constants (journal_config.json figure_specs) -----
FIG_DPI_LINE   <- 1200   # line art (all figures here qualify)
FIG_DPI_OTHER  <- 600    # halftone/combination minimum (kept for reference)
FIG_FORMATS    <- c("tiff")   # 연구자 확정 2026-08-12: 제출용 TIFF + 편집용 PPTX 만 생성 (PNG 중단)

# ---- Colorblind-safe palette (Okabe & Ito, 2008) ---------------------------
okabe_ito <- c(
  black      = "#000000",
  orange     = "#E69F00",
  sky_blue   = "#56B4E9",
  bluish_grn = "#009E73",
  yellow     = "#F0E442",
  blue       = "#0072B2",
  vermillion = "#D55E00",
  reddish_pu = "#CC79A7"
)

# Sex encoding — colour + shape + linetype triplet (never colour alone)
sex_levels    <- c("Male", "Female")   # 데이터 원본 값 (CSV 매칭 전용 — 표시에 쓰지 말 것)
sex_display   <- c("Men", "Women")     # 표시 라벨 (연구자 확정 2026-08-12, AMA: 사람 집단은 명사 Men/Women)
sex_factor    <- function(x) factor(x, levels = sex_levels, labels = sex_display)
sex_colors    <- c(Men = unname(okabe_ito["blue"]),  Women = unname(okabe_ito["vermillion"]))
sex_shapes    <- c(Men = 16,  Women = 17)                   # filled circle / filled triangle
sex_linetypes <- c(Men = "solid", Women = "42")              # solid / dashed
sex_fills     <- c(Men = unname(okabe_ito["blue"]),  Women = unname(okabe_ito["vermillion"]))

# Neutral / stage palette for non-sex-encoded categorical panels (staged
# adjustment, complexity strata) — also Okabe-Ito, colorblind-safe.
stage_palette <- unname(okabe_ito[c("black", "sky_blue", "blue", "vermillion", "bluish_grn", "orange")])

# ---- Base theme -------------------------------------------------------------
theme_kcto <- function(base_size = 8, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      # NOTE(Lead, 2026-08-10, 연구자 요청): 전 figure 에서 세로 보조선 제거로 통일.
      # 격자선이 사라지면 눈금 읽기가 어려워지므로, x축과 동일한 굵기·색의 y축선과
      # 눈금을 함께 부여해 두 축이 같은 무게로 보이도록 한다.
      panel.grid.major.x = element_blank(),
      axis.line.x        = element_line(color = "grey20", linewidth = 0.35),
      axis.line.y        = element_line(color = "grey20", linewidth = 0.35),
      axis.ticks.x        = element_line(color = "grey20", linewidth = 0.3),
      axis.ticks.y        = element_line(color = "grey20", linewidth = 0.3),
      axis.text           = element_text(size = rel(0.95), color = "grey20"),
      strip.background    = element_rect(fill = "grey95", color = NA),
      strip.text           = element_text(face = "bold", size = rel(0.95)),
      plot.title            = element_text(face = "bold", size = rel(1.0), hjust = 0),  # 축 타이틀과 동일 크기 (연구자 지시 2026-08-12)
      plot.subtitle         = element_text(size = rel(0.95), color = "grey30", hjust = 0),
      plot.caption           = element_text(size = rel(0.68), color = "grey35", hjust = 0, lineheight = 1.15),
      legend.position          = "bottom",
      legend.title              = element_text(face = "bold", size = rel(0.9)),
      legend.text                = element_text(size = rel(0.9)),
      panel.spacing               = unit(1.1, "lines"),
      plot.margin                  = margin(4, 8, 4, 4)
    )
}

# Minimal "flow-diagram" theme: no axes/gridlines at all (used for STROBE flow)
theme_flow <- function(base_size = 8, base_family = "Helvetica") {
  theme_void(base_size = base_size, base_family = base_family) +
    theme(
      plot.title    = element_text(face = "bold", size = rel(1.15), hjust = 0),
      plot.subtitle = element_text(size = rel(0.9), color = "grey30", hjust = 0),
      plot.caption  = element_text(size = rel(0.68), color = "grey35", hjust = 0, lineheight = 1.15),
      plot.margin   = margin(4, 4, 4, 4)
    )
}

# ---- Save helper: writes TIFF at 1200 dpi, mm-specified size ---------------
# path_noext: full path WITHOUT extension, e.g. ".../Figure1_STROBE_flow"
save_fig <- function(plot, path_noext, width_mm, height_mm, dpi = FIG_DPI_LINE,
                      formats = FIG_FORMATS, bg = "white") {
  for (fmt in formats) {
    extra_args <- if (fmt == "tiff") list(compression = "lzw") else list()
    do.call(ggplot2::ggsave, c(list(
      filename = paste0(path_noext, ".", fmt),
      plot     = plot,
      width    = width_mm,
      height   = height_mm,
      units    = "mm",
      dpi      = dpi,
      bg       = bg,
      device   = fmt
    ), extra_args))
  }
  invisible(path_noext)
}

# ---- Small formatting helpers ----------------------------------------------
fmt_hr <- function(est, lo, hi, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), est, lo, hi)
}

# Parse "est_fmt" style strings "-0.120 (-0.182 to -0.058)" is avoided; we
# always read numeric lo/hi columns directly from the CSVs (no re-derivation
# of statistics — illustrator role reads, never recomputes, estimates).

# Project root — matches the convention used by scripts/02_descriptive_by_sex.R
# and scripts/04_primary_analysis.R (hardcoded absolute PROJ constant, so that
# `Rscript scripts/07_figures.R` reproduces identically regardless of cwd).
if (!exists("PROJ")) {
  PROJ <- Sys.getenv("KCTO_PROJECT_ROOT", unset = ".")
}

TABLES_DIR  <- file.path(PROJ, "output", "tables")
FIGURES_DIR <- file.path(PROJ, "output", "figures")
if (!dir.exists(FIGURES_DIR)) dir.create(FIGURES_DIR, recursive = TRUE)

cat(sprintf("[_theme.R] PROJ = %s\n", PROJ))
cat(sprintf("[_theme.R] TABLES_DIR = %s\n", TABLES_DIR))
cat(sprintf("[_theme.R] FIGURES_DIR = %s\n", FIGURES_DIR))
