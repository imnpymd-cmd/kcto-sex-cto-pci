# =============================================================================
# 11_figures_pptx.R — 연구자 직접 편집용 editable-vector PPTX 내보내기
#
# 목적: Figure 1–3 + Central Illustration(채택본 v3)을 PowerPoint 개체(도형·선·
#       텍스트)로 직접 수정할 수 있는 .pptx 로 저장한다 (rvg::dml — 비트맵 아님).
# 출력: output/figures/editable_pptx/{Figure1..3, CentralIllustration}.pptx
#       각 파일 1슬라이드, 그림이 슬라이드 좌상단 (0,0)에 실제 제출 치수로 배치.
#       슬라이드 크기 자체를 그림 치수에 맞추는 보정은 후처리(셸)에서 수행한다
#       (officer 는 슬라이드 크기 변경 API 가 없음).
#
# 주의: 07/10 스크립트를 소싱하면 제출용 PNG/TIFF 가 동일 내용으로 재저장된다.
#       바이트 동일성은 셸 래퍼에서 해시로 확인하고, 달라지면 백업으로 복원한다.
#       (figure-craft 규약: 산출물 해시가 의도 없이 움직이면 실결함)
# =============================================================================

suppressMessages({ library(officer); library(rvg) })

source("scripts/07_figures.R")                    # p1, p2, p3 정의 (+제출본 재저장)
source("scripts/10_central_illustration_jacc.R")  # v3, frame() 정의 (채택본만)

OUT_DIR <- "output/figures/editable_pptx"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

mm2in <- function(mm) mm / 25.4

write_fig_pptx <- function(plot, filename, width_mm, height_mm) {
  path <- file.path(OUT_DIR, filename)
  doc <- read_pptx()
  doc <- add_slide(doc, layout = "Blank", master = "Office Theme")
  doc <- ph_with(doc, dml(ggobj = plot, bg = "white"),
                 location = ph_location(left = 0, top = 0,
                                        width  = mm2in(width_mm),
                                        height = mm2in(height_mm)))
  print(doc, target = path)
  cat(sprintf("  %s (%d x %d mm)\n", filename, width_mm, height_mm))
  invisible(path)
}

# 제출 치수 그대로 (scripts/07_figures.R · 10_central_illustration_jacc.R 의 save_fig 인자)
write_fig_pptx(p1,        "Figure1_STROBE_flow.pptx",         183, 108)
write_fig_pptx(p2,        "Figure2_age_attenuation.pptx",     140, 85)
write_fig_pptx(p3,        "Figure3_process_of_care_gap.pptx", 140, 60)
write_fig_pptx(frame(v3), "CentralIllustration.pptx",         210, 210)

cat("editable PPTX 4종 저장 완료\n")
