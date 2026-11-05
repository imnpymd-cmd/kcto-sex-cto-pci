# =============================================================================
# 03_derive_endpoints.R
# 출판된 parent paper (Korean Circ J 2027;57(1):e4, doi:10.4070/kcj.2025.0412)
# 의 복합 endpoint 정의를 그대로 재현하여 파생한다.
#
# 배경: 데이터셋의 기존 `MACE` / `NACE` 컬럼은 출판본 정의와 **다르다**
#       (기존 컬럼은 revascularization 중 TVR만 포함). 그대로 쓰면 안 된다.
#
# 출판본 Methods 원문 정의:
#   MACCE (primary) = all-cause death + MI + CVA + any form of coronary revascularization
#   MACE            = all-cause death + MI + any revascularization
#   NACE            = MACCE + major bleeding (TIMI criteria)
#
# 검증: 전체 9,316명 5년 KM 누적발생률
#   MACE  28.5% (출판본 28.5%) / MACCE 30.5% (30.4%)
#   NACE  31.0% (30.9%)        / Death 14.6% (14.5%)
#   → 반올림 오차 범위 내 완전 재현
#
# 사용법: source("scripts/03_derive_endpoints.R"); d <- derive_endpoints(d)
# =============================================================================

suppressMessages(library(survival))

.num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# 복합 endpoint 생성기 -----------------------------------------------------
# 이벤트 발생 = 구성요소 중 하나라도 발생
# 이벤트 시각 = 발생한 구성요소들의 최소 시각
# 무사건 시각 = 구성요소 추적기간 중 최대값 (가장 긴 추적 정보를 보존)
.make_composite <- function(df, ev_vars, tm_vars) {
  stopifnot(length(ev_vars) == length(tm_vars))
  E <- sapply(ev_vars, function(v) ifelse(is.na(df[[v]]), 0L, as.integer(df[[v]] == 1)))
  T <- sapply(tm_vars, function(v) df[[v]])
  any_ev  <- as.integer(rowSums(E) > 0)
  t_event <- apply(ifelse(E == 1, T, Inf), 1, min)
  t_cens  <- apply(T, 1, function(r) suppressWarnings(max(r, na.rm = TRUE)))
  list(event = any_ev, time = ifelse(any_ev == 1, t_event, t_cens))
}

derive_endpoints <- function(df) {

  ev <- c("Death", "MI", "CVA_re", "TVR_re2", "NLR_re2", "CABG_re", "Bleeding_major")
  tm <- c("death_fudur", "MI_fudur", "CVA_fudur", "TVR_fudur2", "NLR_fudur2",
          "CABG_fudur", "Bleeding_major_fudur")
  for (v in c(ev, tm)) df[[v]] <- .num(df[[v]])

  # --- CVA_re 값 2 (2건) 처리 -----------------------------------------------
  # codebook에 정의되지 않은 값 2가 2건 존재. 연구자 확인 결과 **0(No)의 오기**로
  # 판정됨 (2026-08-08). 명시적으로 0으로 재코딩한다.
  n_recoded <- sum(df$CVA_re == 2, na.rm = TRUE)
  df$CVA_re[!is.na(df$CVA_re) & df$CVA_re == 2] <- 0
  if (n_recoded > 0) message(sprintf("CVA_re: 값 2 %d건을 0(No)으로 재코딩 (연구자 확인)", n_recoded))
  stopifnot(all(df$CVA_re %in% c(0, 1) | is.na(df$CVA_re)))

  # any revascularization = TVR + NLR + CABG
  revasc_ev <- c("TVR_re2", "NLR_re2", "CABG_re")
  revasc_tm <- c("TVR_fudur2", "NLR_fudur2", "CABG_fudur")

  # --- MACE (출판본 정의) ---
  o <- .make_composite(df, c("Death", "MI", revasc_ev),
                           c("death_fudur", "MI_fudur", revasc_tm))
  df$MACE_pub <- o$event; df$MACE_pub_fudur <- o$time

  # --- MACCE (출판본 primary endpoint) ---
  o <- .make_composite(df, c("Death", "MI", "CVA_re", revasc_ev),
                           c("death_fudur", "MI_fudur", "CVA_fudur", revasc_tm))
  df$MACCE_pub <- o$event; df$MACCE_pub_fudur <- o$time

  # --- NACE (출판본 정의 = MACCE + major bleeding) ---
  o <- .make_composite(df, c("Death", "MI", "CVA_re", revasc_ev, "Bleeding_major"),
                           c("death_fudur", "MI_fudur", "CVA_fudur", revasc_tm,
                             "Bleeding_major_fudur"))
  df$NACE_pub <- o$event; df$NACE_pub_fudur <- o$time

  # --- any revascularization 단독 (구성요소 분해 보고용) ---
  o <- .make_composite(df, revasc_ev, revasc_tm)
  df$AnyRevasc <- o$event; df$AnyRevasc_fudur <- o$time

  # === PRIMARY OUTCOME (본 연구, 연구자 확정 2026-08-08) ====================
  # DMS = all-cause Death + MI + Stroke(CVA). **재관류 제외.**
  # 근거 (CLAUDE.md §1-2 "논거 배치 순서" 참조):
  #   1. 순환논리 회피 — 본 연구는 재관류 관련 치료 결정의 성별 차이를 조사 대상으로
  #      삼는다. 조사 대상인 노출 경로를 결과변수에 포함하면 치료 선택과 환자 예후가
  #      한 지표에서 뒤섞인다. (결과와 무관하게 성립하는 설계 논거)
  #   2. 문헌 정합성 — Gallo 2026 (PMID 41753140) 이 동일 구성을 "MACE"로 사용.
  # ⛔ "사망과 재관류가 정반대 방향" 논거는 폐기됨(2026-08-08 QA). 완전보정 후
  #    death 0.93 / AnyRevasc 0.88 로 같은 방향이며, 방향 충돌은 연령이 만든
  #    미보정 현상이다. 이 문구를 되살리지 말 것.
  o <- .make_composite(df, c("Death", "MI", "CVA_re"),
                           c("death_fudur", "MI_fudur", "CVA_fudur"))
  df$DMS <- o$event; df$DMS_fudur <- o$time

  df
}

# 검증 함수: 전체 코호트에서 출판본 수치를 재현하는지 확인 --------------------
verify_against_published <- function(df) {
  chk <- function(e, t, lab, ref) {
    f <- survfit(Surv(df[[t]] / 365.25, df[[e]]) ~ 1)
    s <- summary(f, times = 5, extend = TRUE)
    cat(sprintf("  %-6s 5-yr = %5.1f%%   [published %s]\n", lab, 100 * (1 - s$surv), ref))
  }
  cat("전체 코호트 N =", nrow(df), " — 출판본 대조\n")
  chk("MACE_pub",  "MACE_pub_fudur",  "MACE",  "28.5%")
  chk("MACCE_pub", "MACCE_pub_fudur", "MACCE", "30.4%")
  chk("NACE_pub",  "NACE_pub_fudur",  "NACE",  "30.9%")
  chk("Death",     "death_fudur",     "Death", "14.5%")
  invisible(NULL)
}
