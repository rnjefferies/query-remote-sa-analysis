# =============================================================================
# normality_convention_sensitivity.R
# Chapter 3 (Study 1): sensitivity of parametric/non-parametric test selection
# to the normality convention.
#
# The reported analyses choose between a repeated-measures ANOVA and the
# non-parametric Friedman test using a PER-CONDITION Shapiro-Wilk check of the
# raw scores (see ch3_offline_vs_realtime.R and Methods Section 2.6). The
# conventional ANOVA assumption, however, concerns the normality of the model
# RESIDUALS. This script re-derives the test assignment for every repeated-
# measures variable using a residual-based Shapiro-Wilk check and reports, for
# each variable, the test and significance conclusion under BOTH conventions,
# flagging any variable whose conclusion changes.
#
# Result (n = 24): across all 17 repeated-measures variables the two conventions
# reach the same significant/non-significant conclusion for every variable
# except Temporal Demand (non-significant under the per-condition rule,
# p = .054; significant under residuals, p = .001). Temporal Demand is elevated
# in both query conditions relative to baseline with no Realtime-Offline
# difference, so it neither supports H2 nor alters any other conclusion.
#
# Self-contained: reads only the de-identified CSVs in ../data. Run from the
# repository root:  Rscript ch3_offline_vs_realtime/normality_convention_sensitivity.R
# =============================================================================

suppressMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(rstatix)
})

DATA <- if (dir.exists("data")) "data" else "../data"

# For one variable in long form (P_ID, Condition factor, Score): compute the
# per-condition raw-score Shapiro min-p and the residual Shapiro p, assign a
# test under each convention, run both tests, and compare the conclusions.
assess <- function(df, label) {
  df <- df %>% filter(!is.na(Score)) %>%
    mutate(P_ID = factor(P_ID), Condition = factor(Condition))
  keep <- df %>% count(P_ID) %>% filter(n == nlevels(df$Condition)) %>% pull(P_ID)
  df <- df %>% filter(P_ID %in% keep) %>% droplevels()
  if (nlevels(df$Condition) < 3) return(NULL)

  min_pc <- df %>% group_by(Condition) %>%
    summarise(p = tryCatch(shapiro.test(Score)$p.value, error = function(e) NA), .groups = "drop") %>%
    pull(p) %>% min(na.rm = TRUE)
  res_p <- tryCatch(shapiro.test(residuals(lm(Score ~ Condition + P_ID, data = df)))$p.value,
                    error = function(e) NA)

  fr_p <- tryCatch(friedman_test(data = df, Score ~ Condition | P_ID)$p, error = function(e) NA)
  an_p <- tryCatch(get_anova_table(anova_test(data = df, dv = Score, wid = P_ID, within = Condition))$p[1],
                   error = function(e) NA)

  test_pc <- if (!is.na(min_pc) && min_pc < .05) "Friedman" else "ANOVA"
  test_rs <- if (!is.na(res_p)  && res_p  < .05) "Friedman" else "ANOVA"
  p_pc <- if (test_pc == "Friedman") fr_p else an_p
  p_rs <- if (test_rs == "Friedman") fr_p else an_p
  concl <- function(p) ifelse(!is.na(p) && p < .05, "SIG", "ns")

  data.frame(Variable = label,
             percond_minShapiro = round(min_pc, 3), residShapiro = round(res_p, 3),
             test_percond = test_pc, p_percond = round(p_pc, 3), concl_percond = concl(p_pc),
             test_resid = test_rs, p_resid = round(p_rs, 3), concl_resid = concl(p_rs),
             FLIP = ifelse(concl(p_pc) != concl(p_rs), "*** FLIP ***", ""))
}

# ---- NASA-TLX and SART (de-identified survey table) -------------------------
SND  <- read_csv(file.path(DATA, "SND_data.csv"), show_col_types = FALSE)
labs <- SND[1, ] %>% as.character() %>% str_replace_all("\\n", " ") %>% str_trim() %>% setNames(names(SND))
clean <- SND %>% slice(-1, -2)
ntlx_labs <- labs[grepl("NTLX", names(labs))] %>% sapply(function(x) str_trim(str_split_fixed(x, " How", 2)[, 1]))
sart_labs <- labs[grepl("SART", names(labs))] %>% sapply(function(x) str_trim(str_split_fixed(x, " How", 2)[, 1]))
reshape_block <- function(prefix, kind, L) {
  clean %>% select(P_ID, all_of(paste0(prefix, "Condition")), starts_with(paste0(prefix, kind, "_"))) %>%
    rename(Condition = all_of(paste0(prefix, "Condition"))) %>%
    pivot_longer(starts_with(paste0(prefix, kind, "_")), names_to = "Metric", values_to = "Score") %>%
    mutate(Metric_Label = L[Metric], Condition = trimws(Condition), Score = as.numeric(Score))
}
ntlx <- bind_rows(lapply(c("4_", "5_", "6_"), reshape_block, kind = "NTLX", L = ntlx_labs))
sart <- bind_rows(lapply(c("4_", "5_", "6_"), reshape_block, kind = "SART", L = sart_labs))

rows <- list()
rows[["TLX Overall"]] <- assess(ntlx %>% group_by(P_ID, Condition) %>%
                                  summarise(Score = mean(Score, na.rm = TRUE), .groups = "drop"), "TLX Overall")
for (m in c("Mental Demand", "Physical Demand", "Temporal Demand", "Performance", "Effort", "Frustration"))
  rows[[m]] <- assess(ntlx %>% filter(Metric_Label == m) %>% select(P_ID, Condition, Score), paste("TLX", m))

demand <- c("Instability of Situation", "Complexity of Situation", "Variability of Situation")
supply <- c("Arousal", "Concentration of Attention", "Division of Attention", "Spare Mental Capacity")
under  <- c("Information Quantity", "Information Quality", "Familiarity with Situation")
sart_c <- sart %>%
  mutate(Cat = case_when(Metric_Label %in% demand ~ "Demand",
                         Metric_Label %in% supply ~ "Supply",
                         Metric_Label %in% under  ~ "Understanding")) %>%
  filter(!is.na(Cat)) %>%
  group_by(P_ID, Condition, Cat) %>% summarise(Total = sum(Score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Cat, values_from = Total) %>%
  mutate(SART_Score = Understanding - (Demand - Supply))
for (col in c("Demand", "Supply", "Understanding", "SART_Score"))
  rows[[paste("SART", col)]] <- assess(sart_c %>% transmute(P_ID, Condition, Score = .data[[col]]), paste("SART", col))

# ---- remote-driving performance ---------------------------------------------
cr <- function(path) read_csv(path, show_col_types = FALSE) %>%
  mutate(P_ID = gsub("participant_", "", Participant_ID), Condition = substr(Condition, 1, 1))
rows[["Fine SRR"]]  <- assess(cr(file.path(DATA, "active_steering_reversals_summary.csv"))  %>% transmute(P_ID, Condition, Score = Reversal_rate), "Fine SRR")
rows[["Major SRR"]] <- assess(cr(file.path(DATA, "active_steering_reversals_summary2.csv")) %>% transmute(P_ID, Condition, Score = Reversal_rate), "Major SRR")
rows[["Throttle"]]  <- assess(cr(file.path(DATA, "active_throttle_reversals_summary3.csv")) %>% transmute(P_ID, Condition, Score = Reversal_rate), "Throttle reversals")

sp <- read_csv(file.path(DATA, "speed_metrics_active_driving_FIXED.csv"), show_col_types = FALSE) %>%
  rename(P_ID = "Participant ID", AV = "Average Speed (mph)", SSD = "Speed StdDev (mph)", MX = "Max Speed (mph)") %>%
  mutate(Condition = substr(Condition, 1, 1))
grp <- tibble(P_ID = unique(sp$P_ID)) %>% mutate(k = ifelse(row_number() <= 12, 7.442481, 5.446379))
sp <- sp %>% left_join(grp, by = "P_ID")
rows[["Norm mean speed"]] <- assess(sp %>% transmute(P_ID, Condition, Score = AV / k),  "Normalised mean speed")
rows[["Norm max speed"]]  <- assess(sp %>% transmute(P_ID, Condition, Score = MX / k),  "Normalised max speed")
rows[["Norm speed SD"]]   <- assess(sp %>% transmute(P_ID, Condition, Score = SSD / k), "Normalised speed SD")

out <- bind_rows(rows)
options(width = 200)
cat("\n================ Test selection under per-condition vs residual-based normality ================\n")
print(out, row.names = FALSE)
flips <- out %>% filter(FLIP != "")
cat("\nConclusion changes under residual-based normality:\n")
if (nrow(flips) == 0) cat("  NONE.\n") else
  print(flips %>% select(Variable, concl_percond, p_percond, concl_resid, p_resid), row.names = FALSE)
