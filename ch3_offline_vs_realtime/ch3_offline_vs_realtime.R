# =============================================================================
# ch3_offline_vs_realtime.R
# Chapter 3 (Study 1): Offline versus real-time SA measurement under three
# SA-query methods (No Query = A, Realtime = B, Offline = C), within-subjects,
# 24 operators.
#
# Reproduces the reported analyses: NASA-TLX workload, subjective SA (SART),
# objective SA, and remote-driving performance metrics.
#
# Self-contained: reads only the de-identified CSVs in ../data. Run from the
# repository root:  Rscript ch3_offline_vs_realtime/ch3_offline_vs_realtime.R
# =============================================================================

suppressMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(rstatix)
})

DATA <- if (dir.exists("data")) "data" else "../data"
report <- function(label, value) cat(sprintf("  %-30s %s\n", label, value))

# -----------------------------------------------------------------------------
# 1. Load the de-identified survey table and reshape to long (NTLX + SART)
# -----------------------------------------------------------------------------
SND_data <- read_csv(file.path(DATA, "SND_data.csv"), show_col_types = FALSE)
metric_labels   <- SND_data[1, ]                 # row 1 holds the Qualtrics item wording
SND_data_clean  <- SND_data %>% slice(-1, -2)    # drop the two Qualtrics metadata rows
clean_labels <- metric_labels %>% as.character() %>% str_replace_all("\\n", " ") %>%
  str_trim() %>% setNames(names(SND_data))
ntlx_labels <- clean_labels[grepl("NTLX", names(clean_labels))] %>%
  sapply(function(x) str_trim(str_split_fixed(x, " How", 2)[, 1]), USE.NAMES = TRUE)
sart_labels <- clean_labels[grepl("SART", names(clean_labels))] %>%
  sapply(function(x) str_trim(str_split_fixed(x, " How", 2)[, 1]), USE.NAMES = TRUE)

reshape_block <- function(prefix, condition_col, kind, labs) {
  SND_data_clean %>%
    select(P_ID, all_of(condition_col), starts_with(paste0(prefix, kind, "_"))) %>%
    rename(Condition = all_of(condition_col)) %>%
    pivot_longer(starts_with(paste0(prefix, kind, "_")), names_to = "Metric", values_to = "Score") %>%
    mutate(Metric_Label = labs[Metric], Condition = trimws(Condition))
}
ntlx <- bind_rows(lapply(c("4_", "5_", "6_"), function(p)
  reshape_block(p, paste0(p, "Condition"), "NTLX", ntlx_labels))) %>%
  mutate(P_ID = factor(P_ID), Condition = factor(Condition, levels = c("A", "B", "C")),
         Score = as.numeric(Score))
sart <- bind_rows(lapply(c("4_", "5_", "6_"), function(p)
  reshape_block(p, paste0(p, "Condition"), "SART", sart_labels))) %>%
  mutate(P_ID = factor(P_ID), Condition = factor(Condition, levels = c("A", "B", "C")),
         Score = as.numeric(Score))

# -----------------------------------------------------------------------------
# 2. Workload (NASA-TLX)
# -----------------------------------------------------------------------------
cat("\n== 2. Workload (NASA-TLX) ==\n")
overall <- ntlx %>% group_by(P_ID, Condition) %>%
  summarise(Overall = mean(Score, na.rm = TRUE), .groups = "drop")
ow <- get_anova_table(anova_test(data = overall, dv = Overall, wid = P_ID, within = Condition))
report("Overall workload ANOVA", sprintf("F=%.2f, ges=%.3f", ow$F[1], ow$ges[1]))
overall %>% group_by(Condition) %>% summarise(M = mean(Overall), SD = sd(Overall)) %>% print()

for (m in c("Mental Demand", "Frustration", "Effort", "Physical Demand", "Temporal Demand")) {
  d <- ntlx %>% filter(Metric_Label == m)
  f <- friedman_test(data = d, Score ~ Condition | P_ID)
  report(paste(m, "Friedman"), sprintf("chi2=%.2f, p=%.3f", f$statistic, f$p))
}
perf_sub <- ntlx %>% filter(Metric_Label == "Performance")
pa <- get_anova_table(anova_test(data = perf_sub, dv = Score, wid = P_ID, within = Condition))
report("Performance subscale ANOVA", sprintf("F=%.2f, p=%.3f", pa$F[1], pa$p[1]))

# -----------------------------------------------------------------------------
# 3. Subjective SA (SART): dimensions and total
# -----------------------------------------------------------------------------
cat("\n== 3. Subjective SA (SART) ==\n")
demand_items <- c("Instability of Situation", "Complexity of Situation", "Variability of Situation")
supply_items <- c("Arousal", "Concentration of Attention", "Division of Attention", "Spare Mental Capacity")
under_items  <- c("Information Quantity", "Information Quality", "Familiarity with Situation")
sart_comp <- sart %>%
  mutate(Cat = case_when(Metric_Label %in% demand_items ~ "Demand",
                         Metric_Label %in% supply_items ~ "Supply",
                         Metric_Label %in% under_items  ~ "Understanding")) %>%
  group_by(P_ID, Condition, Cat) %>% summarise(Total = sum(Score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Cat, values_from = Total) %>%
  mutate(SART_Score = Understanding - (Demand - Supply))   # classic SART composite
for (col in c("Demand", "Supply", "Understanding", "SART_Score")) {
  a <- get_anova_table(anova_test(data = sart_comp, dv = !!sym(col), wid = P_ID, within = Condition))
  report(paste("SART", col, "ANOVA"), sprintf("F(%.0f,%.0f)=%.2f, p=%.3f", a$DFn[1], a$DFd[1], a$F[1], a$p[1]))
}

# -----------------------------------------------------------------------------
# 4. Objective SA: percent-correct, Realtime (B) vs Offline (C)
# -----------------------------------------------------------------------------
cat("\n== 4. Objective SA ==\n")
sa <- read_csv(file.path(DATA, "clean_flagged_events.csv"), show_col_types = FALSE) %>%
  rename(P_ID = "Participant ID", Answer_Accuracy = "Answer Accuracy") %>%
  mutate(Condition = substr(Condition, 1, 1)) %>%
  filter(!P_ID %in% c("99_OL", "80_RT", "90_RT"))
sa_scores <- sa %>% group_by(P_ID, Condition) %>%
  summarise(SA_Score = mean(Answer_Accuracy == "Correct") * 100, .groups = "drop")
wide <- sa_scores %>% pivot_wider(names_from = Condition, values_from = SA_Score)
tt <- t.test(wide$B, wide$C, paired = TRUE)
report("Objective SA paired t-test", sprintf("t(%.0f)=%.2f, diff=%.2f", tt$parameter, tt$statistic, tt$estimate))

# -----------------------------------------------------------------------------
# 5. Remote-driving performance
# -----------------------------------------------------------------------------
cat("\n== 5. Performance ==\n")
clean_rev <- function(path) read_csv(path, show_col_types = FALSE) %>%
  mutate(P_ID = gsub("participant_", "", Participant_ID), Condition = substr(Condition, 1, 1),
         P_ID = factor(P_ID), Condition = factor(Condition, levels = c("A", "B", "C")))
frd <- function(d, rate, lab) { f <- friedman_test(data = d, as.formula(paste0("`", rate, "` ~ Condition | P_ID")))
  report(lab, sprintf("chi2=%.2f, p=%.3f", f$statistic, f$p)) }
anv <- function(d, dv, lab) { a <- get_anova_table(anova_test(data = d, dv = !!sym(dv), wid = P_ID, within = Condition))
  report(lab, sprintf("F(%.2f,%.2f)=%.2f, p=%.3f", a$DFn[1], a$DFd[1], a$F[1], a$p[1])) }

fine  <- clean_rev(file.path(DATA, "active_steering_reversals_summary.csv"))
major <- clean_rev(file.path(DATA, "active_steering_reversals_summary2.csv"))
thr   <- clean_rev(file.path(DATA, "active_throttle_reversals_summary3.csv"))
frd(fine,  "Reversal_rate", "Fine steering-reversal rate")
frd(major, "Reversal_rate", "Major steering-reversal rate")
anv(thr,   "Reversal_rate", "Throttle-reversal rate")

# session minutes (from the throttle summary) used to normalise event counts
sess <- thr %>% mutate(mins = as.numeric(sub(":.*", "", Session_Length_min_sec)) +
                              as.numeric(sub(".*:", "", Session_Length_min_sec)) / 60) %>%
  select(P_ID, Condition, mins)

# Speed: normalised by counterbalance-group max speed (first 12 / second 12 operators)
sp <- read_csv(file.path(DATA, "speed_metrics_active_driving_FIXED.csv"), show_col_types = FALSE) %>%
  rename(P_ID = "Participant ID", AV = "Average Speed (mph)", SD = "Speed StdDev (mph)", MX = "Max Speed (mph)") %>%
  mutate(Condition = substr(Condition, 1, 1))
grp <- tibble(P_ID = unique(sp$P_ID)) %>% mutate(k = ifelse(row_number() <= 12, 7.442481, 5.446379))
sp <- sp %>% left_join(grp, by = "P_ID") %>%
  mutate(NAV = AV / k, NSD = SD / k, NMX = MX / k,
         P_ID = factor(P_ID), Condition = factor(Condition, levels = c("A", "B", "C")))
anv(sp, "NAV", "Normalised mean speed")
anv(sp, "NMX", "Normalised max speed")
anv(sp, "NSD", "Normalised speed SD")

# Collision rate (per minute)
col <- read_csv(file.path(DATA, "collisions_clean_full.csv"), show_col_types = FALSE) %>%
  mutate(Condition = substr(Condition, 1, 1)) %>%
  group_by(P_ID = Participant_ID, Condition) %>% summarise(nc = n(), .groups = "drop")
colr <- sess %>% mutate(P_ID = as.character(P_ID)) %>%
  left_join(col, by = c("P_ID", "Condition")) %>%
  mutate(nc = replace_na(nc, 0), rate = nc / mins,
         P_ID = factor(P_ID), Condition = factor(Condition, levels = c("A", "B", "C")))
frd(colr, "rate", "Collision rate (per min)")

# Braking rate (per minute) from the complete preprocessed file
brk <- read_csv(file.path(DATA, "hard_braking_summary_preprocessed.csv"), show_col_types = FALSE) %>%
  rename(P_ID = "Participant ID", ev = "Hard Braking Events", sl = "Session Length (min:sec)") %>%
  mutate(Condition = substr(Condition, 1, 1),
         mins = as.numeric(sub(":.*", "", sl)) + as.numeric(sub(".*:", "", sl)) / 60,
         rate = ev / mins, P_ID = factor(P_ID), Condition = factor(Condition, levels = c("A", "B", "C")))
frd(brk, "rate", "Braking rate (per min)")

# -----------------------------------------------------------------------------
# 6. Demographic summary (aggregate counts)
# -----------------------------------------------------------------------------
cat("\n== 6. Demographics (aggregate; per-participant demographics not distributed) ==\n")
demo_summary <- read_csv(file.path(DATA, "demographics_summary.csv"), show_col_types = FALSE)
print(as.data.frame(demo_summary))

cat("\nDone.\n")
