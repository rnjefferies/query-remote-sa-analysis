#### confidence_pipeline.R ##################################################
# Reproducible, CSV-in / CSV-out pipeline for the CONFIDENCE / OVERCONFIDENCE
# analyses of Chapter 4 (Study 1).
#
# Rebuilds, from the raw CSVs alone, every confidence result reported in the
# chapter: the calibration and resolution indices, the response-category
# decomposition, the latency contrasts, the operator-level driving-style
# correlations, and the real-time vs offline type-2 AUROC. Every intermediate
# is written to outputs/ so the reported numbers can be traced.
#
# INPUTS  : clean_flagged_events.csv, collisions_clean_full.csv,
#           speed_metrics_active_driving_FIXED.csv,
#           active_steering_reversals_summary.csv   (fine SRR, ~4.5 deg gap),
#           active_steering_reversals_summary2.csv  (major SRR, ~45 deg gap),
#           active_throttle_reversals_summary.csv
# OUTPUTS : outputs/confidence_query_level.csv
#           outputs/confidence_calibration_summary.csv       (Sec 3.1/3.2, Table 1)
#           outputs/confidence_calibration_indices.csv       (Sec 3.2)
#           outputs/confidence_calibration_indices_ci.csv    (Sec 3.2, bootstrap CI)
#           outputs/confidence_category_summary.csv          (Sec 3.4, Table 2)
#           outputs/confidence_latency_separation.csv        (Sec 3.4/3.5, median diff + prob. superiority)
#           outputs/confidence_within_level_latency.csv      (Sec 3.5, confident-error slowing per SA level)
#           outputs/confidence_driving_correlations.csv      (Sec 3.6, Spearman + bootstrap CI)
#           outputs/confidence_type2_auroc.csv               (Sec 3.8, AUROC2 real-time vs offline)
#           outputs/confidence_operator_level.csv            (per-operator telemetry + query measures)
#           outputs/confidence_per_operator.csv              (Sec 3.6/3.9 profiles)
#           outputs/confidence_threshold_sweepA_overconf.csv, _sweepB_cautious.csv (Sec 3.7)
#           outputs/confidence_numbers.txt                   (headline numbers)
#
# Run:  Rscript confidence_pipeline.R
#############################################################################

suppressMessages({
    library(readr)
    library(dplyr)
    library(stringr)
    library(purrr)
    library(tidyr)
})

set.seed(123)
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE)
log_lines <- c()
say <- function(...) { line <- paste0(...); cat(line, "\n"); log_lines[[length(log_lines) + 1]] <<- line }

#### 1. Load and clean SA query events ######################################
events <- read_csv("data/clean_flagged_events.csv", show_col_types = FALSE) %>%
    rename(
        P_ID                 = "Participant ID",
        Event_ID             = "Event ID",
        Answer_Accuracy      = "Answer Accuracy",
        Time_Difference_s    = "Time Difference (s)",
        Question_Timestamp_s = "Question Timestamp (s)",
        Answer_Timestamp_s   = "Answer Timestamp (s)"
    ) %>%
    filter(!P_ID %in% c("99_OL", "80_RT", "90_RT"))   # drop dummy IDs

# Realtime condition only (B_L / B_R) -> the queries answered while driving
Real_Query_data <- events %>% filter(Condition %in% c("B_L", "B_R"))

#### 2. Attach collision-proximity flag (+/- 5s of the query) ###############
collisions <- read_csv("data/collisions_clean_full.csv", show_col_types = FALSE) %>%
    rename(P_ID = "Participant_ID") %>%
    filter(Condition %in% c("B_L", "B_R"))

check_collision <- function(p_id, condition, q_time, collision_df) {
    ct <- collision_df %>%
        filter(P_ID == p_id, Condition == condition) %>%
        pull(Collision_Timestamp)
    if (length(ct) == 0) return(FALSE)
    any(abs(ct - q_time) <= 5)
}

Real_Query_data$Collision_PrePost <- pmap_lgl(
    list(Real_Query_data$P_ID, Real_Query_data$Condition, Real_Query_data$Question_Timestamp_s),
    ~ check_collision(..1, ..2, ..3, collisions)
) %>% as.integer()

#### 3. Recode accuracy to 0/1 ##############################################
Real_Query_data <- Real_Query_data %>%
    mutate(Answer_Accuracy = ifelse(Answer_Accuracy == "Correct", 1L, 0L))

#### 4. Four behavioural categories (confidence >= 6 = "confident") ##########
# Confident = high-confidence top-two-box (>=6), matching the overconfidence
# cut and the SA-state indicator. Cautious = all other responses (<=5), the
# complement, so no data are set aside.
Real_Query_data <- Real_Query_data %>%
    mutate(Query_Type = case_when(
        Answer_Accuracy == 1 & Confidence >= 6 ~ "Confident Correct",
        Answer_Accuracy == 0 & Confidence >= 6 ~ "Confident Incorrect",
        Answer_Accuracy == 1 & Confidence <= 5 ~ "Cautious Correct",
        Answer_Accuracy == 0 & Confidence <= 5 ~ "Cautious Incorrect"
    ))
low_correct_pct <- 100 * mean(Real_Query_data$Answer_Accuracy[Real_Query_data$Confidence <= 5] == 1)

#### 5. Headline calibration numbers ########################################
n_q      <- nrow(Real_Query_data)
n_op     <- n_distinct(Real_Query_data$P_ID)
acc      <- mean(Real_Query_data$Answer_Accuracy)
say("==== CONFIDENCE / OVERCONFIDENCE PIPELINE ====")
say(sprintf("Realtime queries: %d across %d operators", n_q, n_op))
say(sprintf("Overall accuracy: %.1f%%", 100 * acc))
say(sprintf("Point-biserial r(confidence, accuracy): %.3f",
            cor(Real_Query_data$Confidence, Real_Query_data$Answer_Accuracy)))
say(sprintf("Pearson r(confidence, latency): %.3f",
            cor(Real_Query_data$Confidence, Real_Query_data$Time_Difference_s)))
# Within-operator mean confidence, correct vs incorrect answers (reported in Sec 3.2)
wo_conf <- Real_Query_data %>% group_by(P_ID) %>%
    summarise(c_correct = mean(Confidence[Answer_Accuracy == 1], na.rm = TRUE),
              c_wrong   = mean(Confidence[Answer_Accuracy == 0], na.rm = TRUE), .groups = "drop")
say(sprintf("Within-operator mean confidence: correct=%.2f, incorrect=%.2f",
            mean(wo_conf$c_correct, na.rm = TRUE), mean(wo_conf$c_wrong, na.rm = TRUE)))

# Overconfidence: high confidence (>=6) yet wrong
high_conf <- Real_Query_data %>% filter(Confidence >= 6)
n_err_tot <- sum(Real_Query_data$Answer_Accuracy == 0)
say(sprintf("High-confidence (>=6) queries: %d (%.1f%% of all)",
            nrow(high_conf), 100 * nrow(high_conf) / n_q))
say(sprintf("  of which INCORRECT (overconfident): %d (%.1f%% of high-conf)",
            sum(high_conf$Answer_Accuracy == 0),
            100 * mean(high_conf$Answer_Accuracy == 0)))
say(sprintf("Overconfident errors as %% of ALL errors: %.1f%%",
            100 * sum(high_conf$Answer_Accuracy == 0) / n_err_tot))

#### 5b. Accuracy and latency by confidence level (Table 1) ##################
# Table 1 reports accuracy and MEDIAN latency per level; the same grouping also
# supplies the p_k / o_k bins for the calibration decomposition below.
Real_Query_data$conf_p <- (Real_Query_data$Confidence - 1) / 6   # map 1-7 scale to p=(c-1)/6
obar <- mean(Real_Query_data$Answer_Accuracy)
Nq_c <- nrow(Real_Query_data)
by_conf <- Real_Query_data %>%
    group_by(Confidence) %>%
    summarise(n = n(), accuracy = mean(Answer_Accuracy),
              median_latency_s = median(Time_Difference_s, na.rm = TRUE),
              mean_latency_s   = mean(Time_Difference_s, na.rm = TRUE),
              p_k = unique((Confidence - 1) / 6),
              o_k = mean(Answer_Accuracy), .groups = "drop")
calib <- by_conf %>% select(Confidence, n, accuracy, median_latency_s, mean_latency_s)

#### 5c. Calibration / meta-SA indices (Lichacz 2008; Murphy 1973) ###########
# Resolution is mapping-independent; reliability and bias depend on p=(c-1)/6.
reliability <- sum(by_conf$n * (by_conf$p_k - by_conf$o_k)^2) / Nq_c   # lower = better
resolution  <- sum(by_conf$n * (by_conf$o_k - obar)^2) / Nq_c          # higher = better
uncertainty <- obar * (1 - obar)
norm_resolution <- resolution / uncertainty                           # 0..1 (Murphy)
bias <- mean(Real_Query_data$conf_p) - obar                           # +over / -under
brier <- mean((Real_Query_data$conf_p - Real_Query_data$Answer_Accuracy)^2)
cal_indices <- data.frame(
    index = c("over_under_confidence_bias", "calibration_reliability",
              "resolution", "normalised_resolution", "uncertainty", "brier"),
    value = round(c(bias, reliability, resolution, norm_resolution, uncertainty, brier), 4))
say("")
say("---- Calibration / meta-SA indices (conf mapped p=(c-1)/6) ----")
say(sprintf("over/under-confidence bias = %+.3f  (%sconfident overall)",
            bias, ifelse(bias > 0, "over", "under")))
say(sprintf("calibration (reliability)  = %.4f  (lower = better)", reliability))
say(sprintf("resolution                 = %.4f ; normalised = %.3f (0-1, higher = better)",
            resolution, norm_resolution))
say(sprintf("Brier score                = %.4f", brier))

#### 6. Category profiles (Table 2) #########################################
cat_summary <- Real_Query_data %>%
    filter(!is.na(Query_Type)) %>%
    group_by(Query_Type) %>%
    summarise(
        n              = n(),
        pct            = 100 * n() / n_q,
        mean_conf      = mean(Confidence, na.rm = TRUE),
        mean_latency_s   = mean(Time_Difference_s, na.rm = TRUE),
        median_latency_s = median(Time_Difference_s, na.rm = TRUE),  # reported in Table 2
        collision_rate = mean(Collision_PrePost, na.rm = TRUE),
        .groups = "drop"
    )
say("")
say("Cut-points: confident >=6 (top-two-box), cautious <=5 (complement, no exclusion).")
say(sprintf("Cautious (<=5) responses correct: %.1f%%", low_correct_pct))
say("---- Behavioural categories (latency: raw seconds) ----")
for (i in seq_len(nrow(cat_summary))) {
    r <- cat_summary[i, ]
    say(sprintf("%-20s n=%3d (%4.1f%%)  conf=%.2f  median_latency=%.2fs  collide=%.1f%%",
                r$Query_Type, r$n, r$pct, r$mean_conf, r$median_latency_s,
                100 * r$collision_rate))
}

#### 7. Threshold sensitivity (robustness, Section 3.7) #####################
# Confirm the findings are not driven by the exact confidence cut-points.
# Sweep A: overconfidence (high) cut in {5,6,7} -> share of errors made at high
# confidence. Sweep B: cautious (low) cut in {3,4,5} -> size and accuracy of the
# cautious group.
sweepA <- do.call(rbind, lapply(c(5, 6, 7), function(hh) {
    hi <- Real_Query_data$Confidence >= hh
    oc <- hi & Real_Query_data$Answer_Accuracy == 0
    data.frame(high_cut = hh, n_high = sum(hi), overconf_errors = sum(oc),
               pct_all_errors = round(100 * sum(oc) / n_err_tot, 1),
               pct_of_high = round(100 * sum(oc) / sum(hi), 2))
}))
sweepB <- do.call(rbind, lapply(c(3, 4, 5), function(ll) {
    ca <- Real_Query_data$Confidence <= ll
    data.frame(low_cut = ll, n_cautious = sum(ca),
               cautious_correct_pct = round(100 * mean(Real_Query_data$Answer_Accuracy[ca] == 1), 1))
}))
say("")
say("---- Threshold sensitivity (Sweep A: overconfidence cut) ----")
say(paste(capture.output(print(sweepA, row.names = FALSE)), collapse = "\n"))
say("---- Threshold sensitivity (Sweep B: cautious cut) ----")
say(paste(capture.output(print(sweepB, row.names = FALSE)), collapse = "\n"))

#### 8. Cluster-bootstrap machinery (resample whole operators) ##############
# All CIs below resample operators with replacement (queries within an operator
# are not independent), 5,000 draws, matching the design of Sections 3.4/3.8.
B_REP    <- 5000
ids_all  <- unique(Real_Query_data$P_ID)
split_rt <- split(Real_Query_data, Real_Query_data$P_ID)

# Rank-based probability of superiority: P(a > b) with ties counted as 0.5.
prob_sup <- function(a, b) {
    na <- length(a); nb <- length(b)
    if (na == 0 || nb == 0) return(NA_real_)
    r <- rank(c(a, b))
    (sum(r[seq_len(na)]) - na * (na + 1) / 2) / (na * nb)
}

#### 8a. Calibration indices with cluster-bootstrap 95% CIs (Section 3.2) ####
cal_index_fun <- function(df) {
    obar <- mean(df$Answer_Accuracy); N <- nrow(df)
    bins <- df %>% group_by(Confidence) %>%
        summarise(n = n(), p_k = unique((Confidence - 1) / 6),
                  o_k = mean(Answer_Accuracy), .groups = "drop")
    unc <- obar * (1 - obar)
    c(bias                  = mean((df$Confidence - 1) / 6) - obar,
      reliability           = sum(bins$n * (bins$p_k - bins$o_k)^2) / N,
      normalised_resolution = (sum(bins$n * (bins$o_k - obar)^2) / N) / unc)
}
set.seed(2024)
cal_reps <- replicate(B_REP, cal_index_fun(bind_rows(split_rt[sample(ids_all, replace = TRUE)])))
cal_ci   <- t(apply(cal_reps, 1, quantile, probs = c(.025, .975), na.rm = TRUE))
cal_indices_ci <- data.frame(
    index = c("over_under_confidence_bias", "calibration_reliability", "normalised_resolution"),
    value = round(c(bias, reliability, norm_resolution), 4),
    ci_lo = round(cal_ci[, 1], 4), ci_hi = round(cal_ci[, 2], 4))
say(""); say("---- Calibration indices with 95% cluster-bootstrap CI ----")
say(paste(capture.output(print(cal_indices_ci, row.names = FALSE)), collapse = "\n"))

#### 8b. Latency separation: median diff + probability of superiority ########
lat_by <- function(df, cat) df$Time_Difference_s[df$Query_Type == cat & !is.na(df$Query_Type)]
boot_stat <- function(fun, B = B_REP) replicate(B, {
    dd <- bind_rows(split_rt[sample(ids_all, replace = TRUE)]); fun(dd) })
contrast_row <- function(c1, c2) {
    md_pt  <- median(lat_by(Real_Query_data, c1)) - median(lat_by(Real_Query_data, c2))
    ps_pt  <- prob_sup(lat_by(Real_Query_data, c1), lat_by(Real_Query_data, c2))
    set.seed(2024); md <- boot_stat(function(dd) median(lat_by(dd, c1)) - median(lat_by(dd, c2)))
    set.seed(2024); ps <- boot_stat(function(dd) prob_sup(lat_by(dd, c1), lat_by(dd, c2)))
    data.frame(contrast = paste(c1, "vs", c2),
               median_diff_s = round(md_pt, 2),
               md_lo = round(quantile(md, .025, na.rm = TRUE), 2),
               md_hi = round(quantile(md, .975, na.rm = TRUE), 2),
               prob_superiority = round(ps_pt, 2),
               pos_lo = round(quantile(ps, .025, na.rm = TRUE), 2),
               pos_hi = round(quantile(ps, .975, na.rm = TRUE), 2), row.names = NULL)
}
latency_separation <- rbind(
    contrast_row("Confident Incorrect", "Confident Correct"),
    contrast_row("Confident Incorrect", "Cautious Correct"),
    contrast_row("Cautious Incorrect",  "Cautious Correct"))   # Sec 3.4: within-cautious check
say(""); say("---- Latency separation: median difference + probability of superiority (95% CI) ----")
say(paste(capture.output(print(latency_separation, row.names = FALSE)), collapse = "\n"))

#### 8c. Confident-error slowing WITHIN each SA level (Section 3.5) ##########
# Rules out the confound that confident errors are slow only because they fall on
# inherently slower probe types. SA level is read from the audio file name in
# Event_ID: Cone queries are perception; other objects carry the suffix _P
# (perception), _C (comprehension), _PR (prediction). Descriptive only: the
# per-level error cells are small (12 / 7 / 10 of the 29 confident errors).
sa_level <- function(eid) {
    m <- regmatches(eid, regexec("/([A-Za-z]+)_([A-Za-z]+)\\.wav", eid))[[1]]
    if (length(m) < 3) return(NA_character_)
    if (tolower(m[2]) == "cone") return("Perception")
    unname(c(P = "Perception", C = "Comprehension", PR = "Prediction")[m[3]])
}
Real_Query_data$SA_Level <- vapply(Real_Query_data$Event_ID, sa_level, character(1))
within_level <- do.call(rbind, lapply(c("Perception", "Comprehension", "Prediction"), function(L) {
    sel <- Real_Query_data$SA_Level == L & Real_Query_data$Confidence >= 6
    cc <- Real_Query_data$Time_Difference_s[sel & Real_Query_data$Answer_Accuracy == 1]
    ci <- Real_Query_data$Time_Difference_s[sel & Real_Query_data$Answer_Accuracy == 0]
    data.frame(SA_level = L,
               n_conf_correct = length(cc),   median_correct_s   = round(median(cc), 2),
               n_conf_incorrect = length(ci), median_incorrect_s = round(median(ci), 2))
}))
say(""); say("---- Confident-error latency WITHIN each SA level (descriptive; small error cells) ----")
say(paste(capture.output(print(within_level, row.names = FALSE)), collapse = "\n"))

#### 8d. Operator-level driving-style correlations (Section 3.6) #############
# Telemetry summarised during ACTIVE driving in the Realtime condition (code B).
# Steering reversal rate at two gaps: fine (~4.5 deg) and major (~45 deg).
norm_id   <- function(x) sub("^participant_", "", x)
tele_rate <- function(f, col) read_csv(f, show_col_types = FALSE) %>%
    filter(Condition == "B") %>% transmute(P_ID = norm_id(Participant_ID), !!col := Reversal_rate)
speed_tbl <- read_csv("data/speed_metrics_active_driving_FIXED.csv", show_col_types = FALSE) %>%
    filter(Condition == "B") %>%
    transmute(P_ID = `Participant ID`, speed = `Average Speed (mph)`, speed_sd = `Speed StdDev (mph)`)
op_level <- Real_Query_data %>% group_by(P_ID) %>%
    summarise(accuracy = mean(Answer_Accuracy),
              overconf_errors = sum(Answer_Accuracy == 0 & Confidence >= 6), .groups = "drop") %>%
    left_join(speed_tbl, "P_ID") %>%
    left_join(tele_rate("data/active_steering_reversals_summary.csv",  "fine_srr"),    "P_ID") %>%
    left_join(tele_rate("data/active_steering_reversals_summary2.csv", "major_srr"),   "P_ID") %>%
    left_join(tele_rate("data/active_throttle_reversals_summary.csv",  "throttle_rr"), "P_ID")

corr_boot <- function(x, y, B = B_REP) {
    ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
    reps <- replicate(B, { i <- sample(n, replace = TRUE)
        suppressWarnings(cor(x[i], y[i], method = "spearman")) })
    c(rho = suppressWarnings(cor(x, y, method = "spearman")),
      lo = quantile(reps, .025, na.rm = TRUE), hi = quantile(reps, .975, na.rm = TRUE))
}
corr_row <- function(lab, x, y) { set.seed(2024); v <- corr_boot(x, y)
    data.frame(pair = lab, rho = round(v[1], 2), ci_lo = round(v[2], 2), ci_hi = round(v[3], 2), row.names = NULL) }
driving_correlations <- rbind(
    corr_row("speed ~ accuracy",        op_level$speed,       op_level$accuracy),
    corr_row("speed ~ overconf_errors", op_level$speed,       op_level$overconf_errors),
    corr_row("speed_sd ~ accuracy",     op_level$speed_sd,    op_level$accuracy),
    corr_row("major_srr ~ accuracy",    op_level$major_srr,   op_level$accuracy),
    corr_row("fine_srr ~ accuracy",     op_level$fine_srr,    op_level$accuracy),
    corr_row("throttle_rr ~ accuracy",  op_level$throttle_rr, op_level$accuracy))
say(""); say("---- Operator-level driving-style Spearman correlations (95% bootstrap CI) ----")
say(paste(capture.output(print(driving_correlations, row.names = FALSE)), collapse = "\n"))

#### 8e. Type-2 ROC / AUROC2, real-time vs offline (Section 3.8) #############
type2_auc <- function(correct, conf) prob_sup(conf[correct == 1], conf[correct == 0])
offline <- events %>% filter(Condition %in% c("C_L", "C_R")) %>%
    mutate(correct = as.integer(grepl("^correct", tolower(Answer_Accuracy)))) %>%
    filter(!is.na(Confidence))
rt_roc <- Real_Query_data %>% transmute(P_ID, correct = Answer_Accuracy, conf = Confidence)
of_roc <- offline         %>% transmute(P_ID, correct, conf = Confidence)
auc_boot <- function(df) { ids <- unique(df$P_ID); sp <- split(df, df$P_ID)
    reps <- replicate(B_REP, { dd <- bind_rows(sp[sample(ids, replace = TRUE)])
        type2_auc(dd$correct, dd$conf) })
    c(auc = type2_auc(df$correct, df$conf),
      lo = quantile(reps, .025, na.rm = TRUE), hi = quantile(reps, .975, na.rm = TRUE)) }
set.seed(2024); a_rt <- auc_boot(rt_roc)
set.seed(2024); a_of <- auc_boot(of_roc)
ids_both <- intersect(unique(rt_roc$P_ID), unique(of_roc$P_ID))
sp_rt <- split(rt_roc, rt_roc$P_ID); sp_of <- split(of_roc, of_roc$P_ID)
set.seed(2024)
diff_reps <- replicate(B_REP, { s <- sample(ids_both, replace = TRUE)
    type2_auc(bind_rows(sp_of[s])$correct, bind_rows(sp_of[s])$conf) -
    type2_auc(bind_rows(sp_rt[s])$correct, bind_rows(sp_rt[s])$conf) })
auroc2 <- data.frame(
    condition = c("real-time", "offline", "offline - real-time"),
    auc = round(c(a_rt[1], a_of[1], a_of[1] - a_rt[1]), 3),
    ci_lo = round(c(a_rt[2], a_of[2], quantile(diff_reps, .025, na.rm = TRUE)), 3),
    ci_hi = round(c(a_rt[3], a_of[3], quantile(diff_reps, .975, na.rm = TRUE)), 3), row.names = NULL) # nolint
say(""); say(sprintf("Offline probe: %d queries; accuracy %.1f%%; %% rated >=6 = %.1f%%",
            nrow(offline), 100 * mean(offline$correct), 100 * mean(offline$Confidence >= 6)))
say(sprintf("Offline point-biserial r(confidence, accuracy): %.3f",
            cor(offline$Confidence, offline$correct)))
say("---- Type-2 AUROC2 (confidence discriminating correct vs incorrect) ----")
say(paste(capture.output(print(auroc2, row.names = FALSE)), collapse = "\n"))

#### 8f. Per-operator profiles + operator 04 case (Sections 3.6, 3.9) ########
per_operator <- Real_Query_data %>% group_by(P_ID) %>%
    summarise(n = n(), accuracy = round(mean(Answer_Accuracy), 3),
              mean_conf = round(mean(Confidence), 2),
              overconf_errors = sum(Answer_Accuracy == 0 & Confidence >= 6),
              median_latency_s = round(median(Time_Difference_s), 2),
              pct_delayed = round(100 * mean(Time_Difference_s >= 3.5), 0), .groups = "drop") %>%
    left_join(select(op_level, P_ID, speed), "P_ID") %>% arrange(desc(overconf_errors))
oc_counts <- per_operator$overconf_errors
say(""); say(sprintf("Overconfident errors per operator: mean %.1f, SD %.1f, range %d-%d; %d of %d operators >=1",
            mean(oc_counts), sd(oc_counts), min(oc_counts), max(oc_counts), sum(oc_counts >= 1), length(oc_counts)))

#### 9. Write explicit outputs ##############################################
write_csv(Real_Query_data, file.path(out_dir, "confidence_query_level.csv"))
write_csv(calib,           file.path(out_dir, "confidence_calibration_summary.csv"))
write_csv(cal_indices,     file.path(out_dir, "confidence_calibration_indices.csv"))
write_csv(sweepA,          file.path(out_dir, "confidence_threshold_sweepA_overconf.csv"))
write_csv(sweepB,          file.path(out_dir, "confidence_threshold_sweepB_cautious.csv"))
write_csv(cat_summary,        file.path(out_dir, "confidence_category_summary.csv"))
write_csv(cal_indices_ci,     file.path(out_dir, "confidence_calibration_indices_ci.csv"))
write_csv(latency_separation, file.path(out_dir, "confidence_latency_separation.csv"))
write_csv(within_level,        file.path(out_dir, "confidence_within_level_latency.csv"))
write_csv(driving_correlations, file.path(out_dir, "confidence_driving_correlations.csv"))
write_csv(auroc2,             file.path(out_dir, "confidence_type2_auroc.csv"))
write_csv(op_level,           file.path(out_dir, "confidence_operator_level.csv"))
write_csv(per_operator,       file.path(out_dir, "confidence_per_operator.csv"))
writeLines(unlist(log_lines), file.path(out_dir, "confidence_numbers.txt"))

say("")
say(sprintf("Wrote outputs to %s/", out_dir))
