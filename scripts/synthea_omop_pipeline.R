## Synthea OMOP (CSV on S3) cohort -> trajectories -> LCA/GBTM -> prediction -> report
##
## Dataset: s3://synthea-omop/synthea1k/*.csv  (no concept tables provided)
## Strategy:
##   - Read OMOP CDM CSVs directly from S3 via DuckDB httpfs + read_csv_auto()
##   - Define cohort using *source_value* text matching (e.g., "metformin", "type 2 diabetes")
##   - Build monthly PDC trajectories (12 months) for new users (180d washout)
##   - Run LCA (poLCA) on binary monthly coverage
##   - Run GBTM (gbmt) on continuous PDC
##   - Predict class membership with baseline covariates (age/sex + optional measurements)
##   - Export tidy outputs + plots; optionally render Quarto report (template provided)

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(tidyr)
  library(lubridate)
  library(stringr)
  library(ggplot2)
  library(rsample)
  library(yardstick)
  library(nnet)
  library(poLCA)
  library(gbmt)
  library(dbplyr)
  library(dplyr)
})

# ---------------------------
# 0) Parameters (edit here)
# ---------------------------
S3_BASE <- "https://synthea-omop.s3.amazonaws.com/synthea1k"  # CSVs live here
TARGET_DRUG_PATTERN <- "metformin"                           # case-insensitive grep on drug_source_value
TARGET_COND_PATTERN <- "type 2 diabetes"                     # grep on condition_source_value
WASHOUT_DAYS <- 180
FOLLOWUP_MONTHS <- 12
REFILL_GAP_DAYS <- 60                                        # for attrition (optional)
SEED <- 123

OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("Output dir: ", normalizePath(OUT_DIR))

# ---------------------------
# 1) DuckDB connection + views
# ---------------------------
con <- dbConnect(duckdb::duckdb(), dbdir = "synthea_omop.duckdb", read_only = FALSE)
dbExecute(con, "INSTALL httpfs;")
dbExecute(con, "LOAD httpfs;")

mk_view_csv <- function(view_name, filename){
  sql <- sprintf(
    "CREATE OR REPLACE VIEW %s AS SELECT * FROM read_csv_auto('%s/%s', HEADER=TRUE, SAMPLE_SIZE=-1);",
    view_name, S3_BASE, filename
  )
  dbExecute(con, sql)
}

# Required tables for this pipeline
mk_view_csv("person", "person.csv")
mk_view_csv("drug_exposure", "drug_exposure.csv")
mk_view_csv("condition_occurrence", "condition_occurrence.csv")
mk_view_csv("visit_occurrence", "visit_occurrence.csv")
# Optional: measurement (big); comment out if slow at first
mk_view_csv("measurement", "measurement.csv")

# Sanity checks
message("Rows(person) = ", dbGetQuery(con, "SELECT COUNT(*) n FROM person")$n)
message("Rows(drug_exposure) = ", dbGetQuery(con, "SELECT COUNT(*) n FROM drug_exposure")$n)

# ---------------------------
# 2) Cohort definition (new users of TARGET_DRUG_PATTERN with TARGET_COND_PATTERN)
#    Uses source_value strings because vocab tables not included.
# ---------------------------

# Find candidate drug rows
# drug_exposure expected columns: person_id, drug_exposure_start_date, days_supply, drug_source_value, ...
# We'll defensively coalesce days_supply to 30 if missing.

TARGET_DRUG_CODE <- "140"      # simvastatin (RxNorm-related)
TARGET_COND_CODE <- "162864005"   # obesity (SNOMED CT)

# drug_target
dbExecute(con, sprintf("
  CREATE OR REPLACE TEMP TABLE drug_target AS
  SELECT
    person_id,
    CAST(drug_exposure_start_date AS DATE) AS start_date,
    CAST(COALESCE(days_supply, 30) AS INTEGER) AS days_supply,
    CAST(drug_source_value AS VARCHAR) AS drug_code
  FROM drug_exposure
  WHERE CAST(drug_source_value AS VARCHAR) = '%s'
", TARGET_DRUG_CODE))

# index date
dbExecute(con, "
  CREATE OR REPLACE TEMP TABLE index_drug AS
  SELECT person_id, MIN(start_date) AS index_date
  FROM drug_target
  GROUP BY person_id
")

# washout
dbExecute(con, sprintf("
  CREATE OR REPLACE TEMP TABLE cohort_drug AS
  SELECT i.person_id, i.index_date
  FROM index_drug i
  WHERE NOT EXISTS (
    SELECT 1
    FROM drug_target d
    WHERE d.person_id = i.person_id
      AND d.start_date < i.index_date
      AND d.start_date >= i.index_date - INTERVAL %d DAY
  )
", WASHOUT_DAYS))

# # condition
# dbExecute(con, sprintf("
#   CREATE OR REPLACE TEMP TABLE cond_target AS
#   SELECT
#     person_id,
#     CAST(condition_start_date AS DATE) AS cond_date,
#     CAST(condition_source_value AS VARCHAR) AS cond_code
#   FROM condition_occurrence
#   WHERE CAST(condition_source_value AS VARCHAR) = '%s'
# ", TARGET_COND_CODE))

dbExecute(con, "
  CREATE OR REPLACE TEMP TABLE cohort AS
  SELECT * FROM cohort_drug
")


cohort_n <- dbGetQuery(con, "SELECT COUNT(*) n FROM cohort")$n
message("Cohort size = ", cohort_n)

# ---------------------------
# 3) Baseline covariates (age/sex) + optional measurement summaries
# ---------------------------
base <- tbl(con, "person") |>
  inner_join(tbl(con, "cohort"), by="person_id") |>
  transmute(
    person_id,
    gender_concept_id = as.integer(gender_concept_id),
    year_of_birth = as.integer(year_of_birth),
    index_date = index_date
  ) |>
  collect()

base <- base |>
  dplyr::mutate(
    sex = factor(gender_concept_id),
    age = as.integer(lubridate::year(index_date) - year_of_birth)
  ) |>
  dplyr::select(person_id, sex, age)

# Optional: baseline measurement features (example: BMI or HbA1c if present in source_value)
# Without concept tables, we use measurement_source_value string patterns.
meas_feats <- function(pattern, feat_name){
  sql <- sprintf("
    SELECT
      c.person_id,
      AVG(CAST(m.value_as_number AS DOUBLE)) AS %s
    FROM measurement m
    JOIN cohort c ON c.person_id = m.person_id
    WHERE CAST(m.measurement_date AS DATE) BETWEEN c.index_date - INTERVAL 180 DAY AND c.index_date - INTERVAL 1 DAY
      AND LOWER(COALESCE(m.measurement_source_value, '')) LIKE '%%%s%%'
      AND m.value_as_number IS NOT NULL
    GROUP BY 1
  ", feat_name, pattern)
  dbGetQuery(con, sql)
}

# Comment these in if you confirm the source_value strings exist in your dataset
# bmi <- meas_feats("bmi", "bmi_baseline")
# a1c <- meas_feats("a1c", "a1c_baseline")
# base <- base |> left_join(bmi, by="person_id") |> left_join(a1c, by="person_id")

# ---------------------------
# 4) Build monthly windows (FOLLOWUP_MONTHS) and compute PDC
# ---------------------------
dbExecute(con, sprintf("
  CREATE OR REPLACE TEMP TABLE win AS
  SELECT
    person_id,
    index_date,
    i AS m,
    index_date + (i || ' months')::INTERVAL AS win_start,
    index_date + ((i+1) || ' months')::INTERVAL AS win_end
  FROM cohort
  CROSS JOIN range(0,%d) r(i)
", FOLLOWUP_MONTHS))

# Exposure intervals for target drug
dbExecute(con, "
  CREATE OR REPLACE TEMP TABLE exp AS
  SELECT
    person_id,
    start_date AS s,
    (start_date + (days_supply || ' days')::INTERVAL) AS e
  FROM drug_target
  WHERE person_id IN (SELECT person_id FROM cohort)
")

# Overlap days and PDC
dbExecute(con, "
  CREATE OR REPLACE TEMP TABLE traj_long AS
  SELECT
    w.person_id,
    w.m,
    CAST(w.win_start AS DATE) AS win_start,
    CAST(w.win_end   AS DATE) AS win_end,
    SUM(
      GREATEST(
        0,
        DATE_DIFF('day',
          GREATEST(CAST(w.win_start AS DATE), CAST(e.s AS DATE)),
          LEAST(CAST(w.win_end   AS DATE), CAST(e.e AS DATE))
        )
      )
    ) AS covered_days,
    DATE_DIFF('day', CAST(w.win_start AS DATE), CAST(w.win_end AS DATE)) AS win_days
  FROM win w
  LEFT JOIN exp e
    ON e.person_id = w.person_id
   AND CAST(e.e AS DATE) > CAST(w.win_start AS DATE)
   AND CAST(e.s AS DATE) < CAST(w.win_end   AS DATE)
  GROUP BY 1,2,3,4,6
")

traj <- tbl(con, "traj_long") |>
  dplyr::mutate(
    pdc = pmin(covered_days / win_days, 1.0),
    covered = as.integer(pdc > 0)
  ) |>
  dplyr::select(person_id, m, pdc, covered) |>
  collect()

write.csv(traj, file.path(OUT_DIR, "traj_long.csv"), row.names = FALSE)

# ---------------------------
# 5) EDA plots
# ---------------------------
p_persist <- traj |>
  group_by(m) |>
  summarise(persistence = mean(covered), .groups="drop") |>
  ggplot(aes(m, persistence)) +
  geom_line() + geom_point() +
  scale_y_continuous(limits=c(0,1)) +
  labs(x="Month since index", y="Proportion covered (PDC>0)",
       title=paste0("Persistence: ", TARGET_DRUG_CODE))

ggsave(file.path(OUT_DIR, "persistence.png"), p_persist, width=7, height=4, dpi=200)

set.seed(SEED)
ids <- sample(unique(traj$person_id), min(200, length(unique(traj$person_id))))
p_spag <- traj |>
  filter(person_id %in% ids) |>
  ggplot(aes(m, pdc, group=person_id)) +
  geom_line(alpha=0.15) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x="Month since index", y="PDC", title="PDC trajectories (sample)")

ggsave(file.path(OUT_DIR, "pdc_spaghetti.png"), p_spag, width=7, height=4, dpi=200)



# ---------------------------
# 6) LCA (poLCA) on binary monthly coverage (robust)
# ---------------------------
traj_wide <- traj |>
  dplyr::mutate(var = paste0("m", m+1)) |>
  dplyr::select(person_id, var, covered) |>
  tidyr::pivot_wider(names_from = var, values_from = covered) |>
  tidyr::drop_na()

month_cols <- names(traj_wide)[-1]
keep_cols <- month_cols[sapply(traj_wide[month_cols], function(x) length(unique(x)) >= 2)]
traj_wide2 <- traj_wide |> dplyr::select(person_id, dplyr::all_of(keep_cols))

if (nrow(traj_wide2) < 50 || length(keep_cols) < 2) {
  message(sprintf("Skip LCA: n=%d, usable_months=%d", nrow(traj_wide2), length(keep_cols)))
  bic_tbl <- tibble(k=NA_integer_, bic=NA_real_, aic=NA_real_)
  class_assign <- tibble(person_id = traj_wide$person_id, lca_class = factor(1))
  K_lca <- 1L
} else {
  
  lca_dat <- traj_wide2 |>
    dplyr::mutate(dplyr::across(dplyr::all_of(keep_cols), ~ .x + 1L))  # 0/1 -> 1/2
  
  f <- as.formula(paste0("cbind(", paste(keep_cols, collapse=","), ") ~ 1"))
  
  set.seed(SEED)
  fits <- lapply(2:5, function(k){
    tryCatch(
      poLCA(f, data=lca_dat, nclass=k, nrep=10, verbose=FALSE),
      error = function(e) e
    )
  })
  
  ok <- sapply(fits, function(x) !inherits(x, "error"))
  if (!any(ok)) stop("All LCA fits failed. Increase cohort or reduce months.")
  
  fits_ok <- fits[ok]
  ks_ok <- (2:5)[ok]
  
  bic_tbl <- tibble(
    k = ks_ok,
    bic = sapply(fits_ok, function(m) m$bic),
    aic = sapply(fits_ok, function(m) m$aic)
  ) |>
    dplyr::arrange(bic)
  
  best_lca <- fits_ok[[which.min(bic_tbl$bic)]]
  K_lca <- bic_tbl$k[which.min(bic_tbl$bic)]
  
  class_assign <- tibble(
    person_id = traj_wide2$person_id,
    lca_class = factor(best_lca$predclass)
  )
}

write.csv(bic_tbl, file.path(OUT_DIR, "lca_model_selection.csv"), row.names = FALSE)
write.csv(class_assign, file.path(OUT_DIR, "lca_classes.csv"), row.names = FALSE)



# ---------------------------
# 7) GBTM (gbmt) on continuous PDC
# ---------------------------
traj_for_gbtm <- traj |>
  dplyr::transmute(id = as.integer(person_id),
                   time = as.integer(m + 1),
                   pdc = as.numeric(pdc)) |>
  dplyr::arrange(id, time) |>
  as.data.frame()

# gbmt can be fragile; fit with tryCatch and keep only successful models
set.seed(SEED)

fit_one <- function(k){
  tryCatch(
    gbmt(x.names = "pdc", unit="id", time="time", ng=k, d=2,
         data=traj_for_gbtm, quiet=TRUE),
    error = function(e) e
  )
}

fits <- lapply(2:5, fit_one)
ok <- sapply(fits, function(x) !inherits(x, "error"))

if (!any(ok)) {
  message("All GBTM fits failed; skipping GBTM outputs.")
} else {
  fits_ok <- fits[ok]
  ks_ok <- (2:5)[ok]
  
  # Coerce loglik to a single numeric value per model
  get_loglik <- function(m){
    ll <- m$loglik
    if (is.list(ll)) ll <- unlist(ll)
    ll <- as.numeric(ll)
    if (length(ll) == 0) NA_real_ else ll[1]
  }
  
  logliks <- tibble(
    k = ks_ok,
    logLik = sapply(fits_ok, get_loglik)
  ) |>
    dplyr::arrange(dplyr::desc(logLik))
  
  write.csv(logliks, file.path(OUT_DIR, "gbtm_model_selection.csv"), row.names = FALSE)
  
  # robust best-model selection even if logLik is missing/NA
  ll_vec <- logliks$logLik
  if (all(is.na(ll_vec))) {
    message("GBTM: logLik not available (all NA). Using the first successful model for plotting.")
    best_gb <- fits_ok[[1]]
  } else {
    best_idx <- which.max(ifelse(is.na(ll_vec), -Inf, ll_vec))
    best_gb <- fits_ok[[best_idx]]
  }
  
  png(file.path(OUT_DIR, "gbtm_plot.png"), width=900, height=600)
  plot(best_gb)
  dev.off()
}





# ---------------------------
# 8) Prediction: multinomial logistic regression for LCA class
# ---------------------------
model_df <- base |>
  inner_join(class_assign, by="person_id") |>
  mutate(lca_class = factor(lca_class)) |>
  filter(!is.na(age), !is.na(sex))

set.seed(SEED)
split <- initial_split(model_df, prop=0.8, strata=lca_class)
train <- training(split)
test  <- testing(split)

fit_mn <- nnet::multinom(lca_class ~ age + sex, data=train, trace=FALSE)
pred <- predict(fit_mn, newdata=test)

f1_macro <- yardstick::f_meas_vec(truth=test$lca_class, estimate=pred, estimator="macro")
acc <- yardstick::accuracy_vec(truth=test$lca_class, estimate=pred)

metrics <- tibble(metric=c("macro_f1","accuracy"),
                  value=c(as.numeric(f1_macro), as.numeric(acc)))
write.csv(metrics, file.path(OUT_DIR, "prediction_metrics.csv"), row.names = FALSE)

message("Done. Outputs written to: ", normalizePath(OUT_DIR))
message(sprintf("Best LCA K = %d | macro-F1 = %.3f | acc = %.3f", K_lca, f1_macro, acc))

# ---------------------------
# 9) Close DB
# ---------------------------
dbDisconnect(con, shutdown=TRUE)







