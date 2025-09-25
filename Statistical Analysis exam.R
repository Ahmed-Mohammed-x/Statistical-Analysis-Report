# Clean end-to-end script with H6 and H7 integrated (single cell)
# Note: Original logic preserved; H6 and H7 are added within Hypotheses section.

# packages we need
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("forcats")
install.packages("lubridate")
install.packages("broom")
install.packages("scales")
install.packages("tibble")

# Load packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(lubridate)
library(broom)
library(scales)
library(tibble)

# 1) Load the dataset that we will work on
sales_df <- read.csv("fashion_retail_sales.csv", check.names = FALSE)
stopifnot(is.data.frame(sales_df))

# 1a) Robust date parsing (tries ymd, then mdy, then dmy)
parse_date_vector <- function(x) {
  y <- suppressWarnings(lubridate::ymd(x, quiet = TRUE))
  ok_y <- !is.na(y)
  if (all(ok_y)) return(as_date(y))
  m <- suppressWarnings(lubridate::mdy(x, quiet = TRUE))
  out <- ifelse(ok_y, y, m)
  ok_out <- !is.na(out)
  if (all(ok_out)) return(as_date(out))
  d <- suppressWarnings(lubridate::dmy(x, quiet = TRUE))
  out2 <- ifelse(ok_out, out, d)
  as_date(out2)
}

# 2) Standardize columns and keep customer id for a dedup key
sales_cleaned <- sales_df %>%
  mutate(
    amount_usd    = `Purchase Amount (USD)`,
    rating_rnum   = suppressWarnings(as.numeric(`Review Rating`)),
    payment_method = factor(`Payment Method`, levels = c("Cash","Credit Card")),
    item          = factor(`Item Purchased`),
    pur_date      = parse_date_vector(`Date Purchase`),
    cust_id       = `Customer Reference ID`
  ) %>%
  select(amount_usd, rating_rnum, payment_method, item, pur_date, cust_id)

# 3) Missingness and deduplication
miss_col1 <- sapply(sales_cleaned, function(x) sum(is.na(x)))
miss_tbl1 <- tibble(var = names(miss_col1), n_miss = as.integer(miss_col1)) %>%
  mutate(pct_miss = round(100 * n_miss / nrow(sales_cleaned), 1)) %>% arrange(desc(pct_miss))

# Composite key: cust_id | pur_date | item
`1key` <- paste(sales_cleaned$cust_id, sales_cleaned$pur_date, sales_cleaned$item, sep = "|")
clean_deduplicate <- sales_cleaned[!duplicated(`1key`), ]

# 4) Analysis subsets
sales_analysis <- clean_deduplicate %>%
  filter(!is.na(amount_usd)) %>%
  mutate(
    zero_amount = amount_usd == 0,
    high_amount = amount_usd >= 1000
  )

sales_dated_dt <- sales_analysis %>% filter(!is.na(pur_date))
rated_sales    <- sales_analysis %>% filter(!is.na(rating_rnum))

# 5) KPIs and descriptive tables
sales_pki       <- sales_analysis %>% summarise(orders = n(), revenue = sum(amount_usd), AOV = mean(amount_usd))
pay_method_mix  <- sales_analysis %>% count(payment_method, name = "orders") %>% mutate(share = percent(orders / sum(orders)))
pay_method_aov  <- sales_analysis %>% group_by(payment_method) %>% summarise(AOV = mean(amount_usd), revenue = sum(amount_usd), .groups = "drop")

top_items <- sales_analysis %>%
  group_by(item) %>%
  summarise(orders = n(), revenue = sum(amount_usd), AOV = mean(amount_usd), avg_rating = mean(rating_rnum, na.rm = TRUE), .groups = "drop") %>%
  slice_max(revenue, n = 10)

monthly_sum <- sales_dated_dt %>%
  mutate(ym = floor_date(pur_date, "month")) %>%
  group_by(ym) %>% summarise(orders = n(), revenue = sum(amount_usd), AOV = mean(amount_usd), .groups = "drop")

# 6) Plots for exploration
plot_amount_hist <- ggplot(sales_analysis, aes(amount_usd)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white") +
  scale_x_continuous(labels = scales::dollar) +
  labs(title = "Purchase Amount Distribution", x = "Amount (USD)", y = "Count") +
  theme_minimal()

h1_plot_aov_by_payment <- pay_method_aov %>%
  ggplot(aes(payment_method, AOV, fill = payment_method)) +
  geom_col(width = 0.6) +
  scale_y_continuous(labels = scales::dollar) +
  labs(title = "Average Order Value by Payment", x = NULL, y = "AOV (USD)") +
  theme_minimal() + theme(legend.position = "none")

plot_monthly_trends <- monthly_sum %>%
  pivot_longer(cols = c(orders, revenue, AOV), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = ym, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(title = "Monthly Orders, Revenue, and AOV (Dated Subset)", x = NULL, y = NULL) +
  theme_minimal()

# 7) Hypotheses and tests
# H1: AOV by payment (Welch); sensitivity excluding high outliers
h1_welch_aov_by_payment_all  <- t.test(amount_usd ~ payment_method, data = sales_analysis, var.equal = FALSE)
h1_welch_aov_by_payment_nohi <- t.test(amount_usd ~ payment_method, data = dplyr::filter(sales_analysis, !high_amount), var.equal = FALSE)

# H2 and H3: Item differences on top-10 items
items_10_vec  <- sales_analysis %>% count(item, sort = TRUE) %>% slice_head(n = 10) %>% pull(item)
sales_topitems <- sales_analysis %>% filter(item %in% items_10_vec)

# H2: rating_rnum ~ item (ANOVA)
h2_anova_rating_item <- aov(rating_rnum ~ item, data = sales_topitems)
# H3: amount_usd ~ item (ANOVA + Kruskal)
h3_anova_amount_item  <- aov(amount_usd ~ item, data = sales_topitems)
h3_kruskal_amount_item <- kruskal.test(amount_usd ~ item, data = sales_topitems)

# Post-hoc for H3 if significant
h3_posthoc_hsd_amount_item <- tryCatch({
  if (summary(h3_anova_amount_item)[[1]][["Pr(>F)"]][1] < 0.05) TukeyHSD(h3_anova_amount_item) else NULL
}, error = function(e) NULL)

# H4: Q4 effect — exclude partial 2023-10 (if present), Welch + Wilcoxon
cutoff_2023_10_01 <- as.Date("2023-10-01")
sales_pre_q4 <- sales_dated_dt %>%
  filter(is.na(pur_date) | pur_date < cutoff_2023_10_01) %>%
  mutate(month = month(pur_date), q4 = factor(month %in% c(10,11,12), levels = c(FALSE, TRUE), labels = c("Non-Q4","Q4")))

h4_welch_q4_vs_nonq4  <- t.test(amount_usd ~ q4, data = sales_pre_q4)
h4_wilcox_q4_vs_nonq4 <- wilcox.test(amount_usd ~ q4, data = sales_pre_q4)

# H5: Multivariate regression — amount_usd ~ rating_rnum + payment_method + item (exclude zeros)
h5_lm_amount_drivers <- lm(amount_usd ~ rating_rnum + payment_method + item, data = sales_analysis %>% filter(!zero_amount))

# H6: Correlation between amount_usd and rating_rnum among rated, non-zero transactions
h6_df <- rated_sales %>% filter(!zero_amount)
# Pearson and Spearman correlations
h6_cor_pearson_test  <- suppressWarnings(cor.test(h6_df$amount_usd, h6_df$rating_rnum, method = "pearson"))
h6_cor_spearman_test <- suppressWarnings(cor.test(h6_df$amount_usd, h6_df$rating_rnum, method = "spearman", exact = FALSE))
# Plot
h6_scatter <- ggplot(h6_df, aes(x = rating_rnum, y = amount_usd)) +
  geom_point(alpha = 0.25, color = "#256fa1") +
  geom_smooth(method = "lm", se = TRUE, color = "#e1361d") +
  scale_y_continuous(labels = dollar) +
  labs(title = "H6: Amount vs Rating", x = "Review rating", y = "Amount (USD)") +
  theme_minimal()

# H7: Simple linear regression focused on rating
# Model 1: amount_usd ~ rating_rnum
h7_lm_simple   <- lm(amount_usd ~ rating_rnum, data = h6_df)
# Model 2: with controls
h7_lm_controls <- lm(amount_usd ~ rating_rnum + payment_method + item, data = h6_df)

# 8) Console summary for professor :)

cat("\n========== SUMMARY ==========\n")

cat("\n======= Missingness (top 10 vars) [preview] ======\n")
print(head(miss_tbl1, 10), row.names = FALSE)

cat("\n======= KPIs ======\n")
print(sales_pki, row.names = FALSE)

cat("\n===== Payment mix ======\n")
print(pay_method_mix, row.names = FALSE)

cat("\n=== =AOV and revenue by payment ====\n")
print(pay_method_aov, row.names = FALSE)

cat("\n==== Top items by revenue (top 10) ====\n")
print(top_items, row.names = FALSE)

cat("\n==== Monthly rollup [first 12 rows] ====\n")
print(head(monthly_sum, 12), row.names = FALSE)

cat("\n======= HYPOTHESES=======\n")

cat("\n==== H1 Welch (all) ====\n")
print(broom::tidy(h1_welch_aov_by_payment_all), row.names = FALSE)

cat("\n==== H1 Welch (no high outliers) ====\n")
print(broom::tidy(h1_welch_aov_by_payment_nohi), row.names = FALSE)

cat("\n==== H2 ANOVA: rating ~ item (top 10) ====\n")
print(broom::tidy(h2_anova_rating_item), row.names = FALSE)

cat("\n==== H3 ANOVA: amount ~ item (top 10) ===\n")
print(broom::tidy(h3_anova_amount_item), row.names = FALSE)

cat("\n==== H3 Kruskal: amount ~ item (top 10) ====\n")
print(h3_kruskal_amount_item)

if (!is.null(h3_posthoc_hsd_amount_item)) {
  cat("\n==== H3 Tukey post-hoc (significant ANOVA) ====\n")
  print(h3_posthoc_hsd_amount_item)
}

cat("\n==== H4 Q4 vs Non-Q4 (Welch) ====\n")
print(broom::tidy(h4_welch_q4_vs_nonq4), row.names = FALSE)

cat("\n=== H4 Q4 vs Non-Q4 (Wilcoxon) ===\n")
print(h4_wilcox_q4_vs_nonq4)

cat("\n==== H5 OLS (multivariate) summary ====\n")
print(summary(h5_lm_amount_drivers))

cat("\n======= H6 Correlation (Pearson) amount_usd ~ rating_rnum ======\n")
print(h6_cor_pearson_test)

cat("\n===== H6 Correlation (Spearman) amount_usd ~ rating_rnum =====\n")
print(h6_cor_spearman_test)

cat("\n====== H7 OLS (simple): amount_usd ~ rating_rnum ======\n")
print(broom::tidy(h7_lm_simple), row.names = FALSE)

cat("\n==== H7 OLS (simple) summary ====\n")
print(summary(h7_lm_simple))

cat("\n=== H7 OLS (with controls): amount_usd ~ rating_rnum + payment_method + item ===\n")
print(broom::tidy(h7_lm_controls), row.names = FALSE)

cat("\n=== H7 OLS (with controls) summary ===\n")
print(summary(h7_lm_controls))

cat("\n========== DIAGNOSTICS ==========\n")

par(mfrow = c(2, 2)); plot(h5_lm_amount_drivers)
cat("\n=-== Shapiro-Wilk on residuals (H5 model) =-==\n")
print(shapiro.test(residuals(h5_lm_amount_drivers)))

par(mfrow = c(2, 2)); plot(h7_lm_simple)
cat("\n==== Shapiro-Wilk on residuals (H7 simple model) -===\n")
print(shapiro.test(residuals(h7_lm_simple)))

cat("\n==== Selected plots ====\n")
plot_amount_hist; h1_plot_aov_by_payment; plot_monthly_trends; h6_scatter