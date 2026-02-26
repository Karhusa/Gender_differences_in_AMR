BMI

## 1.2 Load packages

```r
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
library(mgcv)
library(ggpubr)
library(gt)
```

## 1.2 Load TSE

```r
tse_path <- "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/TSE.rds"
TSE <- readRDS(tse_path)
```

## 1.3 Extract colData

```r
colData_df <- as.data.frame(colData(TSE))

Subset <- colData_df %>%
  select(
    acc,
    ARG_div_shan,
    age_years,
    log10_ARG_load, 
    World_Bank_Income_Group,
    sex, 
    BMI_range_new, 
    Antibiotics_used, 
    UTI_history, 
    precise_age_category, imprecise_age_category
  )
```
## 1.4 Look through the columns:

```r
table(Subset$sex, Subset$BMI_range_new)
table(Subset$sex, Subset$BMI_range_new, Subset$World_Bank_Income_Group)
Subset <- Subset %>% filter(BMI_range_new != "Normal/Overweight (<30)")

```

| Sex    | Normal (18.5–25) | Normal/Overweight (<30) | Obese (>30) | Overweight (25–30) | Underweight (<18.5) |
| ------ | ---------------- | ----------------------- | ----------- | ------------------ | ------------------- |
| Female | 1593             | 32                      | 628         | 556                | 295                 |
| Male   | 1611             | 24                      | 581         | 823                | 170                 |



## 2 Ananlyses of ARG Load, BMI and sex

### 2.1 Boxplot of ARG Load, BMI and sex
```r

Subset <- Subset %>% filter(BMI_range_new != "Normal/Overweight (<30)")

Subset$sex <- recode(Subset$sex,
                     "female" = "Female",
                     "male"   = "Male")

levels_bmi <- c("Underweight (<18.5)", "Normal (18.5-25)", 
                "Overweight (25-30)", "Obese (>30)")
Subset$BMI_range_new <- factor(Subset$BMI_range_new, levels = levels_bmi)

plot_df <- Subset %>% 
  filter(!is.na(log10_ARG_load), !is.na(sex), !is.na(BMI_range_new)) %>%
  mutate(
    precise_age_category = factor(BMI_range_new, levels = levels_bmi),
    sex = factor(sex, levels = c("Female", "Male"))
  )

counts <- plot_df %>%
  group_by(BMI_range_new, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(
    y_pos = max(plot_df$log10_ARG_load, na.rm = TRUE) * 1.05,
    precise_age_category = factor(BMI_range_new, levels = levels_bmi)
  )

bmi_labels <- c("Underweight", "Normal", "Overweight", "Obese")
npg_cols <- pal_npg("nrc")(2)

# Plot
ggplot(plot_df, aes(x = BMI_range_new, y = log10_ARG_load, fill = sex)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    size = 1.2, alpha = 0.25, color = "grey30"
  ) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA, alpha = 0.8,
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    data = counts,
    aes(x = precise_age_category, y = y_pos, label = N, color = sex),
    position = position_dodge(width = 0.6),
    inherit.aes = FALSE,
    size = 2.0,
    fontface = "bold",
    show.legend = FALSE
  ) +
  stat_compare_means(
    aes(group = sex),           # group by sex for comparisons
    method = "t.test",          # or "wilcox.test" for non-parametric
    label = "p.signif",         # shows *, **, ***
    label.size = 4,
    hide.ns = TRUE              # hides ns (not significant)
  ) +
  scale_x_discrete(labels = bmi_labels) +
  scale_fill_npg() +
  scale_color_npg() +
  labs(
    title = "ARG Load by BMI and Sex",
    x = "BMI range",
    y = expression(log[10]*"(ARG load)"),
    fill = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
```
### 2.4.2 Additive linear model of AGR Load, BMI and sex
```r
plot_df$BMI_range_new <- relevel(
  plot_df$BMI_range_new,
  ref = "Normal (18.5-25)"
)

model_add <- lm(log10_ARG_load ~ BMI_range_new + sex, data = plot_df)
summary(model_add)
```

### 2.4.3 Separate by income gruops
```
df_income_ARG <- Subset %>%
  select(
    sex,
    BMI_range_new,            
    World_Bank_Income_Group,
    log10_ARG_load
  ) %>%
  mutate(
    Income_group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c(
        "Low income",
        "Lower middle income",
        "Upper middle income"
      ) ~ "LMIC",
      TRUE ~ NA_character_
    )
  ) %>%
  select(sex, BMI_range_new, Income_group, log10_ARG_load)

df_income_ARG$Income_group <- factor(
  df_income_ARG$Income_group,
  levels = c("HIC", "LMIC")  
)

```
### 2.4.3 Boxplot
```
bmi_labels <- c("Underweight", "Normal", "Overweight", "Obese")

levels_bmi <- c("Underweight (<18.5)", "Normal (18.5-25)", 
                "Overweight (25-30)", "Obese (>30)")

df_income_ARG <- df_income_ARG %>%
  filter(!is.na(log10_ARG_load), !is.na(sex), 
         !is.na(BMI_range_new), !is.na(Income_group)) %>%
  mutate(
    BMI_range_new = factor(BMI_range_new, levels = levels_bmi),
    sex = factor(sex, levels = c("Female", "Male")),
    Income_group = factor(Income_group, levels = c("HIC", "LMIC"))
  )

ggplot(df_income_ARG, aes(x = Income_group, y = log10_ARG_load, fill = sex)) +

  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    size = 1.2, alpha = 0.25, color = "grey30"
  ) +
  
  geom_boxplot(
    width = 0.65,  # wider boxes
    outlier.shape = NA, alpha = 0.8,
    position = position_dodge(width = 0.65)
  ) +
  
  stat_compare_means(
    aes(group = sex),
    method = "t.test",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  
  facet_wrap(
    ~ BMI_range_new,
    nrow = 1,
    scales = "free_x",
    labeller = label_wrap_gen(width = 15)  # wrap long facet labels
  ) +
  
  scale_fill_npg() +
  scale_color_npg() +
  
  labs(
    title = "ARG Load by Income, Sex, and BMI",
    x = "Income group",
    y = expression(log[10]*"(ARG load)"),
    fill = "Sex"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),  # lighter background
    strip.text = element_text(face = "bold", size = 11, lineheight = 0.9),
    panel.grid.major.y = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )
```
Remove ages under 18

```
df_income_BMI_ARG <- Subset %>%
  filter(age_years >= 18) %>%           # Remove ages under 18
  select(
    sex,
    age_years,
    BMI_range_new,            
    World_Bank_Income_Group,
    log10_ARG_load
  ) %>%
  mutate(
    Income_group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c(
        "Low income",
        "Lower middle income",
        "Upper middle income"
      ) ~ "LMIC",
      TRUE ~ NA_character_
    )
  ) %>%
  select(sex, BMI_range_new, Income_group, log10_ARG_load)

df_income_BMI_ARG$Income_group <- factor(
  df_income_BMI_ARG$Income_group,
  levels = c("HIC", "LMIC")  
)


levels_bmi <- c("Underweight (<18.5)", "Normal (18.5-25)", 
                "Overweight (25-30)", "Obese (>30)")

df_income_BMI_ARG <- df_income_BMI_ARG %>%
  filter(!is.na(log10_ARG_load), !is.na(sex), 
         !is.na(BMI_range_new), !is.na(Income_group)) %>%
  mutate(
    BMI_range_new = factor(BMI_range_new, levels = levels_bmi),
    sex = factor(sex, levels = c("Female", "Male")),
    Income_group = factor(Income_group, levels = c("HIC", "LMIC"))
  )

ggplot(df_income_BMI_ARG, aes(x = Income_group, y = log10_ARG_load, fill = sex)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    size = 1.2, alpha = 0.25, color = "grey30"
  ) +
  geom_boxplot(
    width = 0.65,  # wider boxes
    outlier.shape = NA, alpha = 0.8,
    position = position_dodge(width = 0.65)
  ) +
  stat_compare_means(
    aes(group = sex),
    method = "t.test",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  
  facet_wrap(
    ~ BMI_range_new,
    nrow = 1,
    scales = "free_x",
    labeller = label_wrap_gen(width = 15)  # wrap long facet labels
  ) +
  
  scale_fill_npg() +
  scale_color_npg() +
  
  labs(
    title = "ARG Load by Income, Sex, and BMI",
    x = "Income group",
    y = expression(log[10]*"(ARG load)"),
    fill = "Sex"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),  # lighter background
    strip.text = element_text(face = "bold", size = 11, lineheight = 0.9),
    panel.grid.major.y = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )
```



```

df_income_BMI_ARG <- df_income_BMI_ARG %>%
  mutate(BMI_range_new = relevel(BMI_range_new, ref = "Normal (18.5-25)"))

# Linear model for HIC
model_HIC <- lm(
  log10_ARG_load ~ BMI_range_new + sex,
  data = df_income_BMI_ARG %>% filter(Income_group == "HIC")
)
summary(model_HIC)

# Linear model for LMIC
model_LMIC <- lm(
  log10_ARG_load ~ BMI_range_new + sex,
  data = df_income_BMI_ARG %>% filter(Income_group == "LMIC")
)
summary(model_LMIC)

```



Sample counts:

```r

sample_counts <- df_income_BMI_ARG %>%
  group_by(Income_group, BMI_range_new, sex) %>%
  summarise(n_samples = n(), .groups = "drop")

# Display as a separate gt table
sample_counts %>%
  gt(groupname_col = "Income_group") %>%
  tab_header(
    title = "Number of Samples Used in Linear Models",
    subtitle = "Counts per BMI category and sex"
  ) %>%
  cols_label(
    BMI_range_new = "BMI Category",
    sex = "Sex",
    n_samples = "Number of Samples"
  )

```

```
get_signif <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ ".",
    TRUE      ~ ""
  )
}

# Extract HIC model summary
sum_HIC <- summary(model_HIC)
df_HIC <- data.frame(
  term = rownames(sum_HIC$coefficients),
  estimate = sum_HIC$coefficients[, "Estimate"],
  std.error = sum_HIC$coefficients[, "Std. Error"],
  t.value = sum_HIC$coefficients[, "t value"],
  p.value = sum_HIC$coefficients[, "Pr(>|t|)"],
  Income_group = "HIC",
  row.names = NULL
)

# Extract LMIC model summary
sum_LMIC <- summary(model_LMIC)
df_LMIC <- data.frame(
  term = rownames(sum_LMIC$coefficients),
  estimate = sum_LMIC$coefficients[, "Estimate"],
  std.error = sum_LMIC$coefficients[, "Std. Error"],
  t.value = sum_LMIC$coefficients[, "t value"],
  p.value = sum_LMIC$coefficients[, "Pr(>|t|)"],
  Income_group = "LMIC",
  row.names = NULL
)

# Combine, clean names, and add significance
df_combined <- bind_rows(df_HIC, df_LMIC) %>%
  mutate(
    term = recode(term,
                  "(Intercept)" = "Normal (18.5-25) Female",
                  "BMI_range_newUnderweight (<18.5)" = "Underweight (<18.5)",
                  "BMI_range_newOverweight (25-30)" = "Overweight (25-30)",
                  "BMI_range_newObese (>30)" = "Obese (>30)",
                  "sexMale" = "Male"),
    Significance = get_signif(p.value)
  )

# Create gt table
df_combined %>%
  gt(groupname_col = "Income_group") %>%
  fmt_number(
    columns = c(estimate, std.error, t.value, p.value),
    decimals = 3
  ) %>%
  tab_header(
    title = "Linear Model Results: log10(ARG load) ~ BMI + Sex",
    subtitle = "Reference: Normal weight Female"
  ) %>%
  cols_label(
    term = "Predictor",
    estimate = "Estimate",
    std.error = "Std. Error",
    t.value = "t-value",
    p.value = "p-value",
    Significance = "Significance"
  )


```
