# 1. Prepare the data
## 1.1 Load packages
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
## 1.4 Look through columns

```
table(Subset$sex)
table(Subset$sex, Subset$precise_age_category)
Subset <- Subset[Subset$precise_age_category != "Unknown", ]
```
**Sex**
* Female: 7426
* Male: 7349

**Age categories**

| Sex    | Infant | Toddler | Child | Teenage | Young Adult | Middle-Age Adult | Older Adult | Oldest Adult | Total |
|--------|--------|----------|--------|----------|--------------|------------------|-------------|--------------|--------|
| Female | 351    | 46       | 314    | 543      | 1329         | 1664             | 768         | 91           | 5106   |
| Male   | 348    | 71       | 383    | 579      | 1135         | 1642             | 998         | 118          | 5274   |
| **Total** | 699 | 117      | 697    | 1122     | 2464         | 3306             | 1766        | 209          | 10380  |

---

## 2. Ananlyses of ARGlog10, sex and age

### 2.1 New boxplot of ARGlog10, sex and age categories

```r
Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$precise_age_category[Subset$precise_age_category == "" | Subset$precise_age_category == "NA"] <- NA

Subset$sex <- recode(Subset$sex,
                     "female" = "Female",
                     "male"   = "Male")

age_levels <- c(
  "Infant", "Toddler", "Child", "Teenage",
  "Young adult", "Middle-Age Adult",
  "Older Adult", "Oldest Adult"
)

age_labels <- c(
  "Infant" = "Infant\n(<1)",
  "Toddler" = "Toddler\n(1 - 2)",
  "Child" = "Children\n(3 - 11)",
  "Teenage" = "Teenagers\n(12 - 19)",
  "Young adult" = "Young Adult\n(20 - 34)",
  "Middle-Age Adult" = "Middle-Aged Adult\n(35 - 64)",
  "Older Adult" = "Older Adult\n(65 - 80)",
  "Oldest Adult" = "Oldest Adult\n(80<)"
)

plot_df <- Subset %>%
  filter(
    !is.na(log10_ARG_load),
    !is.na(sex),
    !is.na(precise_age_category)
  ) %>%
  mutate(
    precise_age_category = factor(precise_age_category,
                                  levels = age_levels),
    sex = factor(sex, levels = c("Female", "Male"))
  )

plot_df <- plot_df %>%
  filter(precise_age_category != "NA")


plot_df <- plot_df %>%
  mutate(
    precise_age_category = factor(precise_age_category,
                                  levels = age_levels),
    sex = factor(sex, levels = c("Female", "Male"))
  )

counts <- plot_df %>%
  group_by(precise_age_category, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(y_pos = 4.75)  # same level for all counts

ggplot(plot_df, aes(x = precise_age_category, y = log10_ARG_load, fill = sex)) +
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
  scale_x_discrete(labels = age_labels) +  # <-- use descriptive labels here
  scale_fill_npg() +
  scale_color_npg() +
  labs(
    title = "ARG Load by Age Category and Sex",
    x = "Age Category",
    y = expression(log[10]*"(ARG load)"),
    fill = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggsave("Boxplot_log10ARG_by_sex_age_ready.png", width = 8, height = 6, dpi = 300)

```
![ARG Load by Sex and age categories](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_Load_Analyses/Boxplot_log10ARG_by_sex_age_category_ready1.png)

### 2.2 Loess curve of ARG load by sex and age categories
```r

Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$precise_age_category[Subset$precise_age_category == "" |
                            Subset$precise_age_category == "NA"] <- NA

Subset$sex <- recode(Subset$sex,
                     "female" = "Female",
                     "male"   = "Male")

age_levels <- c(
  "Infant", "Toddler", "Child", "Teenage",
  "Young adult", "Middle-Age Adult",
  "Older Adult", "Oldest Adult"
)

plot_df <- Subset %>%
  filter(!is.na(log10_ARG_load),
         !is.na(sex),
         !is.na(precise_age_category),
         precise_age_category != "NA") %>%
  mutate(
    precise_age_category = factor(precise_age_category, levels = age_levels),
    age_num = as.numeric(precise_age_category),
    sex = factor(sex, levels = c("Female", "Male"))
  )

# LOESS curves only
ggplot(plot_df, aes(x = age_num, y = log10_ARG_load, color = sex, fill = sex)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.7, alpha = 0.2, size = 1.2) +
  scale_x_continuous(breaks = 1:length(age_levels), labels = age_levels) +
  scale_color_npg() +
  scale_fill_npg() +
  labs(
    title = "LOESS Curve of ARG Load by Age and Sex",
    x = "Age Category",
    y = expression(log[10]*"(ARG load)"),
    color = "Sex", fill = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("Loess_log10ARG_by_sex_age_category_ready.png", width = 8, height = 6, dpi = 300)

```
![ARG Load by Sex and age categories](https://github.com/Karhusa/F_AMR_project/blob/main/Results/AGR_Load_Analyses/LOESS_log10ARG_by_sex_age_category_ready.png)


### 2.3 Loess curve of ARG load by sex and numeric age

```
Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$age_years[Subset$age_years == "" | is.na(Subset$age_years)] <- NA

Subset$sex <- recode(Subset$sex,
                     "female" = "Female",
                     "male"   = "Male")

plot_df <- Subset %>%
  filter(!is.na(log10_ARG_load),
         !is.na(sex),
         !is.na(age_years)) %>%
  mutate(sex = factor(sex, levels = c("Female", "Male")))


# LOESS curves with numeric age
ggplot(plot_df, aes(x = age_years, y = log10_ARG_load, color = sex, fill = sex)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.7, alpha = 0.2, size = 1.2) +
  scale_color_npg() +
  scale_fill_npg() +
  labs(
    title = "LOESS Curve of ARG Load by Age (Years) and Sex",
    x = "Age (years)",
    y = expression(log[10]*"(ARG load)"),
    color = "Sex",
    fill = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave("Loess_log10ARG_by_sex_age_numeric_ready.png", width = 8, height = 6, dpi = 300)

```
![LOESS ARG Load by Sex and age categories](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_Load_Analyses/Loess_log10ARG_by_sex_age_numeric_ready(2).png)


### 2.4 Linear model
```
plot_df <- Subset %>%
  filter(!is.na(log10_ARG_load),
         !is.na(sex),
         !is.na(age_years)) %>%
  mutate(sex = factor(sex, levels = c("Female", "Male")))
Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$age_years[Subset$age_years == "" | is.na(Subset$age_years)] <- NA

lm_full <- lm(log10_ARG_load ~ age_years + sex, data = plot_df)
summary(lm_full)

```

| Min      | 1Q       | Median  | 3Q      | Max     |
| -------- | -------- | ------- | ------- | ------- |
| -1.00348 | -0.19450 | 0.01204 | 0.19203 | 1.73418 |

Coefficients

| Term | Estimate | Std. Error | t value | Pr(>t) | Significance |
|------------|------------|------------|---------|----------|--------------|
| (Intercept) | 2.6921136 | 0.0064172 | 419.515 | < 2e-16 | *** |
| age_years | 0.0014797 | 0.0001257 | 11.769 | < 2e-16 | *** |
| sexMale | -0.0084758 | 0.0060318 | -1.405 | 0.16 | |

| Item                    | Value                              |
| ----------------------- | ---------------------------------- |
| Formula                 | `log10_ARG_load ~ age_years + sex` |
| Residual Standard Error | 0.3025                             |
| Degrees of Freedom      | 10066                              |
| Multiple R²             | 0.01368                            |
| Adjusted R²             | 0.01348                            |
| F-statistic             | 69.8 (df = 2, 10066)               |
| Model p-value           | < 2.2e-16                          |




### 2.5 Optional GAM model

```r
library(mgcv)

Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$age_years[Subset$age_years == "" | is.na(Subset$age_years)] <- NA

Subset$sex <- recode(Subset$sex,
                     "female" = "Female",
                     "male"   = "Male")

plot_df <- Subset %>%
  filter(!is.na(log10_ARG_load),
         !is.na(sex),
         !is.na(age_years)) %>%
  mutate(sex = factor(sex, levels = c("Female", "Male")))

gam_model <- gam(
  log10_ARG_load ~ sex + s(age_years, by = sex),
  data = plot_df,
  method = "REML"
)

summary(gam_model)
```

Results:

| Section | Term | Estimate / Value | Std. Error | Statistic | p-value | Signif. |
|--------|------|------------------|------------|-----------|---------|---------|
| Model information | Family | Gaussian | – | – | – | – |
| Model information | Link function | Identity | – | – | – | – |
| Model information | Formula | log10_ARG_load ~ sex + s(age_years, by = sex) | – | – | – | – |
| Model information | Sample size (n) | 10069 | – | – | – | – |
| Parametric coefficients | (Intercept) | 2.751428 | 0.004256 | t = 646.462 | <2e-16 | *** |
| Parametric coefficients | sexMale | -0.013276 | 0.005966 | t = -2.225 | 0.0261 | * |
| Smooth terms | s(age_years):sexFemale | – | – | F = 21.80 (edf = 8.301, Ref.df = 8.839) | <2e-16 | *** |
| Smooth terms | s(age_years):sexMale | – | – | F = 34.52 (edf = 6.733, Ref.df = 7.793) | <2e-16 | *** |
| Model fit | Adjusted R-squared | 0.0441 | – | – | – | – |
| Model fit | Deviance explained | 4.56% | – | – | – | – |
| Model fit | −REML | 2127.9 | – | – | – | – |
| Model fit | Scale estimate | 0.088636 | – | – | – | – |
| Significance codes | *** | p < 0.001 | – | – | – | – |
| Significance codes | ** | p < 0.01 | – | – | – | – |
| Significance codes | * | p < 0.05 | – | – | – | – |
| Significance codes | . | p < 0.1 | – | – | – | – |


### 2.6 GAM plot
```r

ggplot(plot_df, aes(x = age_years, y = log10_ARG_load)) +
  geom_point(alpha = 0.2, color = "grey70") +
  geom_ribbon(
    data = newdat,
    aes(x = age_years, ymin = lwr, ymax = upr, fill = sex),
    alpha = 0.2,
    color = NA,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = newdat,
    aes(x = age_years, y = fit, color = sex),
    size = 1.5,
    inherit.aes = FALSE
  ) +
  scale_color_npg() +
  scale_fill_npg() +
  labs(
    x = "Age (years)",
    y = expression(log[10]*"(ARG load)"),
    title = "ARG load vs Age by Sex (GAM)",
    color = "Sex",
    fill = "Sex"
  ) +
  
  theme_minimal()

ggsave("GAM_log10ARG_by_sex_age_ready.png", width = 8, height = 6, dpi = 300)
```

![GAM ARG Load by Sex and numeric age](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_Load_Analyses/GAM_log10ARG_by_sex_age_ready.png)

## 3. Analyses from age categories excluding "Oldest Adults"

### 3.1  Base filtering (age ≤ 80)
```r
df_no_80 <- Subset %>%
  filter(!is.na(age_years), age_years <= 80) %>%
  mutate(
    sex = recode(sex, "female" = "Female", "male" = "Male"),
    sex = factor(sex, levels = c("Female", "Male"))
  )
```
#  3.2 Plot 1: Boxplot by precise_age_category
```
age_levels <- c("Infant", "Toddler", "Child", "Teenage",
                "Young adult", "Middle-Age Adult", "Older Adult")

age_labels <- c(
  "Infant" = "Infant\n(<1)",
  "Toddler" = "Toddler\n(1 - 2)",
  "Child" = "Children\n(3 - 11)",
  "Teenage" = "Teenagers\n(12 - 19)",
  "Young adult" = "Young Adult\n(20 - 34)",
  "Middle-Age Adult" = "Middle-Aged Adult\n(35 - 64)",
  "Older Adult" = "Older Adult\n(65 - 80)"
)

plot_df_cat <- df_no_80 %>%
  filter(
    !is.na(log10_ARG_load), 
    !is.na(sex),                          # <-- remove rows where sex is NA
    precise_age_category %in% age_levels  # <-- keep only valid categories
  ) %>%
  mutate(precise_age_category = factor(precise_age_category, levels = age_levels))

# Counts for annotation
counts <- plot_df_cat %>%
  group_by(precise_age_category, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(y_pos = 4.75)

# Boxplot
ggplot(plot_df_cat, aes(x = precise_age_category, y = log10_ARG_load, fill = sex)) +
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
  scale_x_discrete(labels = age_labels) +
  scale_fill_npg() +
  scale_color_npg() +
  labs(
    title = "ARG Load by Age Category and Sex",
    x = "Age Category",
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
### 3.3 Plot 2: Loess curve by age_years
```r
plot_df_cont <- df_no_80 %>%
  filter(!is.na(log10_ARG_load))

ggplot(plot_df_cont, aes(x = age_years, y = log10_ARG_load, color = sex)) +
  geom_point(alpha = 0.10, size = 1.0) +
  geom_smooth(method = "loess", se = TRUE, span = 0.8, size = 1.2) +
  scale_color_npg() +
  labs(
    title = "Global ARG Load by Age",
    x = "Age (years)",
    y = expression(log[10]*"(ARG load)"),
    color = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "plain")
  )
```
### 3.4 Statistical results
```
summary_table <- plot_df_cont %>%
  group_by(sex) %>%
  summarise(
    N = n(),
    mean_ARG = mean(log10_ARG_load, na.rm = TRUE),
    median_ARG = median(log10_ARG_load, na.rm = TRUE),
    sd_ARG = sd(log10_ARG_load, na.rm = TRUE),
    min_ARG = min(log10_ARG_load, na.rm = TRUE),
    max_ARG = max(log10_ARG_load, na.rm = TRUE),
    .groups = "drop"
  )

# Wilcoxon test
wilcox_p <- wilcox.test(log10_ARG_load ~ sex, data = plot_df_cont)$p.value
summary_table <- summary_table %>% mutate(wilcox_p = wilcox_p)
summary_table

```
# 3.5 GAM: ARG load ~ age_years + sex

```r
library(mgcv)
gam_model <- gam(
  log10_ARG_load ~ sex + s(age_years, by = sex),
  data = plot_df_cont,
  method = "REML"
)
summary(gam_model)
```
# 3.6 Table for statistics
```
poster_table <- summary_table %>%
  mutate(
    mean_ARG = round(mean_ARG, 2),
    median_ARG = round(median_ARG, 2),
    sd_ARG = round(sd_ARG, 2),
    min_ARG = round(min_ARG, 2),
    max_ARG = round(max_ARG, 2),
    wilcox_p = signif(wilcox_p, 3)
  ) %>%
  select(
    Sex = sex, N, Mean = mean_ARG, Median = median_ARG,
    SD = sd_ARG, Min = min_ARG, Max = max_ARG, `Wilcoxon p` = wilcox_p
  ) %>%
  gt() %>%
  tab_header(
    title = md("**ARG Load Statistics by Sex**")
  ) %>%
  arg_number(
    columns = vars(Mean, Median, SD, Min, Max),
    decimals = 2
  ) %>%
  arg_number(
    columns = vars(`Wilcoxon p`),
    decimals = 3
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(columns = everything())
  )

poster_table
```


## 4. Separate by LMIC and HIC and remove ages under 18 and over 80

### 4.1. Prepare the dataset

```r

df_income_age_ARG <- Subset %>%
  select(
    sex,
    World_Bank_Income_Group,        
    age_years,
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
    ),
    sex = recode(sex, "female" = "Female", "male" = "Male"),
    sex = factor(sex, levels = c("Female", "Male"))
  ) %>%
  filter(
    !is.na(age_years),
    age_years >= 18,    # keep adults
    age_years <= 80,    # exclude very old ages
    !is.na(sex),
    !is.na(log10_ARG_load),
    !is.na(Income_group)
  )
```
# 4.2 Plot by Income Group
```

ggplot(df_income_age_ARG, aes(x = age_years, y = log10_ARG_load, color = sex)) +
  geom_point(alpha = 0.1, size = 1) +
  geom_smooth(method = "loess", se = TRUE, span = 0.8, size = 1.2) +
  scale_color_npg() +
  labs(
    title = "Global ARG Load by Age and Income Group",
    x = "Age (years)",
    y = expression(log[10]*"(ARG load)"),
    color = "Sex"
  ) +
  facet_wrap(~Income_group) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "plain")
  )
```
# 4.3 Plot by HIC-countries
```r

ggplot(
  df_income_age_ARG %>% filter(Income_group == "HIC"),
  aes(x = age_years, y = log10_ARG_load, color = sex)
) +
  geom_point(alpha = 0.1, size = 1) +
  geom_smooth(method = "loess", se = TRUE, span = 0.8, size = 1.2) +
  scale_color_npg() +
  labs(
    title = "ARG Load by Age in HIC Countries",
    x = "Age (years)",
    y = expression(log[10]*"(ARG load)"),
    color = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "plain")
  )

```

## 4.4 Convert age into logarithmic scale
```r
df_income_age_ARG_f <- df_income_age_ARG %>%
  mutate(log10_age = log10(age_years)) %>%      # create log10_age
  filter(is.finite(log10_age))                  # remove Inf or -Inf values 
```

## 4.5 Linear models for Log10 Ages
```r
model_HIC <- lm(
  log10_ARG_load ~ log10_age + sex,
  data = df_income_age_ARG_f %>% filter(Income_group == "HIC")
)

summary(model_HIC)

model_LMIC <- lm(
  log10_ARG_load ~ log10_age + sex,
  data = df_income_age_ARG_f %>% filter(Income_group == "LMIC")
)

summary(model_LMIC)

```

## 4.6 Summary and GT table by Income Group and Sex

```
sample_counts <- df_income_age_ARG_f %>%
  group_by(Income_group, sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = sex, values_from = n) %>%
  mutate(Total = Female + Male)

# Create gt table
sample_counts %>%
  gt() %>%
  tab_header(
    title = "Sample Counts by Income Group and Sex"
  ) %>%
  cols_label(
    Income_group = "Income Group",
    Female = "Female",
    Male = "Male",
    Total = "Total"
  ) %>%
  fmt_number(columns = vars(Female, Male, Total), decimals = 0)


```

### LM model results

```
lm_LMIC <- lm(log10_ARG_load ~ log10_age + sex, 
              data = df_income_age_ARG_f %>% filter(Income_group == "LMIC"))
lm_HIC <- lm(log10_ARG_load ~ log10_age + sex, 
             data = df_income_age_ARG_f %>% filter(Income_group == "HIC"))

# Tidy model outputs
tidy_LMIC <- broom::tidy(lm_LMIC) %>%
  mutate(Income_group = "LMIC") %>%
  select(Income_group, term, estimate, std.error, statistic, p.value)

tidy_HIC <- broom::tidy(lm_HIC) %>%
  mutate(Income_group = "HIC") %>%
  select(Income_group, term, estimate, std.error, statistic, p.value)

# Combine
lm_table <- bind_rows(tidy_LMIC, tidy_HIC)

# Add model-level stats
model_stats <- tibble(
  Income_group = c("LMIC", "HIC"),
  Adj_R2 = c(summary(lm_LMIC)$adj.r.squared, summary(lm_HIC)$adj.r.squared),
  RSE = c(summary(lm_LMIC)$sigma, summary(lm_HIC)$sigma)
)

# Merge
lm_table <- lm_table %>%
  left_join(model_stats, by = "Income_group")

# Create gt table
lm_table %>%
  gt() %>%
  tab_header(title = "Linear Model Results for ARG Load") %>%
  cols_label(
    term = "Predictor",
    estimate = "β",
    std.error = "SE",
    statistic = "t",
    p.value = "p-value",
    Adj_R2 = "Adjusted R²",
    RSE = "Residual SE"
  ) %>%
  fmt_number(columns = c(estimate, std.error, statistic, Adj_R2, RSE), decimals = 2) %>%
  fmt_scientific(columns = p.value, decimals = 1)

```
### Reproductive age

```r

df_income_r_age_ARG <- Subset %>%
  select(
    sex,
    World_Bank_Income_Group,        
    age_years,
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
    ),
    sex = recode(sex, "female" = "Female", "male" = "Male"),
    sex = factor(sex, levels = c("Female", "Male"))
  ) %>%
  filter(
    !is.na(age_years),
    !is.na(sex),
    !is.na(log10_ARG_load),
    !is.na(Income_group),
    (sex == "Female" & age_years >= 15 & age_years <= 49) |
    (sex == "Male" & age_years >= 15 & age_years <= 49)
  )

```


```

ggplot(df_income_r_age_ARG, aes(x = Income_group, y = log10_ARG_load, fill = sex)) +
  geom_jitter(aes(color = sex), width = 0.2, alpha = 0.3, size = 1.5, show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  scale_fill_npg() +
  scale_color_npg() +
  labs(
    x = "Income Group",
    y = "log10(ARG load)",
    fill = "Sex"
  ) +
  ggtitle("ARG Load by Income Group and Sex\nReproductive ages: Female 15–49, Male 15–49") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

```
Models:

```
library(dplyr)

model_HIC <- lm(
  log10_ARG_load ~ age_years + sex,
  data = df_income_r_age_ARG %>% filter(Income_group == "HIC")
)

summary(model_HIC)

model_LMIC <- lm(
  log10_ARG_load ~ age_years + sex,
  data = df_income_r_age_ARG %>% filter(Income_group == "LMIC")
)

summary(model_LMIC)
```



tables:

```

sample_counts <- df_income_r_age_ARG %>%
  group_by(sex, Income_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Income_group, values_from = n, values_fill = 0)

sample_counts %>%
  gt() %>%
  tab_header(
    title = "Sample Counts by Sex and Income Group (Reproductive Age)"
  ) %>%
  cols_label(
    sex = "Sex",
    HIC = "HIC",
    LMIC = "LMIC"
  ) %>%
  fmt_number(
    columns = vars(HIC, LMIC),
    decimals = 0
  )

```

```
tidy_LMIC <- broom::tidy(lm_LMIC) %>%
  mutate(Income_group = "LMIC") %>%
  select(Income_group, term, estimate, std.error, statistic, p.value)

tidy_HIC <- broom::tidy(lm_HIC) %>%
  mutate(Income_group = "HIC") %>%
  select(Income_group, term, estimate, std.error, statistic, p.value)

# Combine
lm_table <- bind_rows(tidy_HIC, tidy_LMIC)

# Add model-level stats
model_stats <- tibble(
  Income_group = c("HIC", "LMIC"),
  Adj_R2 = c(summary(lm_HIC)$adj.r.squared, summary(lm_LMIC)$adj.r.squared),
  RSE = c(summary(lm_HIC)$sigma, summary(lm_LMIC)$sigma)
)

# Merge
lm_table <- lm_table %>%
  left_join(model_stats, by = "Income_group")
lm_table %>%
  gt() %>%
  tab_header(title = "Linear Model Results for ARG Load") %>%
  cols_label(
    term = "Predictor",
    estimate = "β",
    std.error = "SE",
    statistic = "t",
    p.value = "p-value",
    Adj_R2 = "Adjusted R²",
    RSE = "Residual SE"
  ) %>%
  fmt_number(columns = c(estimate, std.error, statistic, Adj_R2, RSE), decimals = 3) %>%
  fmt_scientific(columns = p.value, decimals = 2)
