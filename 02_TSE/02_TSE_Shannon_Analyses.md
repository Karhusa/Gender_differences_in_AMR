# ARG Shannon Analysis by Sex

---
## 1. Packages
```r
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
```
## 2. Load TSE

```r
tse_path <- "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/TSE.rds"
TSE <- readRDS(tse_path)
```
## 3. Extract colData

```r
colData_df <- as.data.frame(colData(TSE))

colData_new_subset <- colData_df %>%
  select(
    acc,
    ARG_div_shan,
    age_years,
    log10_ARG_load, 
    sex, 
    BMI_range_new, 
    Antibiotics_used, 
    UTI_history, 
    precise_age_category, imprecise_age_category
  )

```

## 4. Ananlyses of ARG Shannon diversity index and sex

### 4.1.  Boxplot of ARG Shannon diversity index and sex

```r
colData_new_subset$sex[colData_new_subset$sex == "" | colData_new_subset$sex == "NA"] <- NA  # convert empty/NA strings to actual NA
plot_df <- colData_new_subset %>% filter(!is.na(ARG_div_shan), !is.na(sex))

n_df <- plot_df %>%
  group_by(sex) %>%
  summarise(
    n = n(),
    y_max = max(ARG_div_shan, na.rm = TRUE)
  )

ggplot(plot_df, aes(x = sex, y = ARG_div_shan, fill = sex)) +
  geom_jitter(
    width = 0.15,
    size = 1.2,
    alpha = 0.25,
    color = "grey30"
  ) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.8,
    color = "grey25"
  ) +
  geom_text(
    data = n_df,
    aes(
      x = sex,
      y = y_max + 0.05,
      label = paste0("N = ", n)
    ),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = c("#B3A2C7", "#A6D854")) +
  labs(
    title = "ARG Shannon diversity by Sex",
    x = "Sex",
    y = "ARG Shannon diversity"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

ggsave("Boxplot_Shannon_diversity_by_sex.png", width = 8, height = 6, dpi = 300)
```

![Boxplot of ARG Shannon Diversity by Sex](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Boxplot_Shannon_diversity_by_sex.png)


### 4.2. Descriptive statistics (Shannon x Sex)

```r
colData_new_subset %>%
  filter(!is.na(sex), !is.na(ARG_div_shan)) %>% group_by(sex) %>%
  summarise(mean_ARGShannon = mean(ARG_div_shan), median_ARGShannon = median(ARG_div_shan), sd_ARGShannon = sd(ARG_div_shan), n = n())
```


| sex    | Mean Shannon | Median Shannon | sd_ARGShannon | N    |
|--------|----------------|-----------------|---------------|------|
| female | 1.90           | 1.94            | 0.515         | 7426 |
| male   | 1.87           | 1.93            | 0.532         | 7349 |

### 4.3. Wilcoxon test

```r
wilcox_res <- wilcox.test(ARG_div_shan ~ sex, data = plot_df)

p_label <- paste0("Wilcoxon p = ", signif(wilcox_res$p.value, 3))
```

### 4.4 Same Boxplot with wilcoxon test value

```r
ggplot(plot_df, aes(x = sex, y = ARG_div_shan, fill = sex)) +
  geom_jitter(
    width = 0.15,
    size = 1.2,
    alpha = 0.25,
    color = "grey30"
  ) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.8,
    color = "grey25"
  ) +
  geom_text(
    data = n_df,
    aes(
      x = sex,
      y = y_max + 0.05,
      label = paste0("N = ", n)
    ),
    inherit.aes = FALSE,
    size = 4
  ) +
  annotate("text", x = 1.5, y = max(plot_df$ARG_div_shan) + 0.35, label = p_label, size = 4.2, fontface = "italic") +
  scale_fill_manual(values = c("#B3A2C7", "#A6D854")) +
  labs(
    title = "ARG Shannon diversity by Sex",
    x = "Sex",
    y = "ARG Shannon diversity"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

ggsave("wilcoxon_Boxplot_Shannon_diversity_by_sex.png", width = 8, height = 6, dpi = 300)
```

![Boxplot of ARG Shannon Diversity by Sex](https://github.com/Karhusa/F_AMR_project/blob/main/Results/wilcoxon_Boxplot_Shannon_diversity_by_sex.png)

## 5. Analyses of ARG Shannon index by Age and Sex

### 5.1 Boxplot of ARG Shannon index by Category and Sex

```r
colData_new_subset$sex <- factor(colData_new_subset$sex, levels = c("female", "male"))  # make it a factor

counts <- colData_new_subset %>%
  group_by(precise_age_category, sex) %>%
  summarise(
    N = n(),
    y_pos = max(ARG_div_shan, na.rm = TRUE) + 0.3,
    .groups = "drop"
  )

age_levels <- c(
  "Infant", "Toddler", "Child", "Teenage", 
  "Young adult", "Middle-Age Adult", "Older Adult", "Oldest Adult"
)

colData_new_subset <- colData_new_subset %>% mutate(precise_age_category = factor(precise_age_category, levels = age_levels))

ggplot(colData_new_subset, aes(x = precise_age_category, y = ARG_div_shan, fill = sex)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("female" = "#B3A2C7", "male" = "#A6D854")) +
  geom_text(
    data = counts,
    aes(x = precise_age_category, y = y_pos, label = paste0("N=", N)),
    position = position_dodge(width = 0.8),
    size = 3
  ) +
  labs(
    x = "Age Category",
    y = "ARG Shannon diversity",
    title = "ARG Shannon Diversity by Age Category and Sex",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

setwd("/scratch/project_2008149/USER_WORKSPACES/karhula/DATA")

ggsave("Boxplot_Shannon_diversity_by_age_sex.png", width = 8, height = 6, dpi = 300)
```
![Boxplot of ARG Shannon Diversity by Age and Sex](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Boxplot_Shannon_diversity_by_age_sex.png)

* N values can be removed. I left N values there so that it would be easier to interpret results.

```r
age_levels <- c(
  "Infant", "Toddler", "Child", "Teenage", 
  "Young adult", "Middle-Age Adult", "Older Adult", "Oldest Adult"
)

colData_new_subset <- colData_new_subset %>%
  mutate(
    precise_age_category = factor(precise_age_category, levels = age_levels),
    sex = factor(sex, levels = c("female", "male"))
  ) %>%
  filter(!is.na(precise_age_category)) %>%  
  droplevels()

counts <- colData_new_subset %>%
  group_by(precise_age_category, sex) %>%
  summarise(
    N = n(),
    y_pos = max(ARG_div_shan, na.rm = TRUE) + 0.3,
    .groups = "drop"
  )

ggplot(colData_new_subset, aes(x = precise_age_category, y = ARG_div_shan, fill = sex)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("female" = "#B3A2C7", "male" = "#A6D854")) +
  geom_text(
    data = counts,
    aes(x = precise_age_category, y = y_pos, label = paste0("N=", N)),
    position = position_dodge(width = 0.8),
    size = 3
  ) +
  labs(
    x = "Age Category",
    y = "ARG Shannon diversity",
    title = "ARG Shannon Diversity by Age Category and Sex",
    fill = "Sex"
  ) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

setwd("/scratch/project_2008149/USER_WORKSPACES/karhula/DATA")

ggsave("Boxplot_Shannon_diversity_by_age_sex_no_NA.png", width = 8, height = 6, dpi = 300)

```

![Boxplot of ARG Shannon Diversity by Age and Sex](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Boxplot_Shannon_diversity_by_age_sex_no_NA.png)



### 5.2 Scatter plot + separate regression lines of ARG Shannon index by sex and age

```r
colData_sex_clean <- colData_new_subset %>% filter(!is.na(sex))
N_sex <- colData_sex_clean %>%
  filter(!is.na(age_years), !is.na(ARG_div_shan)) %>%
  count(sex)

ggplot(colData_sex_clean,
       aes(x = age_years, y = ARG_div_shan, color = sex)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.4, alpha = 0.15
  ) +
  scale_color_manual(
    values = c("female" = "#B3A2C7", "male" = "#A6D854")
  ) +
  labs(x = "Age (years)", y = "ARG Shannon diversity", title = "ARG Shannon Diversity by Age and Sex", color = "Sex"
  ) +
  theme_minimal()

ggsave("Regression_Shannon_diversity_by_age_sex.png", width = 8, height = 6, dpi = 300)

```

![Regression analysis ARG Shannon diversity by Age and Sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/Regression_Shannon_diversity_by_age_sex.png)

### 5.3 Linear regression with results
```
model <- lm(ARG_div_shan ~ age_years + sex, data = colData_sex_clean)
summary(model)

model_sum <- summary(model)

beta_age <- model_sum$coefficients["age_years", "Estimate"]
p_age    <- model_sum$coefficients["age_years", "Pr(>|t|)"]
r2 <- model_sum$r.squared
n  <- model_sum$df[1] + model_sum$df[2] 

ggplot(colData_sex_clean, aes(x = age_years, y = ARG_div_shan)
  ) +
  geom_point(alpha = 0.2, color = "grey60") +
  geom_smooth(aes(color = sex), method = "lm", se = TRUE,linewidth = 1.5
  ) +
  scale_color_manual(
    values = c("female" = "#B3A2C7", "male" = "#A6D854")
  ) +
  annotate("text", x = Inf, y = Inf,
    hjust = 1.1, vjust = 1.2,
    size = 4,
    label = paste0(
      "Age effect: β = ", round(beta_age, 4),
      "\nP = ", signif(p_age, 3),
      "\nR² = ", round(r2, 3),
      "\nN = ", n
    )
  ) +
  labs(x = "Age (years)", y = "ARG Shannon diversity", title = "ARG Shannon Diversity by Age and Sex", color = "Sex"
  ) +
  theme_minimal()

ggsave("Regression_with_table_ARG_Shannon_by_age_sex.png", width = 8, height = 6, dpi = 300)
```

![Regression analysis with table ARG Shannon index by Age and Sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/Regression_with_table_ARG_Shannon_by_age_sex.png)

| Term       | Estimate| Std. Error | t value | Pr(>t)   | Significance|
|------------|---------|------------|---------|------------|-------------|
| (Intercept)| 1.7123  | 0.0103     | 165.93  | < 2e-16    | ***         |
| age_years  | 0.00603 | 0.00020    | 29.81   | < 2e-16    | ***         |
| sexmale    | -0.06130| 0.00970    | -6.32   | 2.73e-10   | ***         |


* Residual standard error: 0.4864 on 10066 degrees of freedom  
* Multiple R-squared: 0.0836, Adjusted R-squared: 0.0834  
* F-statistic: 459 on 2 and 10066 DF, p-value: < 2.2e-16  
* 4706 observations deleted due to missingness


### 5.4 Generalized Additive Model (GAM)
* A GAM models the outcome as a sum of smooth functions of predictors rather than simple linear effects.

```r
library(mgcv)

colData_sex_clean <- colData_new_subset %>% filter(!is.na(sex))

gam_model <- gam(
  ARG_div_shan ~ s(age_years) + sex,
  data = colData_sex_clean,
  method = "REML"
)

summary(gam_model)

ggplot(colData_sex_clean, aes(x = age_years, y = ARG_div_shan)) +
  geom_point(alpha = 0.2, color = "grey70") +
  geom_smooth(
    aes(color = sex),
    method = "gam",
    formula = y ~ s(x),
    se = TRUE,
    linewidth = 1.5
  ) +
  scale_color_manual(values = c("female" = "#B3A2C7", "male" = "#A6D854")) +
  labs(x = "Age (years)", y = "ARG Shannon", title = "ARG Shannon vs Age by Sex (GAM)", color = "Sex") +
  theme_minimal()

ggsave("GAM_ARG_Shannon_by_age_sex.png", width = 8, height = 6, dpi = 300)

```
![Gam analysis ARG Shannon by Age and Sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/GAM_ARG_Shannon_by_age_sex.png)

| Term       | Estimate | Std. Error | t value | Pr(>t) | Significance|
| ---------- |----------| ---------- | ------- | -------- | ----------- |
| (Intercept)| 1.9466   | 0.00689    | 282.43  | < 2e-16  | ***         |
| sexmale    | -0.06394 | 0.00970    | -6.595  | 4.47e-11 | ***         |


| Term          | edf    | Ref.df | F value | p-value | Significance |
|---------------|-------|--------|---------|---------|-------------|
| s(age_years)  | 7.985 | 8.688  | 115.3   | < 2e-16 | ***         |


* Family: Gaussian, link = identity  
*  R-squared (adj) = 0.093
*  Deviance explained = 9.38%  
* REML = 7000.7  
* Scale estimate = 0.23406  
* Sample size (n) = 10069


---
## 6. Analyses of ARG Shannon index by BMI and Sex

### 6.1 Boxplot of ARG Shannon index by Category and Sex


```r
colData_new_subset$sex[colData_new_subset$BMI_range_new  == "" | colData_new_subset$BMI_range_new  == "NA"] <- NA

colData_subset_clean <- colData_new_subset %>% filter(BMI_range_new != "Normal/Overweight (<30)")

colData_subset_clean$BMI_range_new <- factor(
  colData_subset_clean$BMI_range_new,
  levels = c("Underweight (<18.5)", "Normal (18.5-25)", 
             "Overweight (25-30)", "Obese (>30)")
)

counts <- colData_subset_clean %>%
  group_by(BMI_range_new, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(y_pos = max(colData_subset_clean$ARG_div_shan, na.rm = TRUE) + 0.1)

ggplot(colData_subset_clean, aes(x = BMI_range_new, y = ARG_div_shan, fill = sex)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8), na.rm = TRUE) +
  scale_fill_manual(values = c("female" = "#B3A2C7", "male" = "#A6D854")) +
  geom_text(
    data = counts,
    aes(x = BMI_range_new, y = y_pos, label = paste0("N=", N)),
    position = position_dodge(width = 0.8),
    size = 3) +
  labs(x = "BMI Category", y = "ARG Shannon Diveristy", title = "ARG Shannon by BMI Category and Sex", fill = "Sex") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Shannon_Boxplot_by_BMI_Sex.png", width = 8, height = 6, dpi = 300)

```
![Shannon Boxplot by BMI and Sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/Shannon_Boxplot_by_BMI_Sex.png)

## 6.2 Calclulations

```r
colData_subset_clean$BMI_range_new <- relevel(colData_subset_clean$BMI_range_new, ref = "Normal (18.5-25)")

model_add_norm <- lm(formula = ARG_div_shan ~ BMI_range_new + sex, data = colData_subset_clean)
summary(model_add_norm)

model_int_norm <- lm(ARG_div_shan ~ BMI_range_new * sex, data = colData_subset_clean)
summary(model_int_norm)

```
**Additive model:**

| Term                               | Estimate  | Std. Error | t value | Pr(>t) | Significance |
|------------------------------------|-----------|------------|---------|----------|--------------|
| (Intercept)                        | 2.034651  | 0.010874   | 187.105 | <2e-16   | ***          |
| BMI_range_newUnderweight (<18.5)   | 0.062869  | 0.024795   | 2.536   | 0.01125  | *            |
| BMI_range_newOverweight (25–30)    | -0.047084 | 0.016097   | -2.925  | 0.00346  | **           |
| BMI_range_newObese (>30)           | -0.003632 | 0.016825   | -0.216  | 0.82912  |              |
| sexmale                            | -0.054192 | 0.012691   | -4.270  | 1.98e-05 | ***          |

| Min     | 1Q       | Median   | 3Q       | Max     |
|---------|----------|----------|----------|---------|
| -1.8104 | -0.2862  | 0.0488   | 0.3571   | 1.5165  |

| Metric                    | Value        |
|---------------------------|--------------|
| Residual Std. Error       | 0.4984       |
| Degrees of Freedom        | 6252         |
| Observations Used         | 6257         |
| Observations Removed      | 18292        |
| Multiple R²               | 0.006591     |
| Adjusted R²               | 0.005956     |
| F-statistic               | 10.37        |
| Model p-value             | 2.27e-08     |

**Interaction model:**

| Term                                        | Estimate  | Std. Error | t value | Pr(>t) | Significance |
|---------------------------------------------|-----------|------------|---------|----------|--------------|
| (Intercept)                                 | 2.05754   | 0.01247    | 165.023 | <2e-16   | ***          |
| BMI_range_newUnderweight (<18.5)            | 0.04798   | 0.03154    | 1.521   | 0.128    |              |
| BMI_range_newOverweight (25–30)             | -0.13441  | 0.02451    | -5.483  | 4.34e-08 | ***          |
| BMI_range_newObese (>30)                    | -0.03130  | 0.02345    | -1.335  | 0.182    |              |
| sexmale                                     | -0.09972  | 0.01758    | -5.671  | 1.48e-08 | ***          |
| Underweight (<18.5) × male                  | 0.02365   | 0.05104    | 0.463   | 0.643    |              |
| Overweight (25–30) × male                   | 0.15349   | 0.03249    | 4.725   | 2.36e-06 | ***          |
| Obese (>30) × male                          | 0.05546   | 0.03361    | 1.650   | 0.099    | .            |

| Min     | 1Q       | Median   | 3Q       | Max     |
|---------|----------|----------|----------|---------|
| -1.7877 | -0.2846  | 0.0499   | 0.3536   | 1.5391  |

| Metric                    | Value        |
|---------------------------|--------------|
| Residual Std. Error       | 0.4976       |
| Degrees of Freedom        | 6249         |
| Observations Used         | 6257         |
| Observations Removed      | 18292        |
| Multiple R²               | 0.01018      |
| Adjusted R²               | 0.009071     |
| F-statistic               | 9.181        |
| Model p-value             | 2.43e-11     |



### 6.3 Interaction model plot

```r
colData_subset_clean$BMI_range_new <- factor(colData_subset_clean$BMI_range_new,
  levels = c("Underweight (<18.5)", "Normal (18.5-25)", "Overweight (25-30)", "Obese (>30)"))

counts <- colData_subset_clean %>%
  group_by(BMI_range_new, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(
    y_pos = max(colData_subset_clean$ARG_div_shan, na.rm = TRUE) + 0.05,
    y_offset = ifelse(sex == "female", 0.02, -0.02)  # separate male/female N labels slightly
  )

line_colors <- c("female" = "#B3A2C7", "male" = "#A6D854")  

# Plot Interaction Model
ggplot(colData_subset_clean, aes(x = BMI_range_new, y = ARG_div_shan, color = sex, group = sex)) +
  geom_jitter(width = 0.2, alpha = 0.3, color = "grey60") +  # points in light grey
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, aes(color = sex)) +
  scale_color_manual(values = line_colors) +
  geom_text(
    data = counts,
    aes(x = BMI_range_new, y = y_pos + y_offset, label = paste0("N=", N), color = sex),
    position = position_dodge(width = 0.8),
    size = 3,
    show.legend = FALSE
  ) +
  labs(
    x = "BMI Category",
    y = "ARG Shannon diversity",
    title = "ARG Shannon diversity by BMI × Sex",
    color = "Sex"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Interaction_model_Shannon_ by_BMI_sex.png", width = 8, height = 6, dpi = 300)
```

![Interaction model Shannon by BMI sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/Interaction_model_Shannon_ by_BMI_sex(1).png)

---

## 7. Analyses of ARG Load by UTI and Sex

### 7.1 Boxplot of UTI and sex
```r

colData_subset_clean <- colData_new_subset %>%
  filter(!is.na(UTI_history), !is.na(sex)) %>%
  mutate(
    UTI_history = factor(as.character(UTI_history),
                         levels = c("No", "Yes")),
    sex = factor(as.character(sex),
                 levels = c("female", "male"))
  )

# 2. Counts for annotation
counts_uti <- colData_subset_clean %>%
  group_by(UTI_history, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(
    y_pos = max(colData_subset_clean$ARG_div_shan, na.rm = TRUE) + 0.1
  )

# 3. Plot
ggplot(colData_subset_clean,
       aes(x = UTI_history, y = ARG_div_shan, fill = sex)) +
  geom_boxplot(
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  geom_text(
    data = counts_uti,
    aes(
      x = UTI_history,
      y = y_pos,
      label = paste0("N=", N)
    ),
    position = position_dodge(width = 0.8),
    size = 3
  ) +
  scale_fill_manual(
    values = c("female" = "#B3A2C7", "male" = "#A6D854")
  ) +
  labs(
    x = "UTI History",
    y = "ARG Shannon diversity",
    title = "ARG Shannon diversity by UTI History and Sex",
    fill = "Sex"
  ) +
  theme_minimal()

ggsave("Shannon_Boxplot_by_UTI_Sex.png", width = 7, height = 5, dpi = 300)

```

![Boxplot ARG Load by UTI_history and Sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/Shannon_Boxplot_by_UTI_Sex.png)

### 7.2 Linear model 

```r
lm_uti <- lm(ARG_div_shan ~ UTI_history + sex, data = colData_subset_clean)
summary(lm_uti)
```
| Term             | Estimate   | Std. Error | t value  | Pr(>t) | Significance |
|------------------|-----------|------------|---------|----------|--------------|
| (Intercept)      | 1.903693  | 0.006114   | 311.374 | <2e-16   | ***          |
| UTI_historyYes   | -0.025704 | 0.033403   | -0.769  | 0.442    |              |
| sexmale          | -0.034153 | 0.008616   | -3.964  | 7.42e-05 | ***          |

| Min       | 1Q        | Median   | 3Q       | Max      |
|-----------|----------|----------|----------|----------|
| -1.88005  | -0.32956 | 0.04864  | 0.37961  | 1.62737  |

| Metric                     | Value        |
|----------------------------|--------------|
| Residual Std. Error        | 0.5234       |
| Degrees of Freedom         | 14772        |
| Observations Used          | 14775        |
| Observations Removed       | 9830         |
| Multiple R²                | 0.001091     |
| Adjusted R²                | 0.000956     |
| F-statistic                | 8.065        |
| Model p-value              | 0.0003157    |




### 7.3 Interactive model

```r
lm_uti_int <- lm(ARG_div_shan ~ UTI_history * sex, data = colData_subset_clean)
summary(lm_uti_int)

```

| Term                     | Estimate   | Std. Error | t value  | Pr(>t) | Significance |
|--------------------------|-----------|------------|---------|----------|--------------|
| (Intercept)              | 1.903315  | 0.006138   | 310.064 | <2e-16   | ***          |
| UTI_historyYes           | -0.007583 | 0.042489   | -0.178  | 0.858353 |              |
| sexmale                  | -0.033395 | 0.008686   | -3.845  | 0.000121 | ***          |
| UTI_historyYes:sexmale   | -0.047447 | 0.068753   | -0.690  | 0.490136 |              |

| Min       | 1Q        | Median   | 3Q       | Max      |
|-----------|----------|----------|----------|----------|
| -1.87967  | -0.32917 | 0.04817  | 0.37963  | 1.62699  |

| Metric                     | Value        |
|----------------------------|--------------|
| Residual Std. Error        | 0.5234       |
| Degrees of Freedom         | 14771        |
| Observations Used          | 14774        |
| Observations Removed       | 9830         |
| Multiple R²                | 0.001123     |
| Adjusted R²                | 0.0009201    |
| F-statistic                | 5.535        |
| Model p-value              | 0.0008552    |


---

## 8. Shannon Analyses of Antibiotic Use and Sex

### 8.1 Boxplot of Shannon diversity by Antibiotics Use and Sex

```r
counts_abx <- colData_subset_clean %>%
  group_by(Antibiotics_used, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(y_pos = max(colData_subset_clean$ARG_div_shan) + 0.1)

ggplot(colData_subset_clean, aes(x = Antibiotics_used, y = ARG_div_shan, fill = sex)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = c("female" = "#B3A2C7", "male" = "#A6D854")) +
    geom_text(
        data = counts_abx,
        aes(x = Antibiotics_used, y = y_pos, label = paste0("N=", N)),
        position = position_dodge(width = 0.8),
        size = 3
    ) +
    labs(
        x = "Antibiotics Used",
        y = "ARG Shannon diversity",
        title = "ARG Shannon diversity by Antibiotics Use and Sex",
        fill = "Sex"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


ggsave("Shannon_Boxplot_by_AB_use_Sex.png", width = 7, height = 5, dpi = 300)

```

![Boxplot ARG Load by AB Use and Sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/Shannon_Boxplot_by_AB_use_Sex.png)

### 9.2 Linear model 

```r
lm_ab <- lm(ARG_div_shan ~ Antibiotics_used + sex, data = colData_subset_clean)

summary(lm_ab)
```

| Section                 | Term                  | Estimate / Value | Std. Error | Statistic   | p-value   | Signif. |
|-------------------------|----------------------|----------------|------------|------------|-----------|---------|



```


### 9.3 Interaction model
```
lm_interaction <- lm(ARG_div_shan ~ Antibiotics_used * sex, data = colData_subset_clean)
summary(lm_interaction)
```
