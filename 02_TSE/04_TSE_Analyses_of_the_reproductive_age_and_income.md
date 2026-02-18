
# 1. Prepare the data

### 1.2 Packages
```r
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
library(mgcv)
library(emmeans)
```
---
### 1.2 Load TSE
```r
tse_path <- "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/TSE.rds"
TSE <- readRDS(tse_path)
```
### 1.3 Extract colData
```r
colData_df <- as.data.frame(colData(TSE))

Subset1 <- colData_df %>%
  select(
    acc,
    ARG_div_shan,
    age_years,
    log10_ARG_load, 
    World_Bank_Income_Group,
    sex, 
    BMI_range_new, 
    Antibiotics_used, 
    UTI_history
  ) %>%
  filter(!is.na(age_years)) %>%
  filter(age_years >= 15 & age_years <= 49) %>%
  filter(!is.na(log10_ARG_load), !is.na(sex)) %>%
  mutate(
    # Keep sex as factor
    sex = recode(sex, "female" = "Female", "male" = "Male"),
    sex = factor(sex, levels = c("Female", "Male")),
    
    # Create 5-year AgeGroup for categorical plots
    AgeGroup = cut(
      age_years,
      breaks = seq(15, 50, by = 5),
      right = FALSE,
      labels = paste(seq(15, 45, by = 5), seq(19, 49, by = 5), sep = "–")
    )
  )

```
### 1.4. Inspect the values

```r
table(Subset1$AgeGroup)

table(Subset1$sex)
```
**Results:**

AgeGroup:
* 15–19 803
* 20–24 1115
* 25–29 913
* 30–34 515
* 35–39 493
* 40–44 382
* 45–49 439

Sex:
* Female 2342
* Male 2102

---

# 2. Shannon diveristy

## 2.1 Shannon Diversity by sex and filtered age

### 2.1.1 Boxplot of Shannon Diversity by sex and filtered age

```r

Subset1$sex[Subset1$sex == "" | Subset1$sex == "NA"] <- NA
Subset1$sex[Subset1$ARG_div_shan == "" | Subset1$ARG_div_shan == "NA"] <- NA

Subset1$sex <- recode(Subset1$sex,
                     "female" = "Female",
                     "male"   = "Male")

plot_df <- Subset1 %>% filter(!is.na(sex),!is.na(ARG_div_shan))

n_df <- plot_df %>%
  group_by(sex) %>%
  summarise(
    n = n(),
    y_max = max(ARG_div_shan, na.rm = TRUE)
  )

npg_cols <- pal_npg("nrc")(4)[3:4]

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
  scale_fill_manual(values = npg_cols) +
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

ggsave("Boxplot_Shannon_diversity_by_sex_and_filtered_age.png", width = 8, height = 6, dpi = 300)

```
![Boxplot of Shannon_index with sex and filtered age](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Shannon_Analyses/Boxplot_Shannon_diversity_by_sex_and_filtered_age.png)

### 2.1.2 Do we have repeated samples?

```r
plot_df %>%
  count(acc) %>%
  summarise(
    n_subjects = n(),
    max_samples_per_subject = max(n),
    median_samples_per_subject = median(n),
    subjects_with_repeats = sum(n > 1)
  )

```
**Results:**
* n_subjects: 4444
* max_samples_per_subject: 1
* median_samples_per_subject:1 
* subjects_with_repeats :0

### 2.1.3 Wilcoxon rank-sum test
```r
wilcox.test(ARG_div_shan ~ sex, data = plot_df)
```

**Results:**
* W = 2996053
* p-value < 2.2e-16

### 2.1.4 Cound Cohens d
```r
group1 <- plot_df$ARG_div_shan[plot_df$sex == "Female"]
group2 <- plot_df$ARG_div_shan[plot_df$sex == "Male"]

# Means and SDs
mean1 <- mean(group1, na.rm = TRUE)
mean2 <- mean(group2, na.rm = TRUE)
sd1   <- sd(group1, na.rm = TRUE)
sd2   <- sd(group2, na.rm = TRUE)
n1    <- length(group1)
n2    <- length(group2)

# Pooled SD
pooled_sd <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2 - 2))

# Cohen's d
cohen_d <- (mean1 - mean2) / pooled_sd
cohen_d
```
**Results:**
* 0.4051338

### 2.1.5 Regression analysis of shannon index and sex

```r
plot_df$sex <- factor(plot_df$sex, levels = c("Female", "Male"))
lm_sex <- lm(ARG_div_shan ~ sex, data = plot_df)
summary(lm_sex)

```
**Results:**

Coefficients

| Term | Estimate | Std. Error | t value | Pr(>t) | Significance |
|-------------|---------:|-----------:|--------:|---------:|:------------|
| (Intercept) | 1.94464 | 0.01061 | 183.23 | < 2e-16 | *** |
| sexMale | -0.20809 | 0.01543 | -13.48 | < 2e-16 | *** |


* Residual standard error: 0.5136 on 4442 degrees of freedom
* Multiple R-squared:  0.03932,	Adjusted R-squared:  0.03911 
* F-statistic: 181.8 on 1 and 4442 DF,  p-value: < 2.2e-16


## 3. Ananlyses of ARGlog10 sex and income

### 3.1 Boxbot of ARGlog10 sex and income
```r
Subset1$sex[Subset1$sex == "" | Subset1$sex == "NA"] <- NA
Subset1$World_Bank_Income_Group[Subset1$World_Bank_Income_Group == "" | Subset1$World_Bank_Income_Group == "NA"] <- NA

Subset1$sex <- recode(Subset1$sex,
                     "female" = "Female",
                     "male"   = "Male")

plot_df <- Subset1 %>%
  filter(!is.na(log10_ARG_load),
         !is.na(sex),
         !is.na(World_Bank_Income_Group))

income_levels <- c("Low income", "Lower middle income", "Upper middle income", "High income")
plot_df$World_Bank_Income_Group <- factor(plot_df$World_Bank_Income_Group, levels = income_levels)

n_df <- plot_df %>%
  group_by(World_Bank_Income_Group, sex) %>%
  summarise(n = n(), .groups = "drop")

# Create combined x-axis label with ordered income groups
plot_df <- plot_df %>%
  mutate(x_group = paste0(sex, "\n(", World_Bank_Income_Group, ")"))

n_df <- n_df %>% mutate(x_group = paste0(sex, "\n(", World_Bank_Income_Group, ")"))

# Make x_group a factor with the same order as income_levels for plotting
plot_df$x_group <- factor(plot_df$x_group, levels = n_df$x_group)
n_df$x_group <- factor(n_df$x_group, levels = n_df$x_group)

# Plot
ggplot(plot_df, aes(x = x_group, y = log10_ARG_load, fill = sex)) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.25, color = "grey30") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8) +
  geom_text(data = n_df,
            aes(x = x_group,
                y = max(plot_df$log10_ARG_load, na.rm = TRUE) + 0.1,
                label = paste0("N = ", n)),
            inherit.aes = FALSE,
            size = 4) +
  scale_fill_npg() +
  labs(title = "ARG Load by Sex and Income Group",
       x = "Sex (Income Group)",
       y = expression(log[10]*"(ARG load)")) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


ggsave("Boxplot_ARG_Load_Sex_Income.png", width = 8, height = 6, dpi = 300)

```
![Boxplot of ARG load with sex and income](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_and_Income_analyses/Boxplot_ARG_Load_Sex_Income.png)

## 3.2 Regression analysis of ARGlog10 index, income and sex
```r

plot_df <- plot_df %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    World_Bank_Income_Group = factor(
      World_Bank_Income_Group,
      levels = c("Low income",
                 "Lower middle income",
                 "Upper middle income",
                 "High income")
    )
  )

# Additive mdoel
lm_main <- lm(log10_ARG_load ~ sex + World_Bank_Income_Group, data = plot_df)
summary(lm_main)

# Interaction model
lm_int <- lm(log10_ARG_load ~ sex * World_Bank_Income_Group, data = plot_df)
summary(lm_int)

```
**Results: Additive model**
* Reference group: Female, Low income

Coefficients

| Term | Estimate | Std. Error | t value | p-value | Signif. |
|------|---------:|-----------:|--------:|--------:|:-------:|
| (Intercept) | 2.777709 | 0.021045 | 131.988 | <2e-16 | *** |
| sexMale | -0.025041 | 0.005346 | -4.684 | 2.84e-06 | *** |
| World_Bank_Income_GroupLower middle income | 0.152621 | 0.026204 | 5.824 | 5.87e-09 | *** |
| World_Bank_Income_GroupUpper middle income | 0.126062 | 0.021909 | 5.754 | 8.92e-09 | *** |
| World_Bank_Income_GroupHigh income | -0.042325 | 0.021036 | -2.012 | 0.0442 | * |

Significance codes:  
*** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1

Model Fit

| Metric | Value |
|--------|-------|
| Residual standard error | 0.2982 |
| Degrees of freedom | 12710 |
| R-squared | 0.05172 |
| Adjusted R-squared | 0.05143 |
| F-statistic | 173.3 |
| Model p-value | < 2.2e-16 |


**Results: interaction model**
* Reference group: Female, Low income

Coefficients

| Term | Estimate | Std. Error | t value | p-value | Signif. |
|------|---------:|-----------:|--------:|--------:|:-------:|
| (Intercept) | 2.735527 | 0.031600 | 86.567 | <2e-16 | *** |
| sexMale | 0.049505 | 0.042008 | 1.178 | 0.2386 |  |
| World_Bank_Income_GroupLower middle income | 0.168011 | 0.038325 | 4.384 | 1.18e-05 | *** |
| World_Bank_Income_GroupUpper middle income | 0.162007 | 0.032649 | 4.962 | 7.06e-07 | *** |
| World_Bank_Income_GroupHigh income | 0.003468 | 0.031897 | 0.109 | 0.9134 |  |
| sexMale:World_Bank_Income_GroupLower middle income | -0.016688 | 0.052727 | -0.316 | 0.7516 |  |
| sexMale:World_Bank_Income_GroupUpper middle income | -0.055994 | 0.044330 | -1.263 | 0.2066 |  |
| sexMale:World_Bank_Income_GroupHigh income | -0.081284 | 0.042425 | -1.916 | 0.0554 | . |

Significance codes:  
*** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1

### Model fit

| Metric | Value |
|--------|-------|
| Residual standard error | 0.2981 |
| Degrees of freedom | 12707 |
| R-squared | 0.05243 |
| Adjusted R-squared | 0.05191 |
| F-statistic | 100.4 |
| Model p-value | < 2.2e-16 |

---

## 4. Analysis of ARG load,sex and filtered age (Women of reproductive age (15-49 years))

### 4.1 Boxplot of ARG load, sex and filtered age 

```r
Subset1$sex[Subset1$sex == "" | Subset1$sex == "NA"] <- NA
Subset1$AgeGroup[Subset1$AgeGroup == "" | Subset1$AgeGroup == "NA"] <- NA

# Recode sex
Subset1$sex <- recode(Subset1$sex,
                      "female" = "Female",
                      "male"   = "Male")
age_levels <- c(
  "15–19", "20–24", "25–29", "30–34",
  "35–39", "40–44", "45–49"
)

plot_df <- Subset1 %>%
  filter(
    !is.na(log10_ARG_load),
    !is.na(sex),
    !is.na(AgeGroup)
  ) %>%
  mutate(
    AgeGroup = factor(AgeGroup, levels = age_levels),
    sex = factor(sex, levels = c("Female", "Male"))
  )

counts <- plot_df %>%
  group_by(AgeGroup, sex) %>%
  summarise(
    N = n(),
    y_pos = max(log10_ARG_load, na.rm = TRUE) + 0.15,
    .groups = "drop"
  )

ggplot(plot_df, aes(x = AgeGroup, y = log10_ARG_load, fill = sex)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    size = 1.2, alpha = 0.25, color = "grey30"
  ) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.8,
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    data = counts,
    aes(x = AgeGroup, y = y_pos, label = N, color = sex),
    position = position_dodge(width = 0.6),
    inherit.aes = FALSE,
    size = 2.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_npg() +
  scale_color_npg() +
  labs(
    title = "ARG Load by 5-Year Age Groups (15–49) and Sex",
    x = "Age group (years)",
    y = expression(log[10]*"(ARG load)"),
    fill = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("Boxplot_log10ARG_by_reproductive_sex_age_ready.png", width = 8, height = 6, dpi = 300)
```
![Boxplot of ARG load with sex and income](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_filtered_Age_Analyses/Boxplot_log10ARG_by_reproductive_sex_age_ready.png)


### 4.2 Loess curve of ARG load, sex and filtered age 
```
Subset1$sex[Subset1$sex == "" | Subset1$sex == "NA"] <- NA
Subset1$age_years[Subset1$age_years == "" | Subset1$age_years == "NA"] <- NA

Subset1$sex <- recode(Subset1$sex,
                      "female" = "Female",
                      "male"   = "Male")

plot_df <- Subset1 %>%
  filter(
    !is.na(log10_ARG_load),
    !is.na(sex),
    !is.na(age_years)
  ) %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male"))
  )

sex_counts <- plot_df %>%
  group_by(sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(label = paste0(sex, " (N=", N, ")"))

plot_df <- plot_df %>%
  mutate(sex_label = factor(sex, levels = sex_counts$sex, labels = sex_counts$label))

# --- Plot LOESS curves using continuous age_years ---
ggplot(plot_df, aes(x = age_years, y = log10_ARG_load, color = sex_label, fill = sex_label)) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    span = 0.7,
    alpha = 0.2,
    linewidth = 1.2
  ) +
  scale_color_npg() +
  scale_fill_npg() +
  labs(
    title = "LOESS Curve of ARG Load by Age and Sex",
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


ggsave("Loess_log10ARG_by_reproductive_sex_age_ready.png", width = 8, height = 6, dpi = 300)
```
![Loess of ARG load with sex and filtered age](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_filtered_Age_Analyses/Loess_log10ARG_by_reproductive_sex_age_ready.png)

### 4.3 Linear model of of ARG load, sex and filtered age
```r
* Additive
lm_add <- lm(log10_ARG_load ~ age_years + sex, data = Subset1)
summary(lm_add)

* Interaction
lm_num <- lm(log10_ARG_load ~ age_years * sex, data = Subset1)
summary(lm_num)
```
**Additive model**
lm(formula = log10_ARG_load ~ age_years + sex, data = Subset1)

Residuals: Min: -1.0791, 1Q: -0.1905, Median: 0.0011, 3Q:0.2046, Max: 1.6629 

Coefficients:

| Term | Estimate | Std. Error | t value | Pr(>t) | Significance |
|-------------------|-----------|------------|---------|----------|--------------|
| (Intercept) | 2.5767677 | 0.0146338 | 176.084 | < 2e-16 | *** |
| age_years | 0.0055431 | 0.0004708 | 11.775 | < 2e-16 | *** |
| sexMale | -0.0604523 | 0.0088827 | -6.806 | 1.14e-11 | *** |

* Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
* Residual standard error: 0.2955 on 4441 degrees of freedom (216 observations deleted due to missingness)
* Multiple R-squared:  0.03907,	Adjusted R-squared:  0.03864 
* F-statistic: 90.29 on 2 and 4441 DF,  p-value: < 2.2e-16


**Interaction model** NOT the best for this!
lm(formula = log10(ARG load) ~ AgeNum * sex)

| Term | Estimate | Std. Error | t value | Pr(>t) | Significance |
|-------------------|-----------|------------|---------|----------|--------------|
| (Intercept) | 2.5328 | 0.01955 | 129.549 | < 2e-16 | *** |
| age_years | 0.00710 | 0.000658 | 10.798 | < 2e-16 | *** |
| sexMale | 0.03043 | 0.02826 | 1.077 | 0.282 | |
| age_years:sexMale | -0.00319 | 0.000941 | -3.387 | 0.000711 | *** |

Model statistics:
* Residual standard error: 0.2952 on 4440 df
* Multiple R-squared: 0.04155
* Adjusted R-squared: 0.0409
* F-statistic: 64.16 on 3 and 4440 DF, p-value: < 2.2e-16

### 4.4 GAM model of of ARG load, sex and filtered age
```r
gam_model <- gam(
  log10_ARG_load ~ sex + s(age_years, by = sex),
  data = plot_df,
  method = "REML"
)

summary(gam_model)

```
**Family:** gaussian  
**Link function:** identity  
**Formula:** `log10_ARG_load ~ sex + s(age_years, by = sex)`  
**Sample size:** 4,444  

Parametric Coefficients

| Term        | Estimate | Std. Error | t value | p-value       |
|------------ |---------:|-----------:|--------:|---------------|
| (Intercept) | 2.736069 | 0.006116   | 447.33  | < 2e-16 ***   |
| sexMale     | -0.058478 | 0.008915  | -6.56   | 6.01e-11 ***  |

Smooth Terms (Approximate Significance)

| Term                     | edf   | Ref.df | F value | p-value     |
|-------------------------- |------:|-------:|--------:|------------|
| s(age_years):sexFemale    | 6.307 | 7.429  | 19.61   | < 2e-16 *** |
| s(age_years):sexMale      | 8.111 | 8.762  | 12.13   | < 2e-16 *** |

Model fit:
* Adjusted R² = 0.0637
* Deviance explained = 6.69%
* REML = 866.66
* Scale est. = 0.08506

### Pairwise comparison

```r
emm <- emmeans(lm_add, ~ sex)
pairs(emm)
```

| Contrast      | Estimate | SE      | df   | t.ratio | p.value |
| ------------- | -------- | ------- | ---- | ------- | ------- |
| Female - Male | 0.0605   | 0.00888 | 4441 | 6.806   | <.0001  |

### Adjust for study
---

## 5. Analysis of ARG load, sex, income and filtered age (Women of reproductive age (15-49 years))

### 5.1 Boxplot of ARG load, sex, income, and filtered age 

```r
Subset1$sex[Subset1$sex == "" | Subset1$sex == "NA"] <- NA
Subset1$AgeGroup[Subset1$AgeGroup == "" | Subset1$AgeGroup == "NA"] <- NA
Subset1$World_Bank_Income_Group[Subset1$World_Bank_Income_Group == "" | Subset1$World_Bank_Income_Group == "NA"] <- NA

Subset1$sex <- recode(Subset1$sex,
                      "female" = "Female",
                      "male"   = "Male")

age_levels <- c("15–19", "20–24", "25–29", "30–34","35–39", "40–44", "45–49")

plot_df <- Subset1 %>%
  filter(!is.na(log10_ARG_load),
         !is.na(sex),
         !is.na(AgeGroup),
         !is.na(World_Bank_Income_Group)) %>%
  mutate(
    AgeGroup = factor(AgeGroup, levels = age_levels),
    sex = factor(sex, levels = c("Female", "Male"))
  )

income_levels <- c("Low income", "Lower middle income", "Upper middle income", "High income")
plot_df$World_Bank_Income_Group <- factor(plot_df$World_Bank_Income_Group, levels = income_levels)

ggplot(plot_df, aes(x = AgeGroup, y = log10_ARG_load, fill = sex)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7,
               outlier.shape = NA,
               alpha = 0.8) +
  facet_wrap(~World_Bank_Income_Group, scales = "free_x") +
  scale_fill_npg() +
  labs(
    title = "ARG Load (log10) by AgeGroup, Sex, separated by World Bank Income Group",
    x = "Age Group (years)",
    y = expression(log[10]*"(ARG load)"),
    fill = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", size = 1),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.spacing = unit(1.2, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave("Boxplot_ARG_load_sex_income_filtered_age.png", width = 8, height = 6, dpi = 300)
```
![Boxplot of ARG load, sex, income and filtered age ](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_and_Income_analyses/Boxplot_ARG_load_sex_income_filtered_age.png)



### 5.2 Table for counts for boxplot of ARG load, sex, income, and filtered age 
```r
ggplot(sex_age_income_counts, aes(x = AgeGroup, y = N, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = N),
            position = position_dodge(width = 0.8),
            vjust = -0.5,  # slightly above the bar
            size = 3) +
  scale_fill_npg() +
  facet_wrap(~World_Bank_Income_Group, scales = "free_x") +
  labs(title = "Sample Counts by Age, Sex, and Income Group",
       x = "Age Group",
       y = "Count",
       fill = "Sex") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(face = "bold"))

ggsave("Boxplot_counts_ARG_load_sex_income_filtered_age.png", width = 8, height = 6, dpi = 300)

```
![Boxplot of counts for ARG load, sex, income and filtered age ](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_and_Income_analyses/Boxplot_counts_ARG_load_sex_income_filtered_age.png)

### 5.3 Linear models

```
# Additive
lm_additive <- lm(log10_ARG_load ~ age_years + sex, data = plot_df)
summary(lm_additive)

# Interaction
lm_interaction <- lm(log10_ARG_load ~ age_years * sex, data = plot_df)
summary(lm_interaction)

# Covariete
lm_cov <- lm(log10_ARG_load ~ age_years * sex + World_Bank_Income_Group, data = plot_df)
summary(lm_cov)

```
**Additive Model**

| Term        | Estimate   | Std. Error | t value | p-value      |
|------------|-----------:|-----------:|--------:|-------------|
| (Intercept)| 2.5909485 | 0.0150759 | 171.861 | < 2e-16 *** |
| age_years  | 0.0050173 | 0.0004902 | 10.235  | < 2e-16 *** |
| sexMale    | -0.0702190| 0.0091225 | -7.697  | 1.72e-14 *** |

Model statistics:
- Residual standard error = 0.2947 on 4191 DF  
- Multiple R-squared = 0.03661, Adjusted R-squared = 0.03615  
- F-statistic = 79.63 on 2 and 4191 DF, p-value < 2.2e-16


**Interaction Model:**

| Term            | Estimate   | Std. Error | t value | p-value      |
|-----------------|-----------:|-----------:|--------:|-------------|
| (Intercept)     | 2.5435866 | 0.0199937 | 127.219 | < 2e-16 *** |
| age_years       | 0.0067090 | 0.0006786 | 9.887   | < 2e-16 *** |
| sexMale         | 0.0295811 | 0.0291830 | 1.014   | 0.310812    |
| age_years:sexMale | -0.0035270 | 0.0009798 | -3.600 | 0.000322 *** |

Model statistics: 
- Residual standard error = 0.2942 on 4190 DF  
- Multiple R-squared = 0.03958, Adjusted R-squared = 0.03889  
- F-statistic = 57.56 on 3 and 4190 DF, p-value < 2.2e-16

**Covariete Model**

| Term                                      | Estimate   | Std. Error | t value | p-value     |
|------------------------------------------|-----------:|-----------:|--------:|------------|
| (Intercept)                               | 2.6776108 | 0.0428770 | 62.449  | < 2e-16 *** |
| age_years                                 | 0.0056514 | 0.0006593 | 8.571   | < 2e-16 *** |
| sexMale                                   | 0.0231541 | 0.0282596 | 0.819   | 0.412641    |
| World_Bank_Income_GroupLower middle income | 0.1101279 | 0.0426455 | 2.582   | 0.009845 ** |
| World_Bank_Income_GroupUpper middle income | 0.1847557 | 0.0443868 | 4.162   | 3.21e-05 ***|
| World_Bank_Income_GroupHigh income       | -0.1261621| 0.0380131 | -3.319  | 0.000911 ***|
| age_years:sexMale                         | -0.0034764| 0.0009485 | -3.665  | 0.000250 ***|

Model statistics:
- Residual standard error = 0.2845 on 4187 degrees of freedom  
- Multiple R-squared = 0.1029, Adjusted R-squared = 0.1016  
- F-statistic = 80.01 on 6 and 4187 DF, p-value < 2.2e-16  


### 5.4 LOESS of ARG load, sex, income, and filtered age 

```r
plot_df <- Subset1 %>%
  filter(!is.na(log10_ARG_load),
         !is.na(sex),
         !is.na(age_years),
         !is.na(World_Bank_Income_Group)) %>%
  mutate(
    AgeGroup = factor(AgeGroup, levels = age_levels),
    sex = factor(sex, levels = c("Female", "Male"))
  )

# LOESS curves with individual points and faceted by income group
ggplot(plot_df, aes(x = age_years, y = log10_ARG_load, color = sex, fill = sex)) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.2, size = 1.2, color = "grey30") +  # show points
  geom_smooth(
    method = "loess",
    se = TRUE,
    span = 0.7,
    alpha = 0.2,
    linewidth = 1.2
  ) +
  scale_color_npg() +
  scale_fill_npg() +
  facet_wrap(~World_Bank_Income_Group, scales = "free_x") +
  labs(
    title = "LOESS Curve of ARG Load by Age, Sex, and Income Group",
    x = "Age (years)",
    y = expression(log[10]*"(ARG load)"),
    color = "Sex",
    fill = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", size = 1),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.spacing = unit(1.2, "lines"),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave("LOESS_ARG_load_sex_income_filtered_age.png", width = 8, height = 6, dpi = 300)
```
![LOESS of ARG load, sex, income and filtered age ](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_and_Income_analyses/LOESS_ARG_load_sex_income_filtered_age.png)

### 5.5 GAM model
```
gam_income <- gam(
  log10_ARG_load ~ 
    sex +
    World_Bank_Income_Group +
    s(age_years, by = sex),
  family = gaussian(),
  data = plot_df,
  method = "REML"
)

summary(gam_income)
```

** Generalized Additive Model Results: ARG Load**
Reference groups: 
- Sex: Female  
- World Bank income group: Low income  

Model Information:

| Item | Value |
|-----|------|
| Family | gaussian |
| Link function | identity |
| Formula | `log10_ARG_load ~ sex + World_Bank_Income_Group + s(age_years, by = sex)` |

Parametric Coefficients:

| Term | Estimate | Std. Error | t value | Pr(>t) | Significance |
|-----|---------:|-----------:|--------:|---------:|:------------|
| (Intercept) | 2.835366 | 0.037788 | 75.033 | < 2e-16 | *** |
| sexMale | -0.071526 | 0.008878 | -8.057 | 1.01e-15 | *** |
| Lower middle income | 0.113561 | 0.042316 | 2.684 | 0.00731 | ** |
| Upper middle income | 0.183807 | 0.044000 | 4.177 | 3.01e-05 | *** |
| High income | -0.123971 | 0.037794 | -3.280 | 0.00105 | ** |

**Significance codes:**  
`***` p < 0.001, `**` p < 0.01, `*` p < 0.05, `.` p < 0.1

Smooth Terms:

| Smooth term | edf | Ref.df | F | p-value | Significance |
|------------|----:|-------:|--:|--------:|:------------|
| s(age_years):sexFemale | 7.567 | 8.48 | 14.062 | < 2e-16 | *** |
| s(age_years):sexMale | 8.102 | 8.78 | 9.3_


## 6. Analysis of ARG load,sex and, BMI, and filtered age (Women of reproductive age (15-49 years))

### 6.1 Boxplot of ARG load, sex, BMI and filtered age 

```
colData_subset_clean <- Subset1 %>%
  filter(BMI_range_new != "Normal/Overweight (<30)",
         !is.na(log10_ARG_load),
         !is.na(sex),
         !is.na(AgeGroup),
         !is.na(BMI_range_new)) %>%
  mutate(
    AgeGroup = factor(AgeGroup, levels = c("15–19", "20–24", "25–29", "30–34",
                                           "35–39", "40–44", "45–49")),
    sex = recode(sex,
                 "female" = "Female",
                 "male" = "Male"),
    sex = factor(sex, levels = c("Female", "Male")),
    BMI_range_new = factor(BMI_range_new,
                           levels = c("Underweight (<18.5)",
                                      "Normal (18.5-25)",
                                      "Overweight (25-30)",
                                      "Obese (>30)"))
  )

counts <- colData_subset_clean %>%
  group_by(BMI_range_new, AgeGroup, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  # stagger the y-position of counts by sex
  mutate(
    y_pos = max(colData_subset_clean$log10_ARG_load, na.rm = TRUE) + 
            ifelse(sex == "Female", 0.05, 0.1) # Female lower, Male higher
  )

ggplot(colData_subset_clean, aes(x = AgeGroup, y = log10_ARG_load, fill = sex)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7,
               outlier.shape = NA,
               alpha = 0.8) +
  facet_wrap(~BMI_range_new, scales = "free_x") +
  scale_fill_npg() +
  # Add counts above each box, colored by sex
  geom_text(
    data = counts,
    aes(x = AgeGroup, y = y_pos, label = N, color = sex),
    position = position_dodge(width = 0.8),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold"
  ) +
  scale_color_npg() + # match text color to Sex
  labs(
    title = "ARG Load (log10) by AgeGroup and Sex\nFaceted by BMI Range",
    x = "Age Group (years)",
    y = expression(log[10]*"(ARG load)"),
    fill = "Sex",
    color = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", size = 1),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.spacing = unit(1.2, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )


ggsave("Boxplot_counts_ARG_load_sex_BMI_filtered_age.png", width = 8, height = 6, dpi = 300)

```
![Boxplot of counts for ARG load, sex, income and filtered age ](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_filtered_Age_Analyses/Boxplot_counts_ARG_load_sex_BMI_filtered_age.png)

## 6.2 LOESS curve of ARG load, sex, BMI and filtered age 
```

df <- Subset1 %>%
  filter(BMI_range_new != "Normal/Overweight (<30)",
         !is.na(log10_ARG_load),
         !is.na(sex),
         !is.na(age_years),
         !is.na(BMI_range_new)) %>%
  mutate(
    sex = recode(sex,
                 "female" = "Female",
                 "male" = "Male"),
    sex = factor(sex, levels = c("Female", "Male")),
    BMI_range_new = factor(BMI_range_new,
                           levels = c("Underweight (<18.5)",
                                      "Normal (18.5-25)",
                                      "Overweight (25-30)",
                                      "Obese (>30)"))
  )

ggplot(df, aes(x = age_years, y = log10_ARG_load, color = sex, fill = sex)
) +
  geom_jitter(
    width = 0.3,
    height = 0,
    alpha = 0.2,
    size = 1.2,
    color = "grey30"
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    span = 0.7,
    alpha = 0.2,
    linewidth = 1.2
  ) +
  scale_color_npg() +
  scale_fill_npg() +
  facet_wrap(~ BMI_range_new, scales = "free_x") +
  labs(
    title = "LOESS Curve of ARG Load by Age, Sex, and BMI Category",
    x = "Age (years)",
    y = expression(log[10]*"(ARG load)"),
    color = "Sex",
    fill = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(1.2, "lines"),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave("Loess_ARG_load_sex_BMI_filtered_age.png", width = 8, height = 6, dpi = 300)

```
![Loess of ARG load, sex, BMI and filtered age ](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_filtered_Age_Analyses/Loess_ARG_load_sex_BMI_filtered_age.png)

### 6.3 Linear model of ARG load, sex, BMI and filtered age 

* Three way interaction model

```
lm_3way <- lm(
  log10_ARG_load ~ age_years * sex * BMI_range_new,
  data = df
)

summary(lm_3way)
```

** Linear Model Key Results: ARG Load **

| Term | Estimate | Std. Error | t value | Pr(>t) | Significance | Interpretation |
|------|---------|-----------|--------|---------|:------------|----------------|
| (Intercept) | 2.815 | 0.099 | 28.42 | <2e-16 | *** | Baseline ARG load for Female Underweight at age 0 |
| BMI Normal (18.5–25) | -0.443 | 0.105 | -4.21 | 2.66e-05 | *** | Lower ARG load than Underweight females |
| BMI Overweight (25–30) | -0.366 | 0.144 | -2.54 | 0.011 | * | Lower ARG load than Underweight females |
| BMI Obese (>30) | -0.364 | 0.109 | -3.35 | 0.00081 | *** | Lower ARG load than Underweight females |
| sexMale | -0.239 | 0.171 | -1.40 | 0.162 | | Male vs Female at age 0 (NS) |
| age_years | 0.00190 | 0.00342 | 0.557 | 0.578 | | Age effect for Female Underweight (NS) |
| age:sexMale | 0.0139 | 0.00585 | 2.38 | 0.018 | * | Age slope is steeper for Males vs Females (Underweight) |
| age:BMI Normal | 0.0103 | 0.00363 | 2.85 | 0.004 | ** | Age slope slightly higher for Normal BMI (Females) |
| age:BMI Obese | 0.00852 | 0.00377 | 2.26 | 0.024 | * | Age slope higher for Obese (Females) |
| age:sexMale:BMI Obese | -0.0274 | 0.00624 | -4.39 | 1.20e-05 | *** | Age slope reduced for Male Obese (three-way interaction) |
| age:sexMale:BMI Normal | -0.0154 | 0.00613 | -2.51 | 0.012 | * | Age slope reduced for Male Normal |
| age:sexMale:BMI Overweight | -0.0142 | 0.00686 | -2.07 | 0.039 | * | Age slope reduced for Male Overweight |

* Residual standard error: 0.2915 on 2480 degrees of freedom
* Multiple R-squared:  0.1248,	Adjusted R-squared:  0.1195 
* F-statistic: 23.57 on 15 and 2480 DF,  p-value: < 2.2e-16


## 6.4 Combination Boxplot of ARG load, sex, BMI, income groups, and filtered age 

```
colData_subset_clean <- Subset1 %>%
  filter(BMI_range_new != "Normal/Overweight (<30)",
         !is.na(log10_ARG_load),
         !is.na(sex),
         !is.na(AgeGroup),
         !is.na(BMI_range_new),
         !is.na(World_Bank_Income_Group)) %>%
  mutate(
    AgeGroup = factor(AgeGroup, levels = c("15–19", "20–24", "25–29", "30–34",
                                           "35–39", "40–44", "45–49")),
    sex = recode(sex, "female" = "Female", "male" = "Male"),
    sex = factor(sex, levels = c("Female", "Male")),
    BMI_range_new = factor(BMI_range_new,
                           levels = c("Underweight (<18.5)",
                                      "Normal (18.5-25)",
                                      "Overweight (25-30)",
                                      "Obese (>30)")),
    World_Bank_Income_Group = factor(World_Bank_Income_Group,
                                     levels = c("Low income", "Lower middle income",
                                                "Upper middle income", "High income"))
  )

counts <- colData_subset_clean %>%
  group_by(BMI_range_new, AgeGroup, sex, World_Bank_Income_Group) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(
    y_pos = max(colData_subset_clean$log10_ARG_load, na.rm = TRUE) + 
            ifelse(sex == "Female", 0.05, 0.1)
  )

ggplot(colData_subset_clean, aes(x = AgeGroup, y = log10_ARG_load, fill = sex)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7,
               outlier.shape = NA,
               alpha = 0.8) +
  # Add counts
  geom_text(data = counts,
            aes(x = AgeGroup, y = y_pos, label = N, color = sex),
            position = position_dodge(width = 0.8),
            inherit.aes = FALSE,
            size = 3,
            fontface = "bold") +
  facet_grid(BMI_range_new ~ World_Bank_Income_Group, scales = "free_x") +
  scale_fill_npg() +
  scale_color_npg() +
  labs(
    title = "ARG Load by Age, Sex, BMI, and Income",
    x = "Age Group (years)",
    y = expression(log[10]*"(ARG load)"),
    fill = "Sex",
    color = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", size = 1),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.spacing = unit(1.2, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )


ggsave("Boxplot_ARG_load_sex_BMI_income_filtered_age.png", width = 8, height = 6, dpi = 300)

```
![Boxplot of ARG load, sex, BMI and filtered age ](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_and_Income_analyses/Boxplot_ARG_load_sex_BMI_income_filtered_age.png)

### 6.5 Heatmap of ARG load, sex, BMI, income groups, and filtered age (RE DO!)

```r
# Aggregate data for heatmap
agg_heat <- colData_subset_clean %>%
  group_by(AgeGroup, BMI_range_new, sex, World_Bank_Income_Group) %>%
  summarise(mean_ARG = mean(log10_ARG_load, na.rm = TRUE), .groups = "drop")

# Plot heatmap with continuous NPG-like colors

# Use npg discrete colors as gradient endpoints
npg_colors <- c("#E64B35", "#4DBBD5")  # red and blue from npg palette

ggplot(agg_heat, aes(x = AgeGroup, y = BMI_range_new, fill = mean_ARG)) +
  geom_tile(color = "white") +
  # Add numeric labels on tiles
  geom_text(aes(label = round(mean_ARG, 2)), size = 3, fontface = "bold") +
  # Facet by Sex and Income
  facet_grid(sex ~ World_Bank_Income_Group) +
  # Continuous gradient using NPG-like colors
  scale_fill_gradient(low = npg_colors[2], high = npg_colors[1]) +
  labs(
    title = "ARG Load by Age, Sex, BMI, and Income",
    x = "Age Group (years)",
    y = "BMI Range",
    fill = expression(log[10]*"(ARG load)")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines")
  )

ggsave("Heatmap_ARG_load_sex_income_filtered_age.png", width = 8, height = 6, dpi = 300)
```
![Heatmap of ARG load, sex, income and filtered age ](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_and_Income_analyses/Heatmap_ARG_load_sex_income_filtered_age.png)


## 7. Analysis of ARG load, sex, UTI, and filtered age (Women of reproductive age (15-49 years))

### 7.1 Boxplot of ARG load, sex, UTI, and filtered age 

```
box_df <- Subset1 %>%
  filter(!is.na(log10_ARG_load),
         !is.na(AgeGroup),
         !is.na(UTI_history),
         !is.na(sex)) %>%
  mutate(
    UTI_history = factor(UTI_history, levels = c("No", "Yes")),
    AgeGroup = factor(AgeGroup, levels = c("15–19", "20–24", "25–29", "30–34",
                                           "35–39", "40–44", "45–49")),
    sex = factor(sex, levels = c("Female", "Male"))
  )

counts <- box_df %>%
  group_by(AgeGroup, UTI_history, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(
    y_pos = max(box_df$log10_ARG_load, na.rm = TRUE) + 0.1
  )

custom_fill_colors <- c("No" = "#332288", 
                        "Yes" = "#D55E00")
custom_text_colors <- custom_fill_colors  # same for labels

ggplot(box_df, aes(x = AgeGroup, y = log10_ARG_load, fill = UTI_history)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7,
               outlier.shape = NA,
               alpha = 0.8) +
  geom_text(data = counts,
            aes(x = AgeGroup, y = y_pos, label = N, color = UTI_history, group = sex),
            position = position_dodge2(width = 0.8, preserve ="single"),
            inherit.aes = FALSE,
            size = 3,
            fontface = "bold") +
  facet_wrap(~sex) +
  scale_fill_manual(values = custom_fill_colors) +
  scale_color_manual(values = custom_text_colors) +
  labs(
    title = "ARG Load (log10) by Age, UTI History, and Sex",
    x = "Age Group (years)",
    y = expression(log[10]*"(ARG load)"),
    fill = "UTI History",
    color = "UTI History"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )

ggsave("Boxplot_counts_ARG_load_sex_UTI_filtered_age.png", width = 8, height = 6, dpi = 300)

```
![Boxplot of counts for ARG load, sex, UTI, and filtered age ](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_filtered_Age_Analyses/Boxplot_counts_ARG_load_sex_UTI_filtered_age.png)
