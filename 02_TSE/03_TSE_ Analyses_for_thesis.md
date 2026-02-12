
# 1. Prepare the data
## 1.1 Load packages
```r
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
library(mgcv)
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
# 2. ARG load analyses

## 2.2 Ananlyses of ARGlog10 and sex

### 2.2.1 Boxplot of ARGlog10 and sex
```r

Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$sex[Subset$log10_ARG_load == "" | Subset$log10_ARG_load == "NA"] <- NA

Subset$sex <- recode(Subset$sex,
                     "female" = "Female",
                     "male"   = "Male")

plot_df <- Subset %>% filter(!is.na(log10_ARG_load), !is.na(sex))

n_df <- plot_df %>%
  group_by(sex) %>%
  summarise(
    n = n(),
    y_max = max(ARG_div_shan, na.rm = TRUE)
  )

npg_cols <- pal_npg("nrc")(2)

ggplot(plot_df, aes(x = sex, y = log10_ARG_load, fill = sex)) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.25, color = "grey30") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8) +
  geom_text(data = n_df,
            aes(x = sex,
                y = max(plot_df$log10_ARG_load, na.rm = TRUE) + 0.1,
                label = paste0("N = ", n)),
            inherit.aes = FALSE, size = 4) +
  labs(title = "ARG Load by Sex",
       x = "Sex",
       y = expression(log[10]*"(ARG load)")) +
  scale_fill_npg() +  # <-- NPG palette
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))) +
    scale_x_discrete(labels = c(female = "Female", male = "Male"))

ggsave("Boxplot_log10ARG_by_sex_ready.png", width = 8, height = 6, dpi = 300)

```
![ARG Load by Sex](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_Load_Analyses/Boxplot_log10ARG_by_sex_ready.png)

### 2.2.2 Do we have repeated samples?

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
**Results:***
* n_subjects: 14775
* max_samples_per_subject: 1
* median_samples_per_subject:1 
* subjects_with_repeats :0



### 2.2.3 Wilcoxon rank-sum test
```r


```
**Results:**


### 2.2.4 Count Cohen's d
```r

### Cound Cohens d
```r
group1 <- plot_df$log10_ARG_load[plot_df$sex == "Female"]
group2 <- plot_df$log10_ARG_load[plot_df$sex == "Male"]

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
# -> 0.1136723
**Results**
Cohen's D: 0.1136723

---

## 2.3 Ananlyses of ARGlog10, sex and age

### 2.3.1 New boxplot of ARGlog10, sex and age categories

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
  "Child" = "Children\n(3 -11)",
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
![ARG Load by Sex and age categories](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_Load_Analyses/Boxplot_log10ARG_by_sex_age_ready.png)

### 2.3.2 Loess curve of ARG load by sex and age categories
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
![ARG Load by Sex and age categories](https://github.com/Karhusa/F_AMR_project/blob/main/Results/AGR_Load_Analyses/Boxplot_log10ARG_by_sex_age_category_ready.png)


### 2.3.3 Loess curve of ARG load by sex and numeric age

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


### 2.3.4 Linear model
```

lm_full <- lm(log10_ARG_load ~ age_years * sex, data = plot_df)
summary(lm_full)

```

| Section | Term | Estimate / Value | Std. Error | Statistic | p-value | Signif. |
|--------|------|-----------------|------------|-----------|---------|---------|
| Model information | Formula | log10_ARG_load ~ age_years * sex | – | – | – | – |
| Model information | Residual standard error | 0.3024 (df = 10065) | – | – | – | – |
| Model information | Observations deleted | 311 | – | – | – | – |
| Model information | Multiple R-squared | 0.01433 | – | – | – | – |
| Model information | Adjusted R-squared | 0.01404 | – | – | – | – |
| Model information | F-statistic | 48.79 (3, 10065 DF) | – | – | <2.2e-16 | – |
| Parametric coefficients | (Intercept) | 2.7052081 | 0.0081720 | t = 331.032 | <2e-16 | *** |
| Parametric coefficients | age_years | 0.0011343 | 0.0001834 | t = 6.185 | 6.45e-10 | *** |
| Parametric coefficients | sexMale | -0.0336205 | 0.0114390 | t = -2.939 | 0.0033 | ** |
| Parametric coefficients | age_years:sexMale | 0.0006515 | 0.0002519 | t = 2.587 | 0.0097 | ** |
| Significance codes | *** | p < 0.001 | – | – | – | – |
| Significance codes | ** | p < 0.01 | – | – | – | – |
| Significance codes | * | p < 0.05 | – | – | – | – |
| Significance codes | . | p < 0.1 | – | – | – | – |



## 2.3.5 Optional GAM model

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



```r

ggplot(plot_df, aes(x = age_years, y = log10_ARG_load)) +
  # raw points
  geom_point(alpha = 0.2, color = "grey70") +
  
  # ribbon from predictions
  geom_ribbon(
    data = newdat,
    aes(x = age_years, ymin = lwr, ymax = upr, fill = sex),
    alpha = 0.2,
    color = NA,
    inherit.aes = FALSE
  ) +
  
  # predicted lines
  geom_line(
    data = newdat,
    aes(x = age_years, y = fit, color = sex),
    size = 1.5,
    inherit.aes = FALSE
  ) +
  
  # ggsci palettes
  scale_color_npg() +
  scale_fill_npg() +
  
  # labels
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




# 3. Shannon analyses

## 3.1. Analyses of shannon index and sex

### 3.1.1 Boxplot of shannon index and sex
```r

Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$ARG_div_shan[Subset$ARG_div_shan == "" | Subset$ARG_div_shan == "NA"] <- NA

Subset$sex <- recode(Subset$sex,
                     "female" = "Female",
                     "male"   = "Male")

plot_df <- Subset %>% filter(!is.na(sex), !is.na(ARG_div_shan))

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

ggsave("Boxplot_Shannon_diversity_by_sex_age.png", width = 8, height = 6, dpi = 300)
```

![Shannon diversity by sex by Sex](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Shannon_Analyses/Boxplot_Shannon_diversity_by_sex_age.png)

### 3.1.2 Regression analysis of shannon index and sex
```r
plot_df$sex <- factor(plot_df$sex, levels = c("Female", "Male"))
lm_sex <- lm(ARG_div_shan ~ sex, data = plot_df)
summary(lm_sex)

```
**Results:**

| Outcome | Predictor | Reference group | Effect estimate | p-value | Variance explained (R²) |
|--------|-----------|-----------------|-----------------|---------|--------------------------|
| ARG Shannon diversity | Sex (Male vs Female) | Female | −0.034 | 8.1 × 10⁻⁵ | 0.001 |


| Model | Formula | Observations | Residual SE | DF |
|------|--------|--------------|-------------|----|
| Linear regression | ARG_div_shan ~ sex | 14,775 | 0.523 | 14,773 |

| Term | Estimate | Std. Error | t value | p-value | Significance |
|------|---------:|-----------:|--------:|--------:|:------------:|
| Intercept (Female) | 1.903 | 0.006 | 313.34 | < 2 × 10⁻¹⁶ | *** |
| Sex (Male vs Female) | −0.034 | 0.009 | −3.94 | 8.12 × 10⁻⁵ | *** |

| Metric | Value |
|-------|-------|
| R² | 0.00105 |
| Adjusted R² | 0.00098 |
| F-statistic | 15.54 |
| Model p-value | 8.12 × 10⁻⁵ |

| Min | 1Q | Median | 3Q | Max |
|----:|---:|-------:|---:|----:|
| −1.88 | −0.33 | 0.05 | 0.38 | 1.63 |



### 3.1.3 Cohen's d of shannon index categorized age and, sex

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
**Results**
Cohen's d: 0.06485994

---
## 3.2. Analyses of shannon index, sex and, categorized age

### 3.2.1 Boxplot of of shannon index, sex and, categorized age

```r
Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$ARG_div_shan[Subset$ARG_div_shan == "" | Subset$ARG_div_shan == "NA"] <- NA
Subset$precise_age_category[
  Subset$precise_age_category == "" | Subset$precise_age_category == "NA"
] <- NA

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
  "Child" = "Children\n(3 -11)",
  "Teenage" = "Teenagers\n(12 - 19)",
  "Young adult" = "Young Adult\n(20 - 34)",
  "Middle-Age Adult" = "Middle-Aged Adult\n(35 - 64)",
  "Older Adult" = "Older Adult\n(65 - 80)",
  "Oldest Adult" = "Oldest Adult\n(80<)"
)

plot_df <- Subset %>%
  filter(
    !is.na(ARG_div_shan),
    !is.na(sex),
    !is.na(precise_age_category)
  ) %>%
  mutate(
    precise_age_category = factor(precise_age_category,
                                  levels = age_levels),
    sex = factor(sex, levels = c("Female", "Male"))
  )

counts <- plot_df %>%
  group_by(precise_age_category, sex) %>%
  summarise(
    N = n(),
    y_pos = max(ARG_div_shan, na.rm = TRUE) + 0.15,
    .groups = "drop"
  )

npg_cols <- pal_npg("nrc")(4)[3:4]

ggplot(plot_df, aes(x = precise_age_category, y = ARG_div_shan, fill = sex)) +
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
    aes(x = precise_age_category, y = y_pos,
        label = paste0("N=", N), color = sex),
    position = position_dodge(width = 0.6),
    inherit.aes = FALSE,
    size = 2.0,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_x_discrete(labels = age_labels) +
  scale_fill_manual(values = npg_cols) +
  scale_color_manual(values = npg_cols) +
  labs(
    title = "ARG Shannon diversity by Age category and Sex",
    x = "Age Category",
    y = "ARG Shannon diversity",
    fill = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("Boxplot_Shannon_diversity_by_age_categry_sex.png", width = 8, height = 6, dpi = 300)

```

![Shannon diversity by sex by Sex](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Shannon_Analyses/Boxplot_Shannon_diversity_by_age_categry_sex.png)

### 3.2.2 Regression analysis of shannon index, sex and, categorized age
```r
#Additive model
lm_sex_age <- lm(ARG_div_shan ~ sex + precise_age_category, data = plot_df)
summary(lm_sex_age)

#Interaction model
lm_interact <- lm(ARG_div_shan ~ sex * precise_age_category, data = plot_df)
summary(lm_interact)

```
**Results:**

**Table X.** Linear regression analysis of ARG Shannon diversity as a function of sex and age category. Female sex and the youngest age category were used as reference levels.

| Term                              | Estimate (β) | Std. Error | t value | p value        |
|-----------------------------------|--------------|------------|---------|----------------|
| Intercept                         | 1.840        | 0.019      | 97.10   | < 2.2 × 10⁻¹⁶ |
| Sex (Male vs Female)              | −0.057       | 0.010      | −5.92   | 3.24 × 10⁻⁹   |
| Toddler                           | −0.164       | 0.048      | −3.38   | 7.25 × 10⁻⁴   |
| Child                             | −0.069       | 0.026      | −2.66   | 7.76 × 10⁻³   |
| Teenage                           | −0.088       | 0.023      | −3.77   | 1.62 × 10⁻⁴   |
| Young adult                       | −0.004       | 0.021      | −0.18   | 0.854         |
| Middle-aged adult                 | 0.211        | 0.020      | 10.46   | < 2.2 × 10⁻¹⁶ |
| Older adult                       | 0.283        | 0.022      | 13.05   | < 2.2 × 10⁻¹⁶ |
| Oldest adult                      | 0.394        | 0.038      | 10.31   | < 2.2 × 10⁻¹⁶ |


| Predictor            | β (Estimate) | p value        |
|----------------------|--------------|----------------|
| Sex (Male)           | −0.057       | 3.24 × 10⁻⁹   |
| Age category         | —            | < 2.2 × 10⁻¹⁶ |


**Model statistics**

- Significance codes:  
  `***` p < 0.001, `**` p < 0.01, `*` p < 0.05, `.` p < 0.1

**Residual standard error:** 0.485 on 103371 degrees of freedom  
**Multiple R-squared:** 0.084, Adjusted R-squared: 0.084  
**F-statistic:** 119.6 (df= 15, and 10371), p-value: < 2.2 × 10⁻¹⁶

**Additive model**

| Term                                  | Estimate (β) | Std. Error | t value | p value        |
|---------------------------------------|--------------|------------|---------|----------------|
| Intercept                             | 1.827        | 0.026      | 71.32   | < 2.2 × 10⁻¹⁶ |
| Sex (Male vs Female)                  | −0.030       | 0.036      | −0.82   | 0.411          |
| Toddler                               | −0.184       | 0.075      | −2.44   | 1.46 × 10⁻²   |
| Child                                 | −0.081       | 0.037      | −2.17   | 3.03 × 10⁻²   |
| Teenage                               | −0.031       | 0.033      | −0.96   | 0.339          |
| Young adult                           | 0.116        | 0.029      | 4.03    | 5.67 × 10⁻⁵   |
| Middle-aged adult                      | 0.169        | 0.028      | 6.01    | 1.94 × 10⁻⁹   |
| Older adult                            | 0.222        | 0.031      | 7.19    | 6.84 × 10⁻¹³  |
| Oldest adult                           | 0.381        | 0.056      | 6.76    | 1.50 × 10⁻¹¹  |
| SexMale:Toddler                        | 0.028        | 0.098      | 0.29    | 0.773          |
| SexMale:Child                           | 0.019        | 0.052      | 0.36    | 0.717          |
| SexMale:Teenage                         | −0.111       | 0.046      | −2.40   | 1.65 × 10⁻²   |
| SexMale:Young adult                     | −0.258       | 0.041      | −6.27   | 3.81 × 10⁻¹⁰  |
| SexMale:Middle-Age Adult                | 0.084        | 0.040      | 2.11    | 3.46 × 10⁻²   |
| SexMale:Older Adult                      | 0.104        | 0.043      | 2.42    | 1.54 × 10⁻²   |
| SexMale:Oldest Adult                     | 0.020        | 0.076      | 0.26    | 0.796          |

**Significance codes:** 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1  

**Residual standard error:** 0.480 (df = 10,364)  
**Multiple R²:** 0.104  
**Adjusted R²:** 0.103  
**F-statistic:** 80.51 (df = 15, 10,364), p < 2.2 × 10⁻¹⁶




## LOESS curve 

```r


N_sex <- colData_sex_clean %>%
  filter(!is.na(age_years), !is.na(ARG_div_shan)) %>%
  count(sex)

ggplot(colData_sex_clean,
       aes(x = age_years, y = log10_ARG_load, color = sex)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.4, alpha = 0.15
  ) +
  scale_color_manual(
    values = c("female" = "#D55E00", "male" = "#0072B2")
  ) +
  labs(x = "Age (years)", y = "Log10 ARG load", title = "ARG load vs Age by Sex", color = "Sex"
  ) +
  theme_minimal()
```


LOESS

```

Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$ARG_div_shan[Subset$ARG_div_shan == "" | Subset$ARG_div_shan == "NA"] <- NA
Subset$precise_age_category[
  Subset$precise_age_category == "" | Subset$precise_age_category == "NA"
] <- NA

plot_df_loess <- Subset %>%
  dplyr::filter(
    !is.na(age_years),
    !is.na(ARG_div_shan),
    !is.na(sex)
  )
npg_cols <- pal_npg("nrc")(4)[3:4]

ggplot(plot_df_loess, aes(x = age_years, y = ARG_div_shan)) +
  geom_point(alpha = 0.3, size = 1.5, color = "grey50") +  # grey dots
  geom_smooth(aes(color = sex), method = "loess", span = 0.3, se = TRUE) +
  scale_color_manual(values = npg_cols) +  # LOESS curves in NPG colors
  labs(
    x = "Age (years)",
    y = "ARG Shannon diversity",
    title = "LOESS Curve of ARG Shannon diversity across age by sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

#Without dots:

ggplot(plot_df_loess, aes(x = age_years, y = ARG_div_shan, color = sex, fill = sex)) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, size = 1.2) +
  scale_color_manual(values = npg_cols) +
  scale_fill_manual(values = npg_cols) +
  labs(
    x = "Age (years)",
    y = "ARG Shannon diversity",
    title = "LOESS curve of ARG Shannon diversity across age by sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )
  
ggsave("Loess_Shannon_by_sex_age_numeric_ready.png", width = 8, height = 6, dpi = 300)
```

![LOESS curve of Shannon diversity by sex and numeric age](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Shannon_Analyses/Loess_Shannon_by_sex_age_numeric_ready.png)


# GAM Model

```r

gam_model <- gam(ARG_div_shan ~ sex + s(age_years, k = 10), 
                 data = plot_df_loess, 
                 method = "REML")

# Summary of the model
summary(gam_model)


```

| Term               | Type       | Estimate / edf | Std. Error / Ref.df | t / F  | p-value       | Significance |
| ------------------ | ---------- | -------------- | ------------------- | ------ | ------------- | ------------ |
| (Intercept)        | Parametric | 1.9466         | 0.0069              | 282.43 | < 2.2 × 10⁻¹⁶ | ***          |
| sexMale            | Parametric | -0.0639        | 0.0097              | -6.60  | 4.47 × 10⁻¹¹  | ***          |
| s(age_years)       | Smooth     | 7.99 (edf)     | 8.69                | 115.3  | < 2.2 × 10⁻¹⁶ | ***          |
| **Model fit**      | —          | —              | —                   | —      | —             | —            |
| Adjusted R²        | Fit        | 0.093          | —                   | —      | —             | —            |
| Deviance explained | Fit        | 9.38%          | —                   | —      | —             | —            |
| Scale estimate     | Fit        | 0.234          | —                   | —      | —             | —            |
| −REML              | Fit        | 7000.7         | —                   | —      | —             | —            |
| Sample size (n)    | Fit        | 10,069         | —                   | —      | —             | —            |


























