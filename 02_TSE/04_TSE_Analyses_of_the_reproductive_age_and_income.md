
# 1. Prepare the data

### 1.2 Packages
```r
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
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
  filter(age_years >= 15 & age_years <= 49) %>%   # only ages 15–49
  mutate(
    # Create 5-year age groups
    AgeGroup = cut(age_years,
                   breaks = seq(15, 50, by = 5),   # 15-19, 20-24, ..., 45-49
                   right = FALSE,                  # include left, exclude right
                   labels = paste(seq(15, 45, by = 5),
                                  seq(19, 49, by = 5),
                                  sep = "–"))
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
**Results:***
* n_subjects: 4444
* max_samples_per_subject: 1
* median_samples_per_subject:1 
* subjects_with_repeats :0

### 2.1.3Wilcoxon rank-sum test
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
**Results:***
* 0.4051338

### Regression analysis of shannon index and sex

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
**Results:**

Linear model: log10(ARG load) ~ sex + World Bank income group

Reference group: Female, Low income

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

### Model fit

| Metric | Value |
|--------|-------|
| Residual standard error | 0.2982 |
| Degrees of freedom | 12710 |
| R-squared | 0.05172 |
| Adjusted R-squared | 0.05143 |
| F-statistic | 173.3 |
| Model p-value | < 2.2e-16 |


## Linear model: log10(ARG load) ~ sex × World Bank income group

Reference group: Female, Low income

### Coefficients

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
  summarise(N = n(), .groups = "drop") %>%
  mutate(
    y_pos = max(AgeGroup, na.rm = TRUE) + 0.15,
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


# 4.2 Loess curve of ARG load, sex and filtered age 
```
Subset1$sex[Subset1$sex == "" | Subset1$sex == "NA"] <- NA
Subset1$AgeGroup[Subset1$AgeGroup == "" | Subset1$AgeGroup == "NA"] <- NA

Subset1$sex <- recode(Subset1$sex,
                      "female" = "Female",
                      "male"   = "Male")
age_levels <- c("15–19", "20–24", "25–29", "30–34", "35–39", "40–44", "45–49")

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

age_midpoints <- c(
  "15–19" = 17,
  "20–24" = 22,
  "25–29" = 27,
  "30–34" = 32,
  "35–39" = 37,
  "40–44" = 42,
  "45–49" = 47
)

plot_df <- plot_df %>%
  mutate(AgeNum = age_midpoints[as.character(AgeGroup)])

sex_counts <- plot_df %>%
  group_by(sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(label = paste0(sex, " (N=", N, ")"))

plot_df <- plot_df %>%
  mutate(sex_label = factor(sex, levels = sex_counts$sex, labels = sex_counts$label))

# --- Plot LOESS curves only ---
ggplot(plot_df, aes(x = AgeNum, y = log10_ARG_load, color = sex_label, fill = sex_label)) +
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


```

### 4.4 GAM model of of ARG load, sex and filtered age
```r


```


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

### 5.2 Table for counts for 
```r
sex_age_income_counts <- plot_df %>%
  group_by(World_Bank_Income_Group, AgeGroup, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  arrange(World_Bank_Income_Group, AgeGroup, sex)

sex_age_income_counts <- plot_df %>%
  group_by(World_Bank_Income_Group, AgeGroup, sex) %>%
  summarise(N = n(), .groups = "drop")

ggplot(sex_age_income_counts, aes(x = AgeGroup, y = N, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
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


## 6. Analysis of ARG load,sex and, BMI, and filtered age (Women of reproductive age (15-49 years))

### 5.1 Boxplot of ARG load, sex, BMI and filtered age 

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
![Boxplot of counts for ARG load, sex, income and filtered age ](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_filtered_Age_Analyses/Boxplot_counts_ARG_load_sex_BMI_filtered_age.png.png)



## 5.2 LOESS curve
```
age_midpoints <- c(
  "15–19" = 17,
  "20–24" = 22,
  "25–29" = 27,
  "30–34" = 32,
  "35–39" = 37,
  "40–44" = 42,
  "45–49" = 47
)

plot_df <- colData_subset_clean %>%
  mutate(AgeNum = age_midpoints[as.character(AgeGroup)])

ggplot(plot_df, aes(x = AgeNum, y = log10_ARG_load, color = sex, fill = sex)) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    span = 0.7,
    linewidth = 1.2,
    alpha = 0.2
  ) +
  facet_wrap(~BMI_range_new, scales = "free_x") +
  scale_color_npg() +
  scale_fill_npg() +
  labs(
    title = "LOESS Curves of ARG Load by Age, Sex, and BMI Range",
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
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("Loess_ARG_load_sex_BMI_filtered_age.png", width = 8, height = 6, dpi = 300)

```
![Loess of counts for ARG load, sex, BMI and filtered age ](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Loess_ARG_load_sex_BMI_filtered_age.png)

### 5.3 Linear model

* Three way interaction model
  
```


```
