
## 1. Prepare the data
### 1.1 Load packages
```r
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
library(mgcv)
library(ggpubr)
```

### 1.2 Load TSE

```r
tse_path <- "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/TSE.rds"
TSE <- readRDS(tse_path)
```

### 1.3 Extract colData

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
### 1.4 Do we have repeated samples?

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


## 2 ARG load analyses

### 2.1 Sample distribution (ARG log10 and sex)

```r
Subset_clean <- Subset %>%
  filter(!is.na(log10_ARG_load), !is.na(sex))

ggplot(Subset_clean, aes(x = log10_ARG_load, fill = sex)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.8, position = "identity") +
  facet_wrap(~ sex) +
  scale_fill_npg() +  # Nature Publishing Group palette from ggsci
  labs(
    title = "Distribution of ARG Load by Sex",
    x = "log10 ARG Load",
    y = "Sample Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#Statistics summary
ARG_summary <- Subset_clean %>%
  group_by(sex) %>%
  summarise(
    n = n(),                                # sample count
    mean_ARG = mean(log10_ARG_load),        # mean
    sd_ARG = sd(log10_ARG_load),            # standard deviation
    median_ARG = median(log10_ARG_load),    # median
    min_ARG = min(log10_ARG_load),          # minimum
    max_ARG = max(log10_ARG_load)           # maximum
  )
ARG_summary
```
![Sample distribution ARG Load](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Thesis/Sample_distribution_ARG.png)


### 2.3 Boxplot of ARGlog10 and sex

```r
Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$sex[Subset$log10_ARG_load == "" | Subset$log10_ARG_load == "NA"] <- NA

Subset$sex <- recode(Subset$sex,"female" = "Female", "male"   = "Male")

plot_df <- Subset %>% filter(!is.na(log10_ARG_load), !is.na(sex))

n_df <- plot_df %>%
  group_by(sex) %>%
  summarise(
    n = n(),
    y_max = max(log10_ARG_load, na.rm = TRUE)
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

### 324.3 Linear regression

```r
lm_ARG_sex <- lm(log10_ARG_load ~ sex, data = Subset_clean)

```
**Results:**

Residuals:
     Min       1Q   Median       3Q      Max 
-1.03367 -0.20111  0.00798  0.20251  1.79600 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept)  2.750388   0.003608 762.201  < 2e-16 ***
sexMale     -0.035347   0.005117  -6.908  5.1e-12 ***

---
Signif. codes:  
0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.311 on 14773 degrees of freedom
Multiple R-squared:  0.00322,	Adjusted R-squared:  0.003153 
F-statistic: 47.73 on 1 and 14773 DF,  p-value: 5.099e-12


---

## 4 Ananlyses of ARGlog10, sex and age

### 2.3.1 Sample distribution(ARG x sex x age)

```r
plot_df <- Subset %>%
  filter(
    !is.na(log10_ARG_load),
    !is.na(sex),
    !is.na(precise_age_category),
    !tolower(sex) %in% "Unknown",
  )

table_df <- plot_df %>%
  count(precise_age_category, sex) %>%
  tidyr::pivot_wider(names_from = sex, values_from = n)

table_df <- table_df %>%
  mutate(Total = Female + Male)

table_df

table_image <- table_df %>%
  gt() %>%
  cols_label(
    precise_age_category = "Age Category",
    Female = "Female (n)",
    Male = "Male (n)",
    Total = "Total (n)"
  ) %>%
  tab_header(
    title = "Sample Distribution by Age Category and Sex"
  )

table_image
```

KUVA!!


### 2.3.2 Boxplot (ARG, sex and age categories)

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

plot_df <- plot_df %>% filter(precise_age_category != "NA")

plot_df <- plot_df %>%
  mutate(
    precise_age_category = factor(precise_age_category,
                                  levels = age_levels),
    sex = factor(sex, levels = c("Female", "Male"))
  )

ggplot(plot_df, aes(x = precise_age_category, y = log10_ARG_load, fill = sex)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    size = 1.2, alpha = 0.25, color = "grey30"
  ) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA, alpha = 0.8,
    position = position_dodge(width = 0.6)
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
  plot.title = element_text(face = "plain"),  # not bold
  axis.text.x = element_text(angle = 45, hjust = 1)
)
ggsave("Boxplot_log10ARG_by_sex_age_ready.png", width = 8, height = 6, dpi = 300)

```
![ARG Load by Sex and age categories](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Thesis/Boxplot_ARG_by_sex_age.png)

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

ggplot(plot_df, aes(x = age_years, y = log10_ARG_load, color = sex, fill = sex)) +
  geom_point(alpha = 0.08, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE, span = 0.7, alpha = 0.2, size = 1.2) +
  scale_color_npg() +
  scale_fill_npg() +
  labs(
    title = "ARG Load Across Age by Sex",
    x = "Age (years)",
    y = expression(log[10]*"(ARG load)"),
    color = "Sex",
    fill = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "plain")
  )

ggsave("Loess_log10ARG_by_sex_age_numeric_ready.png", width = 8, height = 6, dpi = 300)

```
![LOESS ARG Load by Sex and age categories](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_Load_Analyses/Loess_log10ARG_by_sex_age_numeric_ready(2).png)


### 2.3.4 Linear model
```
n_samples <- Subset %>%
  filter(
    !is.na(log10_ARG_load),
    !is.na(sex),
    !is.na(age_years)
  ) %>%
  nrow()

n_samples

lm_full <- lm(log10_ARG_load ~ sex + age_years, data = plot_df)
summary(lm_full)

```
** Results**

lm(formula = log10_ARG_load ~ age_years + sex, data = plot_df)

Residuals:
     Min       1Q   Median       3Q      Max 
-1.00348 -0.19450  0.01204  0.19203  1.73418 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)  2.6921136  0.0064172 419.515   <2e-16 ***
age_years    0.0014797  0.0001257  11.769   <2e-16 ***
sexMale     -0.0084758  0.0060318  -1.405     0.16    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.3025 on 10066 degrees of freedom
  (311 observations deleted due to missingness)
Multiple R-squared:  0.01368,	Adjusted R-squared:  0.01348 
F-statistic:  69.8 on 2 and 10066 DF,  p-value: < 2.2e-16




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


---

## 4 Ananlyses of ARG Load, BMI and sex
* Remove ages < 18

### 2.4.1. Sample distribution (ARG x BMI x sex)

```
Subset_BMI <- Subset %>%
  filter(
    age_years >= 18,
    BMI_range_new != "Normal/Overweight (<30)"
  ) %>%
  filter(
    !is.na(log10_ARG_load),
    !is.na(BMI_range_new),
    !is.na(sex)
  )

table_df_B <- Subset_BMI %>%
  count(BMI_range_new, sex) %>%
  tidyr::pivot_wider(names_from = sex, values_from = n)

table_df_B <- table_df_B %>%
  mutate(Total = Female + Male)

table_df_B

table_image <- table_df_B %>%
  gt() %>%
  cols_label(
    BMI_range_new = "BMI Range",
    Female = "Female (n)",
    Male = "Male (n)",
    Total = "Total (n)"
  ) %>%
  tab_header(
    title = "Sample Distribution by Age Category and Sex"
  )

table_image

```

### 2.4.2 Boxplot of ARG Load, BMI and sex
```r

Subset_BMI <- Subset %>%
  filter(
    age_years >= 18,
    BMI_range_new != "Normal/Overweight (<30)"
  ) %>%
  filter(
    !is.na(log10_ARG_load),
    !is.na(BMI_range_new),
    !is.na(sex)
  )

levels_bmi <- c("Underweight (<18.5)", "Normal (18.5-25)", 
                "Overweight (25-30)", "Obese (>30)")

Subset_BMI$BMI_range_new <- factor(
  Subset_BMI$BMI_range_new,
  levels = levels_bmi
)
sex_comparisons <- list(
  c("Female", "Male")  
)

ggplot(Subset_BMI, aes(x = sex, y = log10_ARG_load, fill = sex)) +
  geom_jitter(aes(color = sex),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.2, size = 0.5) +
  geom_boxplot(position = position_dodge(width = 0.75),
               outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ BMI_range_new, nrow = 1) +
  stat_compare_means(
    comparisons = sex_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    position = position_dodge(width = 0.75)
  ) +
  scale_fill_npg(name = "Sex") +
  scale_color_npg(guide = "none") +
  labs(
    x = "BMI Category",
    y = "log10 ARG Load",
    title = "Sex Differences in ARG Load Across BMI Categories"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_blank(),   # remove sex labels
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )
```
### 2.4.3 Linear model (AGR Load x BMI and sex)
```r
Subset_BMI$BMI_range_new <- relevel(
  Subset_BMI$BMI_range_new, 
  ref = "Normal (18.5-25)"
)

lm_bmi <- lm(log10_ARG_load ~ BMI_range_new + sex, data = Subset_BMI)
summary(lm_bmi)

```
Call:
lm(formula = log10_ARG_load ~ BMI_range_new + sex, data = Subset_BMI)

Residuals:
     Min       1Q   Median       3Q      Max 
-1.02651 -0.20462  0.01754  0.20880  1.75810 

Coefficients:
                                  Estimate Std. Error t value
(Intercept)                       2.741763   0.007035 389.706
BMI_range_newUnderweight (<18.5)  0.113370   0.020590   5.506
BMI_range_newOverweight (25-30)  -0.019072   0.010161  -1.877
BMI_range_newObese (>30)         -0.009937   0.011676  -0.851
sexMale                           0.021112   0.008397   2.514
                                 Pr(>|t|)    
(Intercept)                       < 2e-16 ***
BMI_range_newUnderweight (<18.5) 3.84e-08 ***
BMI_range_newOverweight (25-30)    0.0606 .  
BMI_range_newObese (>30)           0.3948    
sexMale                            0.0120 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.3011 on 5245 degrees of freedom
Multiple R-squared:  0.007836,	Adjusted R-squared:  0.007079 
F-statistic: 10.36 on 4 and 5245 DF,  p-value: 2.366e-08


---


# 5. Shannon analyses

## 5.1. ARG diversity and sex

### 2.2 Sample distribution (ARG div x sex)

```r
Subset_shan <- Subset %>%
  filter(!is.na(ARG_div_shan), !is.na(sex))

npg_cols <- pal_npg("nrc")(4)[3:4]

ggplot(Subset_shan, aes(x = ARG_div_shan, fill = sex)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.8, position = "identity") +
  facet_wrap(~ sex, scales = "fixed") +  # separate panels for each sex
    scale_fill_manual(values = npg_cols) +
  labs(
    title = "Distribution of Shannon Diversity Index by Sex",
    x = "Shannon Diversity Index",
    y = "Sample Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

Shan_summary <- Subset_clean %>%
  group_by(sex) %>%
  summarise(
    n = n(),                                # sample count
    mean_shan = mean(ARG_div_shan),        # mean
    sd_shan = sd(ARG_div_shan),            # standard deviation
    median_shan = median(ARG_div_shan),    # median
    min_shan = min(ARG_div_shan),          # minimum
    max_shan = max(ARG_div_shan)           # maximum
  )
Shan_summary
```
![Sample distribution ARG diversity](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Thesis/Sample_distribution_ARG_diversity.png)







### 5.1.1 Boxplot of shannon index and sex
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
plot_df <- Subset %>% filter(!is.na(sex), !is.na(ARG_div_shan))

plot_df$sex <- factor(plot_df$sex, levels = c("Female", "Male"))

lm_sex <- lm(ARG_div_shan ~ sex, data = plot_df)

summary(lm_sex)

```
---
## 3.2. Analyses of shannon index, sex and categorized age

### 3.2.1 Sample distribution (ARG diveristy x sex x categorized age )

```r
plot_df_S <- Subset %>%
  filter(
    !is.na(ARG_div_shan),
    !is.na(sex),
    !is.na(precise_age_category),
    tolower(sex) != "unknown",
    tolower(precise_age_category) != "unknown"
  )

table_df_S <- plot_df_S %>%
  count(precise_age_category, sex) %>%
  tidyr::pivot_wider(names_from = sex, values_from = n)

table_df_S <- table_df_S %>%
  mutate(Total = Female + Male)

table_df_S

table_image <- table_df_S %>%
  gt() %>%
  cols_label(
    precise_age_category = "Age Category",
    Female = "Female (n)",
    Male = "Male (n)",
    Total = "Total (n)"
  ) %>%
  tab_header(
    title = "Sample Distribution by Age Category and Sex"
  )

table_image
```
### 3.2.2 Boxplot of of shannon index, sex and, categorized age

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
    !is.na(ARG_div_shan),
    !is.na(sex),
    !is.na(precise_age_category)
  ) %>%
  mutate(
    precise_age_category = factor(precise_age_category,
                                  levels = age_levels),
    sex = factor(sex, levels = c("Female", "Male"))
  )

plot_df <- plot_df %>% filter(precise_age_category != "NA")

plot_df <- plot_df %>%
  mutate(
    precise_age_category = factor(precise_age_category,
                                  levels = age_levels),
    sex = factor(sex, levels = c("Female", "Male"))
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

### 3.2.2 Regression analysis of shannon index, sex and age
```r
#Additive model
lm_sex_age <- lm(ARG_div_shan ~ sex + age_years, data = plot_df)
summary(lm_sex_age)

```
**Results:**

lm(formula = ARG_div_shan ~ sex + age_years, data = plot_df)

Residuals:
     Min       1Q   Median       3Q      Max 
-1.75710 -0.29371  0.03775  0.33452  1.56527 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)  1.7122904  0.0103191  165.93  < 2e-16 ***
sexMale     -0.0612997  0.0096994   -6.32 2.73e-10 ***
age_years    0.0060274  0.0002022   29.81  < 2e-16 ***
---
Signif. codes:  
0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4864 on 10066 degrees of freedom
  (311 observations deleted due to missingness)
Multiple R-squared:  0.08358,	Adjusted R-squared:  0.0834 
F-statistic:   459 on 2 and 10066 DF,  p-value: < 2.2e-16



```



















