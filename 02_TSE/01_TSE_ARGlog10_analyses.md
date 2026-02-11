ARG Load Analysis by Sex

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

colData_subset <- colData_df %>%
  select(
    acc,
    age_years,
    log10_ARG_load, 
    sex, 
    BMI_range_new, 
    Antibiotics_used, 
    UTI_history, 
    GI_disease_history, 
    Cancers_and_adenomas,
    Vancomycin, Ampicillin, Cefotaxime, TicarcillinClavulanate, Gentamicin,
    Clindamycin, Cefazolin, Oxacillin, Cefoxitin, Cefepime, 
    Meropenem, PenicillinG, Amoxacillin, AmpicillinSulbactam, SulfamethoxazoleTrimethoprim,
    precise_age_category, imprecise_age_category
  )
```
---
## 4.Look through specific columns

### 4.1 Precise_age_category
```
colData_subset$precise_age_category[colData_subset$precise_age_category == "Unknown"] <- NA
table(colData_subset$precise_age_category)
sum(!is.na(colData_subset$precise_age_category))
```
**Results:**
* Infant 1406
* Toddler 185
* Child 837
* Teeanage 1127
* Young Adult 2583
* Middle-Age Adult 5666
* Older Adult 1941
* Oldest adult 212
* All together 13957


### 4.2 age_years

```r
table(colData_subset$age_years)
# Many values (some are babys ages like 0.00XXXXXXXXXXXXX, so rounding up might be the best way to make sense to these numbers)

colData_subset$age_years <- round(colData_subset$age_years)

#   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
# 1098  235   34   48   48   53   53   59  113   80  117  123  100   82  101  160  232  119  159 
#  19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37 
# 170  114  268  276  268  189  203  187   91  289  141  129   91  137   76   82   96  138  161 
#  38   39   40   41   42   43   44   45   46   47   48   49   50   51   52   53   54   55   56 
#  62   37  102   59   58   78   86  101   84  110   90   59  130  151   71   98   92  138  108 
#  57   58   59   60   61   62   63   64   65   66   67   68   69   70   71   72   73   74   75 
# 107  148  160  179  173  279  144  167  194  145  146  156  173  169  130   92  119  125   89 
#  76   77   78   79   80   81   82   83   84   85   86   87   88   90   91   92   98   99  100 
#  88   81   63   38   28   30   19   13   13    8   22    8    3    1    6    1    1    4    5 
# 101  102  104  105  106  107  108  109 
#   3    1    2   10    4    4    1    3 

sum(!is.na(colData_subset$age_years))
```
**Results:**
* Ages 0 to 109, all ages presented
* All together 11189 

### 4.3 sex 

```r
colData_subset$sex[colData_subset$sex == "" | colData_subset$sex == "NA"] <- NA  # convert empty/NA strings to actual NA
table(colData_subset$sex)
sum(!is.na(colData_subset$sex))
```
**Results:**
* Female 7426
* Male 7249
* All together 14775

### 4.4. UTI_history

```r
table(colData_subset$UTI_history)
```
**Results**
* No 24355
* Yes 250


### 4.4  BMI_range_new
```r
colData_subset$BMI_range_new <- as.character(colData_subset$BMI_range_new)  
colData_subset$BMI_range_new[colData_subset$BMI_range_new == "" | colData_subset$BMI_range_new == "NA"] <- NA
table(colData_subset$BMI_range_new)
sum(!is.na(colData_subset$BMI_range_new))
```
**Results:**
* Underweight (<18.5)  467
* Normal (18.5-25)    3294
* Overweight (25-30)  1478
* Obese (>30)         1231
* Normal/Overweight (<30) # We will leave this out
* Alltogether 6526

### 4.5 Antibiotics used
```r
table(colData_subset$Antibiotics_used)

```

**Results**
* No 23611
* Yes 994

---
## 5. Ananlyses of ARG load and sex

### 5.1 Boxplot ARG load and sex
```r
plot_df <- colData_subset %>% filter(!is.na(log10_ARG_load), !is.na(sex))

n_df <- plot_df %>% count(sex)

ggplot(plot_df, aes(x = sex, y = log10_ARG_load, fill = sex)) + geom_jitter(width = 0.15, size = 1.2, alpha = 0.25,color = "grey30") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8) +
  geom_text(data = n_df, aes(x = sex, y = max(plot_df$log10_ARG_load, na.rm = TRUE) + 0.1,label = paste0("N = ", n)),
    inherit.aes = FALSE,size = 4) +
  labs(title = "ARG Load by Sex", x = "Sex",y = expression(log[10]*"(ARG load)")) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

ggsave("ARG_load_by_sex.png", width = 8, height = 6, dpi = 300)

```
![ARG Load by Sex](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_load_by_sex.png)

### 5.2. Descriptive statistics
```r
colData_subset %>% filter(!is.na(sex), !is.na(log10_ARG_load)) %>% group_by(sex) %>%
  summarise(mean_ARG = mean(log10_ARG_load), median_ARG = median(log10_ARG_load), sd_ARG = sd(log10_ARG_load),n = n())
```

| Sex | Mean | Median | SD | N |
|-----|-----:|-------:|---:|--:|
| Female | 2.75 | 2.76 | 0.311 | 7,426 |
| Male   | 2.72 | 2.72 | 0.311 | 7,349 |


### 5.3 Wilcoxon test

```r
plot_df <- colData_subset %>% filter(!is.na(log10_ARG_load), !is.na(sex))

# Sample size per sex
n_df <- plot_df %>% count(sex)

# Wilcoxon test
wilcox_res <- wilcox.test(log10_ARG_load ~ sex, data = plot_df)

p_label <- paste0("Wilcoxon p = ", signif(wilcox_res$p.value, 3))

```
## 5.4 Boxplot with wilcoxon test value
```r

ggplot(plot_df, aes(x = sex, y = log10_ARG_load, fill = sex)) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.25, color = "grey30") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8) +
  geom_text(data = n_df, aes(x = sex, y = max(plot_df$log10_ARG_load) + 0.15, label = paste0("N = ", n)),
    inherit.aes = FALSE, size = 4) +
  
  # p-value 
  annotate("text", x = 1.5, y = max(plot_df$log10_ARG_load) + 0.35, label = p_label, size = 4.2, fontface = "italic") +
  
  labs(title = "ARG Load by Sex", x = "Sex", y = expression(log[10]*"(ARG load)")) +
  
  theme_minimal(base_size = 13) + theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave("wilcoxon_ARG_load_by_sex.png", width = 8, height = 6, dpi = 300)
```
![ARG Load by Sex](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Wilcoxon_ARG_load_by_sex.png)


--
## 6. Analyses of ARG Load by Age and Sex

### 6.1 Boxplot of ARG Load by Category and Sex

```r
colData_subset$sex <- factor(colData_subset$sex, levels = c("female", "male"))  # make it a factor

# Compute counts per age category × sex
counts <- colData_subset %>%
  group_by(precise_age_category, sex) %>%
  summarise(
    N = n(),
    y_pos = max(log10_ARG_load, na.rm = TRUE) + 0.2,  # position above box
    .groups = "drop"
  )

age_levels <- c(
  "Infant", "Toddler", "Child", "Teenage", 
  "Young adult", "Middle-Age Adult", "Older Adult", "Oldest Adult"
)

colData_subset <- colData_subset %>%
  mutate(precise_age_category = factor(precise_age_category, levels = age_levels))

ggplot(colData_subset, aes(x = precise_age_category, y = log10_ARG_load, fill = sex)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("female" = "#FF9999", "male" = "#9999FF")) +
  geom_text(
    data = counts,
    aes(x = precise_age_category, y = y_pos, label = paste0("N=", N)),
    position = position_dodge(width = 0.8),
    size = 3
  ) +
  labs(
    x = "Age Category",
    y = "Log10 ARG Load",
    title = "ARG Load by Age Category and Sex",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

setwd("/scratch/project_2008149/USER_WORKSPACES/karhula/DATA")
ggsave("ARG_load_by_age_sex.png", width = 8, height = 6, dpi = 300)
```
![ARG Load by Age and Sex](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_Load_Analyses/ARG_load_by_age_sex.png)

* N values can be removed. I left N values there so that it would be easier to interpret results.

### 6.2 Scatter plot + separate regression lines by sex

```r
colData_sex_clean <- colData_subset %>% filter(!is.na(sex))
N_sex <- colData_sex_clean %>%
  filter(!is.na(age_years), !is.na(log10_ARG_load)) %>%
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

![Regression analysis ARG Load by Age and Sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/Regression_ARG_load_by_age_sex.png)

### 6.3 Linear regression with results
```
model <- lm(log10_ARG_load ~ age_years + sex, data = colData_sex_clean)
summary(model)

model_sum <- summary(model)

beta_age <- model_sum$coefficients["age_years", "Estimate"]
p_age    <- model_sum$coefficients["age_years", "Pr(>|t|)"]

r2 <- model_sum$r.squared
n  <- model_sum$df[1] + model_sum$df[2] + 

ggplot(colData_sex_clean, aes(x = age_years, y = log10_ARG_load)
  ) +
  geom_point(alpha = 0.2, color = "grey60") +
  geom_smooth(aes(color = sex), method = "lm", se = TRUE,linewidth = 1.5
  ) +
  scale_color_manual(
    values = c("female" = "#D55E00", "male" = "#0072B2")
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
  labs(x = "Age (years)", y = "Log10 ARG load", title = "ARG load vs Age by Sex", color = "Sex"
  ) +
  theme_minimal()

ggsave("Regression_with_table_ARG_load_by_age_sex.png", width = 8, height = 6, dpi = 300)
```

![Regression analysis with table ARG Load by Age and Sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/Regression_with_table_ARG_load_by_age_sex.png)

| Section                 | Term        | Estimate / Value | Std. Error | Statistic   | p-value | Signif. |
|-------------------------|------------|----------------|------------|------------|---------|---------|
| Parametric coefficients | (Intercept)| 2.6921235      | 0.0064149  | t = 419.667 | <2e-16 | ***     |
| Parametric coefficients | age_years  | 0.0014798      | 0.0001257  | t = 11.775  | <2e-16 | ***     |
| Parametric coefficients | sexmale    | -0.0084796     | 0.0060318  | t = -1.406  | 0.16   | –       |



| Metric | Value |
|-------|-------|
| Residual standard error | 0.3024 |
| Degrees of freedom | 10066 |
| Observations removed (missingness) | 4706 |
| Multiple R-squared | 0.01369 |
| Adjusted R-squared | 0.0135 |
| F-statistic | 69.87 (2, 10066 DF) |
| Model p-value | < 2.2e-16 |

### 6.4 Generalized Additive Model (GAM)
* A GAM models the outcome as a sum of smooth functions of predictors rather than simple linear effects.

```r

library(mgcv)

gam_model <- gam(
  log10_ARG_load ~ s(age_years) + sex,
  data = colData_sex_clean,
  method = "REML"
)

summary(gam_model)

ggplot(colData_sex_clean, aes(x = age_years, y = log10_ARG_load)) +
  geom_point(alpha = 0.2, color = "grey70") +
  geom_smooth(
    aes(color = sex),
    method = "gam",
    formula = y ~ s(x),
    se = TRUE,
    linewidth = 1.5
  ) +
  scale_color_manual(values = c("female" = "#D55E00", "male" = "#0072B2")) +
  labs(x = "Age (years)", y = "Log10 ARG load", title = "ARG load vs Age by Sex (GAM)", color = "Sex") +
  theme_minimal()

ggsave("GAM_ARG_load_by_age_sex.png", width = 8, height = 6, dpi = 300)
```
![Regression analysis with table ARG Load by Age and Sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/GAM_ARG_load_by_age_sex.png)

| Section | Term | Estimate / Value | Std. Error | Statistic | p-value | Signif. |
|--------|------|------------------|------------|-----------|---------|---------|
| Model information | Family | Gaussian | – | – | – | – |
| Model information | Link function | Identity | – | – | – | – |
| Model information | Formula | log10_ARG_load ~ s(age_years) + sex | – | – | – | – |
| Model information | Sample size (n) | 10069 | – | – | – | – |
| Parametric coefficients | (Intercept) | 2.751674 | 0.004265 | t = 645.156 | <2e-16 | *** |
| Parametric coefficients | sexmale | -0.013133 | 0.005999 | t = -2.189 | 0.0286 | * |
| Smooth terms | s(age_years) | – | – | F = 41.25 (edf = 7.73, Ref.df = 8.392) | <2e-16 | *** |
| Model fit | Adjusted R-squared | 0.033 | – | – | – | – |
| Model fit | Deviance explained | 3.39% | – | – | – | – |
| Model fit | REML | 2171.1 | – | – | – | – |
| Model fit | Scale estimate | 0.089663 | – | – | – | – |

---

## 7. Analyses of ARG Load by BMI and Sex

### 7.1 Boxplot of ARG Load by Category and Sex


```r
colData_subset_clean <- colData_subset %>%
  filter(BMI_range_new != "Normal/Overweight (<30)")

# Reorder BMI categories
colData_subset_clean$BMI_range_new <- factor(
  colData_subset_clean$BMI_range_new,
  levels = c("Underweight (<18.5)", "Normal (18.5-25)", 
             "Overweight (25-30)", "Obese (>30)")
)

counts <- colData_subset_clean %>%
  group_by(BMI_range_new, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(y_pos = max(colData_subset_clean$log10_ARG_load, na.rm = TRUE) + 0.1)

ggplot(colData_subset_clean, aes(x = BMI_range_new, y = log10_ARG_load, fill = sex)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8), na.rm = TRUE) +
  scale_fill_manual(values = c("female" = "#FF9999", "male" = "#9999FF")) +
  geom_text(
    data = counts,
    aes(x = BMI_range_new, y = y_pos, label = paste0("N=", N), color = sex),
    position = position_dodge(width = 0.8),
    size = 3) +
  labs(x = "BMI Category", y = "Log10 ARG Load", title = "ARG Load by BMI Category and Sex", fill = "Sex") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Boxplot_by_BMI_sex.png", width = 8, height = 6, dpi = 300)

```
![Boxplot by BMI and Sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/ARG_Load_Analyses/Boxplot_by_BMI_sex.png)

## 7.2 Calclulations

```r

colData_subset_clean$BMI_range_new <- relevel(colData_subset_clean$BMI_range_new, ref = "Normal (18.5-25)")

model_add_norm <- lm(formula = log10_ARG_load ~ BMI_range_new + sex, data = colData_subset_clean)

summary(model_add_norm)

model_int_norm <- lm(log10_ARG_load ~ BMI_range_new * sex, data = colData_subset_clean)

summary(model_int_norm)

```
**Additive model:**
| Term                              | Estimate  | Std. Error | t value | p-value     | Significance |
|----------------------------------|----------|------------|---------|------------|-------------|
| (Intercept)                       | 2.760985 | 0.006611   | 417.631 | < 2e-16    | ***         |
| BMI_range_newUnderweight (<18.5)  | 0.073392 | 0.015074   | 4.869   | 1.15e-06   | ***         |
| BMI_range_newOverweight (25-30)   | -0.035234| 0.009786   | -3.601  | 0.00032    | **          |
| BMI_range_newObese (>30)          | -0.051683| 0.010229   | -5.053  | 4.48e-07   | ***         |
| sexmale                            | 0.019259 | 0.007715   | 2.496   | 0.01258    | *           |

* Residual standard error: 0.303 on 6252 df
* Multiple R-squared: 0.0117, Adjusted R-squared: 0.0111
* F-statistic: 18.5 on 4 and 6252 df, p = 3.961e-15


**Interaction model:**
| Term                                   | Estimate  | Std. Error | t value | p-value     | Significance |
|---------------------------------------|----------|------------|---------|------------|-------------|
| (Intercept)                            | 2.751587 | 0.007578   | 363.110 | < 2e-16    | ***         |
| BMI_range_newUnderweight (<18.5)       | 0.085481 | 0.019171   | 4.459   | 8.38e-06   | ***         |
| BMI_range_newOverweight (25-30)        | -0.042835| 0.014898   | -2.875  | 0.004051   | **          |
| BMI_range_newObese (>30)               | -0.004661| 0.014251   | -0.327  | 0.743622   |             |
| sexmale                                | 0.037950 | 0.010687   | 3.551   | 0.000386   | ***         |
| BMI_range_newUnderweight:sexmale       | -0.026054| 0.031022   | -0.840  | 0.401031   |             |
| BMI_range_newOverweight:sexmale        | 0.009792 | 0.019745   | 0.496   | 0.619973   |             |
| BMI_range_newObese:sexmale             | -0.096983| 0.020428   | -4.747  | 2.11e-06   | ***         |

* Residual standard error: 0.3024 on 6249 df
* Multiple R-squared: 0.01584, Adjusted R-squared: 0.01474
* F-statistic: 14.37 on 7 and 6249 df, p < 2.2e-16

### 7.3 Interaction model plot

```

colData_subset_clean$BMI_range_new <- factor(
  colData_subset_clean$BMI_range_new,
  levels = c("Underweight (<18.5)", "Normal (18.5-25)", 
             "Overweight (25-30)", "Obese (>30)")
)

counts <- colData_subset_clean %>%
  group_by(BMI_range_new, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(
    y_pos = max(colData_subset_clean$log10_ARG_load, na.rm = TRUE) + 0.05,
    y_offset = ifelse(sex == "female", 0.02, -0.02)  # separate male/female N labels slightly
  )

# Colors for regression lines
line_colors <- c("female" = "#FF6666", "male" = "#6666FF")  

# Plot Interaction Model
ggplot(colData_subset_clean, aes(x = BMI_range_new, y = log10_ARG_load, color = sex, group = sex)) +
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
    y = "Log10 ARG Load",
    title = "ARG Load by BMI × Sex",
    color = "Sex"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Interaction_model_by_BMI_sex.png", width = 8, height = 6, dpi = 300)

```

![Interaction model by BMI sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/Interaction_model_by_BMI_sex.png)


---

## 8. Analyses of ARG Load by UTI and Sex

### 8.1 Boxplot of UTI and sex
```r

colData_subset_clean <- colData_subset %>%
  filter(!is.na(UTI_history) & !is.na(sex))

colData_subset_clean$UTI_history <- factor(
  colData_subset_clean$UTI_history,
  levels = c("No", "Yes"))

colData_subset_clean$sex <- factor(
  colData_subset_clean$sex,
  levels = c("female", "male"))

counts_uti <- colData_subset_clean %>%
  group_by(UTI_history, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(
    y_pos = max(colData_subset_clean$log10_ARG_load, na.rm = TRUE) + 0.1)

ggplot(
  colData_subset_clean,
  aes(x = UTI_history, y = log10_ARG_load, fill = sex)
) +
  geom_boxplot(
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(
    values = c(
      "female" = "#E69F00",  # muted orange
      "male"   = "#56B4E9"   # muted blue
    )
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
  labs(
    x = "UTI History",
    y = "Log10 ARG Load",
    title = "ARG Load by UTI History and Sex",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave("Boxplot_by_UTI_sex_muted.png", width = 7, height = 5, dpi = 300)

```

![Boxplot ARG Load by UTI_history and Sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/Boxplot_by_UTI_sex_muted.png)

### 8.2 Linear model 

```r
lm_uti <- lm(log10_ARG_load ~ UTI_history + sex, data = colData_subset_clean)

summary(lm_uti)

```
| Section                 | Term           | Estimate / Value | Std. Error | Statistic   | p-value   | Signif. |
|-------------------------|----------------|----------------|------------|------------|-----------|---------|
| Parametric coefficients | (Intercept)    | 2.748995       | 0.003631   | t = 757.110 | <2e-16   | ***     |
| Parametric coefficients | UTI_historyYes | 0.066736       | 0.019838   | t = 3.364   | 0.00077  | ***     |
| Parametric coefficients | sexmale        | -0.034817      | 0.005117   | t = -6.804  | 1.06e-11 | ***     |
| Model fit               | Residual standard error | 0.3108 | – | – | – | – |
| Model fit               | Multiple R-squared      | 0.003983 | – | – | – | – |
| Model fit               | Adjusted R-squared      | 0.003849 | – | – | – | – |
| Model fit               | F-statistic             | 29.54 on 2 and 14772 DF | – | – | 1.574e-13 | – |




### 8.3 Interactive model

```r
lm_uti_int <- lm(
  log10_ARG_load ~ UTI_history * sex,
  data = colData_subset_clean
)
summary(lm_uti_int)

```
| Section                 | Term                  | Estimate / Value | Std. Error | Statistic   | p-value   | Signif. |
|-------------------------|----------------------|----------------|------------|------------|-----------|---------|
| Parametric coefficients | (Intercept)           | 2.7489996      | 0.0036456  | t = 754.062 | <2e-16   | ***     |
| Parametric coefficients | UTI_historyYes        | 0.0665069      | 0.0252336  | t = 2.636   | 0.00841  | **      |
| Parametric coefficients | sexmale               | -0.0348266     | 0.0051587  | t = -6.751  | 1.52e-11 | ***     |
| Parametric coefficients | UTI_historyYes:sexmale| 0.0005994      | 0.0408320  | t = 0.015   | 0.98829  | –       |
| Model fit               | Residual standard error | 0.3109       | –          | –          | –         | –       |
| Model fit               | Multiple R-squared      | 0.003983     | –          | –          | –         | –       |
| Model fit               | Adjusted R-squared      | 0.003781     | –          | –          | –         | –       |
| Model fit               | F-statistic             | 19.69 on 3 and 14771 DF | – | – | 9.814e-13 | – |


---

## 9. Analyses of antibiotic use

### 9.1 Boxplot Antibiotics Use and Sex

```r
counts_abx <- colData_subset_clean %>%
  group_by(Antibiotics_used, sex) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(y_pos = max(colData_subset_clean$log10_ARG_load) + 0.1)  # position above max ARG load

# Plot ARG load by Antibiotics_used and Sex with N values
ggplot(colData_subset_clean, aes(x = Antibiotics_used, y = log10_ARG_load, fill = sex)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("female" = "#D55E00",  # muted red-orange
                               "male"   = "#0072B2")) + # muted blue
  geom_text(
    data = counts_abx,
    aes(x = Antibiotics_used, y = y_pos, label = paste0("N=", N)),
    position = position_dodge(width = 0.8),
    size = 3
  ) +
  labs(
    x = "Antibiotics Used",
    y = "Log10 ARG Load",
    title = "ARG Load by Antibiotics Use and Sex",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave("Boxplot_by_AB_Use_sex.png", width = 7, height = 5, dpi = 300)

```

![Boxplot ARG Load by AB Use and Sex](https://github.com/Karhusa/Gender_differences_in_AMR/blob/main/Results/Boxplot_by_AB_Use_sex.png)

### 9.2 Linear model 

```r
lm_ab <- lm(log10_ARG_load ~ Antibiotics_used + sex, data = colData_subset_clean)

summary(lm_ab)
```

| Section                 | Term                  | Estimate / Value | Std. Error | Statistic   | p-value   | Signif. |
|-------------------------|----------------------|----------------|------------|------------|-----------|---------|
| Parametric coefficients | (Intercept)           | 2.742171       | 0.003666   | t = 747.967 | <2e-16   | ***     |
| Parametric coefficients | Antibiotics_usedYes   | 0.116893       | 0.010358   | t = 11.286  | <2e-16   | ***     |
| Parametric coefficients | sexmale               | -0.034034      | 0.005096   | t = -6.678  | 2.5e-11  | ***     |
| Model fit               | Residual standard error | 0.3096       | –          | –           | –         | –       |
| Model fit               | Multiple R-squared      | 0.01174      | –          | –           | –         | –       |
| Model fit               | Adjusted R-squared      | 0.01161      | –          | –           | –         | –       |
| Model fit               | F-statistic             | 87.75 on 2 and 14772 DF | – | – | < 2.2e-16 | – |


```


### 9.3 Interaction model
```
lm_interaction <- lm(log10_ARG_load ~ Antibiotics_used * sex, data = colData_subset_clean)
summary(lm_interaction)
```
