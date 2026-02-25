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

## 2.4 Analyses from age categories excluding "Oldest Adults"

### 2.4.1  Base filtering (age ≤ 80)
```r
df_no_80 <- Subset %>%
  filter(!is.na(age_years), age_years <= 80) %>%
  mutate(
    sex = recode(sex, "female" = "Female", "male" = "Male"),
    sex = factor(sex, levels = c("Female", "Male"))
  )
```
#  2.4.2 Plot 1: Boxplot by precise_age_category
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
### 2.4.3 Plot 2: Loess curve by age_years
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
### 2.4.4 Statistical results
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
# 2.4.5 GAM: ARG load ~ age_years + sex

```r
library(mgcv)
gam_model <- gam(
  log10_ARG_load ~ sex + s(age_years, by = sex),
  data = plot_df_cont,
  method = "REML"
)
summary(gam_model)
```
# 2.4.6 Table for statistics
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
  fmt_number(
    columns = vars(Mean, Median, SD, Min, Max),
    decimals = 2
  ) %>%
  fmt_number(
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
