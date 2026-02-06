


## 1. Packages
```r
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
```
---
## 2. Load TSE
```r
tse_path <- "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/TSE.rds"
TSE <- readRDS(tse_path)
```
## 3. Extract colData
```r
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
## 4. Inspect the values

```r
table(Subset1$AgeGroup)

table(Subset1$sex)

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

# 6. Ananlyses of ARGlog10 sex and income

## 6.1 Boxbot of ARGlog10 sex and income
```r
Subset1$sex[Subset$sex == "" | Subset1$sex == "NA"] <- NA
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

## 6.1 Linear model of ARGlog10 index, income and sex

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

lm_main <- lm(
  log10_ARG_load ~ sex + World_Bank_Income_Group, data = plot_df)

summary(lm_main)

```
## Linear model: log10(ARG load) ~ sex + World Bank income group

Reference group: Female, Low income

### Coefficients

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


## 5.2 Interaction model of ARGlog10 index, income and sex
```r
lm_int <- lm(
  log10_ARG_load ~ sex * World_Bank_Income_Group,
  data = plot_df
)

summary(lm_int)


```

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

## 4. Analysis of sex and filtered age (Women of reproductive age (15-49 years))

## 

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
    y_pos = max(plot_df$log10_ARG_load, na.rm = TRUE) + 0.1
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

# Loess curve
```
loess_df <- Subset1 %>%
  filter(
    !is.na(age_years),
    age_years >= 15,
    age_years <= 49,
    !is.na(log10_ARG_load),
    !is.na(sex)
  ) %>%
  mutate(
    sex = recode(sex,
                 "female" = "Female",
                 "male"   = "Male"),
    sex = factor(sex, levels = c("Female", "Male"))
  )



```
