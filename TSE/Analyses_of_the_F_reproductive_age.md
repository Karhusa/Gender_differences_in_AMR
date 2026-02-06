---
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
