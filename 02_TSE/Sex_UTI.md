## 1.2 Load packages

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
  )
```

```
df_income_uti_ARG <- Subset %>%
  select(
    sex,
    UTI_history,            
    World_Bank_Income_Group,
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
    )
  ) %>%
  select(sex, UTI_history, Income_group, log10_ARG_load)

df_income_uti_ARG$Income_group <- factor(
  df_income_uti_ARG$Income_group,
  levels = c("HIC", "LMIC")  
)
```


```

df_income_uti_ARG_clean <- df_income_uti_ARG %>%
  filter(!is.na(UTI_history),
         !is.na(log10_ARG_load),
         !is.na(sex),
         !is.na(Income_group))

# Set factors
df_income_uti_ARG_clean$UTI_history <- factor(df_income_uti_ARG_clean$UTI_history, levels = c("No", "Yes"))
df_income_uti_ARG_clean$sex <- factor(df_income_uti_ARG_clean$sex, levels = c("Female", "Male"))

# Compute counts per group
counts <- df_income_uti_ARG_clean %>%
  group_by(Income_group, UTI_history, sex) %>%
  summarise(n = n(), .groups = "drop")

# Plot
ggplot(df_income_uti_ARG_clean,
       aes(x = UTI_history,
           y = log10_ARG_load,
           fill = sex)) +
  # Jitter points behind
  geom_jitter(aes(color = sex),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.3,
              size = 1,
              show.legend = FALSE) +
  # Boxplot on top
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = NA,
               alpha = 0.8) +
  # Add sample counts
  geom_text(data = counts,
            aes(x = UTI_history, y = max(df_income_uti_ARG_clean$log10_ARG_load) + 0.2,
                label = paste0("n=", n),
                group = sex),
            position = position_dodge(width = 0.8),
            size = 3) +
  facet_wrap(~ Income_group) +
  scale_fill_npg() +
  scale_color_npg() +
  labs(
    x = "UTI history",
    y = "log10 ARG load",
    fill = "Sex"
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
```
## Linear model:
```

df_HIC <- df_income_uti_ARG_clean %>% filter(Income_group == "HIC")
df_LMIC <- df_income_uti_ARG_clean %>% filter(Income_group == "LMIC")

# Linear model for HIC
lm_HIC <- lm(log10_ARG_load ~ sex + UTI_history, data = df_HIC)
summary(lm_HIC)

# Linear model for LMIC
lm_LMIC <- lm(log10_ARG_load ~ sex + UTI_history, data = df_LMIC)
summary(lm_LMIC)
