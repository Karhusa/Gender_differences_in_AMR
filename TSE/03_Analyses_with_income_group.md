# 1. Load packages
```r
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
```

# 2. Load TSE

```r
tse_path <- "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/TSE.rds"
TSE <- readRDS(tse_path)
```

# 3. Extract colData

```r
colData_df <- as.data.frame(colData(TSE))

colData_new_subset <- colData_df %>%
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

# 4. Ananlyses of ARGlog10 index and sex

```r

Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA

plot_df <- Subset %>% filter(!is.na(log10_ARG_load), !is.na(sex))

n_df <- plot_df %>% count(sex)

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
    axis.title.y = element_text(margin = margin(r = 10))
  )
```
## 4.1 Linear model

```r
model <- lm(log10_ARG_load ~ sex, data = )
summary(model)
```

# 5. Ananlyses of ARGlog10 index and sex

```r
Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$World_Bank_Income_Group[Subset$World_Bank_Income_Group == ""|Subset$World_Bank_Income_Group == "NA"] 

plot_df <- Subset %>% 
  filter(!is.na(log10_ARG_load), 
         !is.na(sex),
         !is.na(World_Bank_Income_Group))

n_df <- plot_df %>% 
  group_by(World_Bank_Income_Group, sex) %>% 
  summarise(n = n(), .groups = "drop")

ggplot(plot_df, aes(x = sex, y = log10_ARG_load, fill = sex)) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.25, color = "grey30") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8) +
  geom_text(data = n_df,
            aes(x = sex,
                y = max(plot_df$log10_ARG_load, na.rm = TRUE) + 0.1,
                label = paste0("N = ", n)),
            inherit.aes = FALSE, size = 4) +
  labs(title = "ARG Load by Sex and Income Group",
       x = "Sex",
       y = expression(log[10]*"(ARG load)")) +
  scale_fill_npg() +   # NPG colors
  facet_wrap(~World_Bank_Income_Group) +  # separate plots by income group
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


ggsave("Boxplot_Shannon_diversity_by_sex.png", width = 8, height = 6, dpi = 300)

```

Boxplot of ARG Shannon Diversity by Sex
