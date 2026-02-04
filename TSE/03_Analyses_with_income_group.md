## 1. Load packages
```r
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
library(effectsize)
```

## 2. Load TSE

```r
tse_path <- "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/TSE.rds"
TSE <- readRDS(tse_path)
```

## 3. Extract colData

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

## 4. Ananlyses of ARGlog10 index and sex

### 5.1 Boxplot of ARGlog10 index and sex
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
    axis.title.y = element_text(margin = margin(r = 10))) +
    scale_x_discrete(labels = c(female = "Female", male = "Male"))

ggsave("Boxplot_log10ARG_by_sex_ready.png", width = 8, height = 6, dpi = 300)

```
## 4.1 Linear model

```r

# Confirmation for the analysis
any(duplicated(plot_df$acc))
# False
table(table(plot_df$acc))
# 1 
# 14775

model <- lm(log10_ARG_load ~ sex, data = plot_df)
summary(model)
```
| Term | Estimate | Std. Error | t value | Pr(>t) |
|------------|-----------|------------|----------|-------------|
| (Intercept) | 2.750388 | 0.003608 | 762.201 | < 2e-16 *** |
| sexmale | -0.035347 | 0.005117 | -6.908 | 5.1e-12 *** |

**Model summary:**

* Residual standard error: 0.311 on 14773 degrees of freedom
* Multiple R-squared: 0.00322, Adjusted R-squared: 0.003153
* F-statistic: 47.73 on 1 and 14773 DF, p-value: 5.099e-12


### 4.2 Count Cohen's d
```r

# Group means and SDs
group_stats <- plot_df %>%
  group_by(sex) %>%
  summarise(
    mean_ARG = mean(log10_ARG_load, na.rm = TRUE),
    sd_ARG   = sd(log10_ARG_load, na.rm = TRUE),
    n        = n()
  )

group_stats

# Pooled standard deviation
mean_f <- 2.750
mean_m <- 2.715
sd_f   <- 0.311
sd_m   <- 0.312
n_f    <- 7426
n_m    <- 7349

# Pooled SD
sd_pooled <- sqrt( ((n_f - 1)*sd_f^2 + (n_m - 1)*sd_m^2) / (n_f + n_m - 2) )
sd_pooled
# -> 0.3114978

d <- (mean_m - mean_f) / sd_pooled
d
# -> -0.1123603
```



# 5. Ananlyses of ARGlog10 index and sex

```r

# Clean data
Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$World_Bank_Income_Group[Subset$World_Bank_Income_Group == "" | Subset$World_Bank_Income_Group == "NA"] <- NA

# Filter complete cases
plot_df <- Subset %>%
  filter(!is.na(log10_ARG_load),
         !is.na(sex),
         !is.na(World_Bank_Income_Group))

# Set income group order
income_levels <- c("Low income", "Lower middle income", "Upper middle income", "High income")
plot_df$World_Bank_Income_Group <- factor(plot_df$World_Bank_Income_Group, levels = income_levels)

# Sample sizes per group
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


ggsave(".png", width = 8, height = 6, dpi = 300)

```

Boxplot of ARG Shannon Diversity by Sex
