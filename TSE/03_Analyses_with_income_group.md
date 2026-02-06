## 1. Load packages
```r
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
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
**Results**
Cohen's D: -0.1123603


# 5. Ananlyses of ARGlog10 index and sex

## 5.1 Boxbot of ARGlog10 index, income and sex
```r
Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$World_Bank_Income_Group[Subset$World_Bank_Income_Group == "" | Subset$World_Bank_Income_Group == "NA"] <- NA

Subset$sex <- recode(Subset$sex,
                     "female" = "Female",
                     "male"   = "Male")

plot_df <- Subset %>%
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

## 5.1 Linear model of ARGlog10 index, income and sex

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
Call:
lm(formula = log10_ARG_load ~ sex + World_Bank_Income_Group, 
    data = plot_df)

Residuals:
     Min       1Q   Median       3Q      Max 
-1.14215 -0.19192  0.00643  0.18917  1.80070 

Coefficients:
                                            Estimate Std. Error t value Pr(>|t|)    
(Intercept)                                 2.777709   0.021045 131.988  < 2e-16 ***
sexMale                                    -0.025041   0.005346  -4.684 2.84e-06 ***
World_Bank_Income_GroupLower middle income  0.152621   0.026204   5.824 5.87e-09 ***
World_Bank_Income_GroupUpper middle income  0.126062   0.021909   5.754 8.92e-09 ***
World_Bank_Income_GroupHigh income         -0.042325   0.021036  -2.012   0.0442 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.2982 on 12710 degrees of freedom
Multiple R-squared:  0.05172,	Adjusted R-squared:  0.05143 
F-statistic: 173.3 on 4 and 12710 DF,  p-value: < 2.2e-16

## 5.2 Interaction model of ARGlog10 index, income and sex
```r
lm_int <- lm(
  log10_ARG_load ~ sex * World_Bank_Income_Group,
  data = plot_df
)

summary(lm_int)


```

lm(formula = log10_ARG_load ~ sex * World_Bank_Income_Group, 
    data = plot_df)

Residuals:
     Min       1Q   Median       3Q      Max 
-1.11536 -0.19214  0.00655  0.19053  1.80383 

Coefficients:
                                                    Estimate Std. Error t value Pr(>|t|)
(Intercept)                                         2.735527   0.031600  86.567  < 2e-16
sexMale                                             0.049505   0.042008   1.178   0.2386
World_Bank_Income_GroupLower middle income          0.168011   0.038325   4.384 1.18e-05
World_Bank_Income_GroupUpper middle income          0.162007   0.032649   4.962 7.06e-07
World_Bank_Income_GroupHigh income                  0.003468   0.031897   0.109   0.9134
sexMale:World_Bank_Income_GroupLower middle income -0.016688   0.052727  -0.316   0.7516
sexMale:World_Bank_Income_GroupUpper middle income -0.055994   0.044330  -1.263   0.2066
sexMale:World_Bank_Income_GroupHigh income         -0.081284   0.042425  -1.916   0.0554
                                                      
(Intercept)                                        ***
sexMale                                               
World_Bank_Income_GroupLower middle income         ***
World_Bank_Income_GroupUpper middle income         ***
World_Bank_Income_GroupHigh income                    
sexMale:World_Bank_Income_GroupLower middle income    
sexMale:World_Bank_Income_GroupUpper middle income    
sexMale:World_Bank_Income_GroupHigh income         .  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.2981 on 12707 degrees of freedom
Multiple R-squared:  0.05243,	Adjusted R-squared:  0.05191 
F-statistic: 100.4 on 7 and 12707 DF,  p-value: < 2.2e-16
