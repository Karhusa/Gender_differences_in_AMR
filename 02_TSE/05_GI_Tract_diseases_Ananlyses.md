ARG Load Analysis by Sex

---
## 1. Packages
```r
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
librarry(ggsci)
```
## 2. Load TSE

```r
tse_path <- "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/TSE.rds"
TSE <- readRDS(tse_path)
```
## 3. Extract colData

```r
colData_df <- as.data.frame(colData(TSE))

subset2 <- colData_df %>%
  select(
    acc,
    age_years,
    log10_ARG_load,
    ARG_div_shan, 
    sex, 
    BMI_range_new, 
    Antibiotics_used, 
    UTI_history, 
    GI_disease_history
  )
```

4. Inspect the GI_disease history

```r
table(subset2$GI_disease_history)

subset2$GI_disease_binary <- ifelse(
  subset2$GI_disease_history == "None",
  "No",
  "Yes"
)

```r
table(subset2$GI_disease_history)
```
* Cholera 89
* Clostridium difficile infection 292
* Crohn's disease 936
* Helminthiasis 298 
* Indeterminate colitis 29
* Intestinal_disease 59
* Intestinal_disease;Colitis 70
* Intestinal_disease;Colitis;Crohns_disease  1
* Intestinal_disease;Colitis;IBD 5
* Intestinal_disease;Crohns_disease;IBD 1
* Intestinal_disease;IBD 10
*  None  22157 
* Schistosomiasis 18
* ulcerative colitis 640
* None   22157

## 5. Combine variables for nicer table

* Crohn's disease	--> IBD
* Ulcerative colitis	--> IBD
* Indeterminate colitis	-->IBD
* Intestinal_disease;Colitis -->	IBD
* Intestinal_disease;IBD -->	IBD
* Cholera --> Bacterial_GI
* Clostridium difficile infection --> Bacterial_GI
* Schistosomiasis --> Remove
* Helminthiasis --> Remove

```r
subset2 <- subset2 %>% separate_rows(GI_disease_history, sep=";")

subset2 <- subset2 %>%
  mutate(
    Disease_group = case_when(
      str_detect(GI_disease_history, regex("crohn|colitis|ibd", ignore_case = TRUE)) ~ "IBD",
      str_detect(GI_disease_history, regex("clostridium", ignore_case = TRUE)) ~ "Cdiff",
      str_detect(GI_disease_history, regex("cholera", ignore_case = TRUE)) ~ "Cholera",
      GI_disease_history == "None" ~ "None",
      TRUE ~ "Other"
    )
  )

subset2 <- subset2 %>%
  filter(Disease_group != "None") %>%
  mutate(
    Disease_group2 = case_when(
      Disease_group == "IBD" ~ "IBD",
      Disease_group %in% c("Cdiff", "Cholera") ~ "Bacterial_GI"
    )
  )

table(subset2$Disease_group2)
```
5 Boxplot of GI disease history

```r

subset2$sex[subset2$sex == ""] <- NA
subset2$log10_ARG_load[subset2$log10_ARG_load == ""] <- NA
subset2$GI_disease_binary[subset2$GI_disease_binary == ""] <- NA

plot_df <- subset2 %>%filter(!is.na(sex), !is.na(log10_ARG_load), !is.na(GI_disease_binary))

plot_df$sex <- recode(plot_df$sex, "female" = "Female", "male" = "Male")

n_df <- plot_df %>%
  group_by(GI_disease_binary, sex) %>%
  summarise(
    n = n(),
    y_max = max(log10_ARG_load, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y_label = y_max + 0.05 + ifelse(sex == "Female", 0.02, 0))

dodge_width <- 0.6

sex_colors <- c(
  "Female" = "#D95F02",  
  "Male"   = "#1B9E77"  
)

ggplot(plot_df, aes(x = GI_disease_binary, y = log10_ARG_load, fill = sex)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_width),
    size = 1.2,
    alpha = 0.2,        # very light points
    color = "grey60"
  ) +
  geom_boxplot(
    position = position_dodge(width = dodge_width),
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.9,       # solid fill
    color = "grey25"
  ) +
  geom_text(
    data = n_df,
    aes(
      x = GI_disease_binary,
      y = y_label,
      label = paste0("N = ", n),
      color = sex
    ),
    position = position_dodge(width = dodge_width),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = sex_colors) +
  scale_color_manual(values = sex_colors) +
  labs(
    title = "log10(ARG load) by GI Disease History and Sex",
    x = "GI Disease History",
    y = "log10(ARG load)",
    fill = "Sex",
    color = "Sex"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave("Boxplot_ARG_load_sex_GI.png", width = 8, height = 6, dpi = 300)

```

![Boxplot ARG load sex GI age](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_and_GI_analyses/Boxplot_ARG_load_sex_GI_age.png)


## 6. Plot Disease_Group2 (Only IBD and Bacterial infections)
```r

subset2$sex[subset2$sex == ""] <- NA
subset2$log10_ARG_load[subset2$log10_ARG_load == ""] <- NA
subset2$Disease_group2[subset2$Disease_group2 == ""] <- NA

plot_df <- subset2 %>%
  filter(!is.na(sex),
         !is.na(log10_ARG_load),
         !is.na(Disease_group2))

counts_df <- subset2 %>% group_by(Disease_group2, sex) %>%  summarise(N = n(), .groups = "drop")

y_max <- max(subset2$log10_ARG_load, na.rm = TRUE)

counts_df <- plot_df %>%
  group_by(Disease_group2, sex) %>%
  summarise(N = n(), .groups = "drop")

y_max <- max(plot_df$log10_ARG_load, na.rm = TRUE)

ggplot(plot_df,
       aes(x = Disease_group2,
           y = log10_ARG_load,
           fill = sex)) +
  
  geom_jitter(aes(color = sex),
              position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width = 0.8),
              alpha = 0.12,
              size = 0.6) +
  
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = NA,
               alpha = 0.85) +
  
  geom_text(data = counts_df,
            aes(x = Disease_group2,
                y = y_max + 0.2,
                label = paste0("n=", N),
                group = sex),
            position = position_dodge(width = 0.8),
            inherit.aes = FALSE,
            size = 3.5) +
  
  scale_fill_manual(values = c("Female" = "#D95F02",
                               "Male" = "#1B9E77")) +
  scale_color_manual(values = c("Female" = "#D95F02",
                                "Male" = "#1B9E77")) +
  
  expand_limits(y = y_max + 0.4) +
  
  theme_bw() +
  labs(
    x = "Disease Group",
    y = "Log10 ARG Load",
    fill = "Sex",
    color = "Sex",
    title = "ARG Load by GI Disease Group and Sex"
  ) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

```


### Linear regression

```r
model1 <- lm(log10_ARG_load ~ IBD + sex, data = subset2)
summary(model1)
```

* lm(formula = log10_ARG_load ~ IBD + sex, data = subset2) 

**Coefficients:** 

| Term | Estimate | Std. Error | t value | Pr(>t) | Significance |
|-------------|---------:|-----------:|--------:|---------:|:------------|
| (Intercept) | 2.81511 | 0.03566 | 78.941| < 2e-16 | *** |
| IBD | -0.06461 | 0.03614 | -1.788 | 0.07407| |
| sexMale | -0.04594 | 0.01547 | 2.969 | < 2e-16 | ** |

* Residuals: Min 1Q Median 3Q Max -0.99232 -0.18817 0.02285 0.17113 0.99336
  
* Residual standard error: 0.2838 on 1348 degrees of freedom (875 observations deleted due to missingness)
* Multiple R-squared: 0.009332
* Adjusted R-squared: 0.007862
* F-statistic: 6.349 on 2 and 1348 DF, p-value: 0.001801




