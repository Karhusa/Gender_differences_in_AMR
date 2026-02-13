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
table(Subset2$GI_disease_history)

#Modify

```


5 Boxplot of GI disease history

```r
Subset2$sex[Subset1$sex == "" | Subset2$sex == "NA"] <- NA
Subset2$sex[Subset1$log10_ARG_load == "" | Subset2$log10_ARG_load == "NA"] <- NA
Subset2$sex[Subset1$GI_disease_history == "" | Subset2$GI_disease_history == "NA"] <- NA

Subset2$sex <- recode(Subset21$sex,
                     "female" = "Female",
                     "male"   = "Male")

plot_df <- Subset1 %>% filter(!is.na(sex),!is.na(ARG_div_shan))

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



```
