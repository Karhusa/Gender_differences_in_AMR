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

5. Combine variables
* Crohn's disease	-->IBD
*  Ulcerative colitis	--> IBD
* Indeterminate colitis	-->IBD
* Intestinal_disease;Colitis -->	IBD
* Intestinal_disease;IBD -->	IBD
* Cholera 
* Clostridium difficile infection 292

```r
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
  mutate(
    Disease_group2 = case_when(
      Disease_group == "IBD" ~ "IBD",
      Disease_group %in% c("Cdiff", "Cholera") ~ "Bacterial_GI",
      Disease_group == "None" ~ "None"
    )
  )
```
5 Boxplot of GI disease history

```r

subset2$sex[subset2$sex == ""] <- NA
subset2$log10_ARG_load[subset2$log10_ARG_load == ""] <- NA
subset2$GI_disease_binary[subset2$GI_disease_binary == ""] <- NA

plot_df <- subset2 %>%
  filter(!is.na(sex), !is.na(log10_ARG_load), !is.na(GI_disease_binary))


plot_df$sex <- recode(plot_df$sex, "female" = "Female", "male" = "Male")

n_df <- plot_df %>%
  group_by(GI_disease_binary, sex) %>%
  summarise(
    n = n(),
    y_max = max(log10_ARG_load, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # small vertical nudge for labels so they don't overlap
  mutate(y_label = y_max + 0.05 + ifelse(sex == "Female", 0.02, 0))

dodge_width <- 0.6

sex_colors <- c(
  "Female" = "#1B9E77",  # muted teal
  "Male"   = "#D95F02"   # muted orange
)

ggplot(plot_df, aes(x = GI_disease_binary, y = log10_ARG_load, fill = sex)) +
  
  # Light jittered points
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_width),
    size = 1.2,
    alpha = 0.2,        # very light points
    color = "grey60"
  ) +
  
  # Darker boxplots
  geom_boxplot(
    position = position_dodge(width = dodge_width),
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.9,       # solid fill
    color = "grey25"
  ) +
  
  # N labels above boxes, colored by sex
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
  
  # Fill and label colors
  scale_fill_manual(values = sex_colors) +
  scale_color_manual(values = sex_colors) +
  
  # Labels
  labs(
    title = "log10(ARG load) by GI Disease History and Sex",
    x = "GI Disease History",
    y = "log10(ARG load)",
    fill = "Sex",
    color = "Sex"
  ) +
  
  # Theme
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave("Boxplot_ARG_load_sex_GI_age.png", width = 8, height = 6, dpi = 300)

```

![Boxplot ARG load sex GI age](https://github.com/Karhusa/F_AMR_project/blob/main/Results/ARG_and_GI_analyses/Boxplot_ARG_load_sex_GI_age.png)

