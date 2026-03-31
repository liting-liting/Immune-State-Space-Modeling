############################################################
# TCR Immune Repertoire Dynamics Analysis
# Author: Li Ting
# Description:
#   PCA-based immune state-space analysis including:
#   - State radius
#   - Temporal step distance
#   - Stability index
#   - Effect size (Cohen's d)
############################################################

## ========================
## 1. Load packages
## ========================
library(ggplot2)
library(dplyr)
library(ggpubr)
library(effsize)

## ========================
## 2. Load data
## ========================
df <- read.csv("data/PCA_input.csv")

# Standardize column names
df <- df %>%
  rename(
    PC1 = X1,
    PC2 = X2,
    group = Group3,
    individual = lables
  )

## ========================
## 3. Theme
## ========================
cns_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

## ========================
## 4. PCA trajectory plot
## ========================
p_pca <- ggplot(df, aes(x = PC1, y = PC2, color = group)) +
  geom_path(aes(group = individual), alpha = 0.5) +
  geom_point(size = 2) +
  labs(title = "Immune State-Space Trajectories") +
  cns_theme

ggsave("results/figures/PCA_trajectory.png", p_pca, width = 8, height = 6)

## ========================
## 5. State radius
## ========================
radius_df <- df %>%
  group_by(individual, group) %>%
  summarise(
    mean_PC1 = mean(PC1),
    mean_PC2 = mean(PC2),
    R = sqrt(mean((PC1 - mean_PC1)^2 + (PC2 - mean_PC2)^2)),
    .groups = "drop"
  )

p_radius <- ggplot(radius_df, aes(x = group, y = R, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(y = "Immune State Radius") +
  cns_theme

ggsave("results/figures/state_radius.png", p_radius, width = 6, height = 6)

## ========================
## 6. Step distance
## ========================
step_df <- df %>%
  arrange(individual, time) %>%
  group_by(individual, group) %>%
  mutate(
    step_dist = sqrt((PC1 - lag(PC1))^2 + (PC2 - lag(PC2))^2)
  ) %>%
  summarise(
    D = mean(step_dist, na.rm = TRUE),
    .groups = "drop"
  )

p_step <- ggplot(step_df, aes(x = group, y = D, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(y = "Mean Temporal Step Distance") +
  cns_theme

ggsave("results/figures/step_distance.png", p_step, width = 6, height = 6)

## ========================
## 7. Stability index
## ========================
stability_df <- read.csv("data/stability_score.csv")

p_stability <- ggplot(stability_df, aes(x = group, y = Stability_score, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(y = "Stability Index") +
  cns_theme

ggsave("results/figures/stability.png", p_stability, width = 6, height = 6)

## ========================
## 8. Effect size (Cohen's d)
## ========================
compute_cohens_d <- function(df, metric){
  groups <- unique(df$group)
  res <- list()
  
  for(i in 1:(length(groups)-1)){
    for(j in (i+1):length(groups)){
      g1 <- groups[i]
      g2 <- groups[j]
      
      d_val <- cohen.d(
        df[[metric]][df$group == g1],
        df[[metric]][df$group == g2]
      )
      
      res[[paste(g1, "vs", g2)]] <- d_val
    }
  }
  return(res)
}

radius_d <- compute_cohens_d(radius_df, "R")
step_d <- compute_cohens_d(step_df, "D")
stability_d <- compute_cohens_d(stability_df, "Stability_score")

print(radius_d)
print(step_d)
print(stability_d)