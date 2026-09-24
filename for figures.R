library(tidyverse)
library(ggpubr)    
library(corrplot)  
library(ggplot2)

data_path <- "C:/Users/user/OneDrive/Desktop/thesis complete data.csv"
df <- read.csv(data_path, stringsAsFactors = FALSE)

df$Land.Type <- trimws(as.character(df$Land.Type))

if (any(is.na(df$Land.Type)) || any(df$Land.Type == "") || length(unique(df$Land.Type)) < 2) {
  # Fallback: assign categories by plot order (10 Forest, 10 Agriculture, 10 Grassland)
  df$Land.Type <- factor(rep(c("Protected Forest", "Peri-Urban Agriculture", "Degraded Grassland"), each = 10),
                         levels = c("Protected Forest", "Peri-Urban Agriculture", "Degraded Grassland"))
} else {
  
  df$Land.Type <- factor(df$Land.Type)
}

df$k <- df$k..day.1.

cat("--- Land Use Category Breakdown ---\n")
print(table(df$Land.Type))

#  FIGURE 2: Boxplots of k and S across Land Types
# Plot A: Decomposition Rate Constant (k)
p_k <- ggplot(df, aes(x = Land.Type, y = k, fill = Land.Type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("#2e7d32", "#f9a825", "#c62828")) +
  theme_bw(base_size = 12) +
  labs(x = "", y = expression(paste("Initial Decomposition Rate ", k, " (day"^-1, ")")),
       title = "A) Decomposition Kinetics (k)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 1))

# Plot B: Soil Stabilization Factor (S)
p_S <- ggplot(df, aes(x = Land.Type, y = S, fill = Land.Type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("#2e7d32", "#f9a825", "#c62828")) +
  theme_bw(base_size = 12) +
  labs(x = "", y = "Soil Stabilization Factor (S)",
       title = "B) Carbon Stabilization Factor (S)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 1))

# Combine into Figure 1 & Save
fig1 <- ggarrange(p_k, p_S, ncol = 2, nrow = 1)
ggsave("Figure1_TBI_Metrics.png", fig1, width = 9, height = 4.5, dpi = 300)
print(fig1)

# One-Way ANOVA & Post-Hoc Tests
cat("\n--- ANOVA: Decomposition Rate (k) ---\n")
fit_k <- aov(k ~ Land.Type, data = df)
print(summary(fit_k))
print(TukeyHSD(fit_k))

cat("\n--- ANOVA: Carbon Stabilization (S) ---\n")
fit_S <- aov(S ~ Land.Type, data = df)
print(summary(fit_S))
print(TukeyHSD(fit_S))

# Multiple Linear Regressions
cat("\n--- OLS Regression: k ---\n")
model_k <- lm(k ~ OM.... + N.... + P2O5.kg.ha. + K..kg.ha. + pH, data = df)
print(summary(model_k))

cat("\n--- OLS Regression: S ---\n")
model_S <- lm(S ~ OM.... + N.... + P2O5.kg.ha. + K..kg.ha. + pH, data = df)
print(summary(model_S))

#  FIGURE 3: Pearson Correlation Heatmap

numeric_df <- df %>% select(pH, Clay = Clay...., OM = OM...., N = N...., 
                            P2O5 = P2O5.kg.ha., K = K..kg.ha., k, S)
cor_matrix <- cor(numeric_df, use = "complete.obs")

png("Figure3_Correlation_Heatmap.png", width = 2000, height = 2000, res = 300)
corrplot(cor_matrix, method = "color", type = "upper", 
         addCoef.col = "black", number.cex = 0.7,
         tl.col = "black", tl.srt = 45, title = "Pearson Correlation Matrix",
         mar = c(0,0,1,0))
dev.off()

# For figure 4
library(ggplot2)

# 1. Define nodes with multiline titles for longer names
nodes <- data.frame(
  id = 1:6,
  x = c(1, 2, 3, 1, 2, 3),
  y = c(2, 2, 2, 1, 1, 1),
  title = c("Protected Forest", "Peri-Urban\nAgriculture", "Degraded Grassland",
            "High Activity", "Balanced Pool", "Low Activity"),
  subtitle = c("(High OM & pH)", "(Nutrient Enriched)", "(Depleted & Acidic)",
               "Fast Kinetic (k)\nLow Stabilization (S)", "Intermediate\nMetrics", "Slow Kinetic (k)\nHigh Stabilization (S)"),
  category = c("Forest", "Agriculture", "Grassland", "Forest", "Agriculture", "Grassland")
)

# 2. Expanded Box Dimensions
nodes$xmin <- nodes$x - 0.44
nodes$xmax <- nodes$x + 0.44
nodes$ymin <- nodes$y - 0.25
nodes$ymax <- nodes$y + 0.25

# Colors
category_colors <- c("Forest" = "#2e6f40", "Agriculture" = "#d97706", "Grassland" = "#dc2626")
category_fill   <- c("Forest" = "#f0fdf4", "Agriculture" = "#fffbeb", "Grassland" = "#fef2f2")

# Flow Arrows
arrows <- data.frame(
  x = c(1, 2, 3), xend = c(1, 2, 3),
  y = c(1.73, 1.73, 1.73), yend = c(1.27, 1.27, 1.27),
  category = c("Forest", "Agriculture", "Grassland")
)

p <- ggplot() +
  # Box Cards
  geom_rect(data = nodes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category, color = category),
            linewidth = 0.8) +
  
  # Arrows
  geom_segment(data = arrows, aes(x = x, xend = xend, y = y, yend = yend, color = category),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"), linewidth = 1.1) +
  
  # Primary Header
  geom_text(data = nodes, aes(x = x, y = y + 0.09, label = title, color = category),
            fontface = "bold", size = 3.6, lineheight = 0.95) +
  
  # Subtitles
  geom_text(data = nodes, aes(x = x, y = y - 0.11, label = subtitle),
            color = "#374151", size = 3.0, lineheight = 0.95) +
  
  scale_color_manual(values = category_colors) +
  scale_fill_manual(values = category_fill) +
  scale_x_continuous(limits = c(0.4, 3.6)) +
  scale_y_continuous(limits = c(0.5, 2.5)) +
  
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "#ffffff", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

# 4. Save the plot
ggsave("ecosystem_diagram.png", plot = p, width = 8, height = 6, dpi = 300)

