# Load necessary library for plotting
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to perform Kruskal-Wallis test and plot the results
compare_categorical_continuous <- function(df, categorical_col, continuous_col) {

  # Perform Kruskal-Wallis test
  kruskal_result <- kruskal.test(df[[continuous_col]] ~ df[[categorical_col]])
  p_value <- format(kruskal_result$p.value, scientific = TRUE)

  # Calculate median continuous values for each categorical group
  median_values <- df %>% 
    group_by(!!sym(categorical_col)) %>% 
    summarise(median_value = median(!!sym(continuous_col))) %>% 
    arrange(median_value)

  # Order categorical variable levels based on median continuous values
  df[[categorical_col]] <- factor(df[[categorical_col]], 
                                   levels = median_values[[categorical_col]])

  # Create boxplot with p-value as subtitle
  plot <- ggplot(df, aes(x = factor(df[[categorical_col]]), y = df[[continuous_col]])) +
    geom_boxplot(notch = TRUE) +
    labs(x = categorical_col, y = continuous_col, 
         subtitle = paste("Kruskal-Wallis p-value =", p_value)) +
    coord_flip() +  # Flip coordinates for horizontal boxplot
    theme_bw()

  # Save the plot
  ggsave(paste0("boxplot_", categorical_col, "_vs_", continuous_col, ".png"), 
         plot = plot,
         width = 6, height = 8)
}

# Create a sorted combination of HTO_maxID and HTO_secondID
sc_total$max.second.ID <- apply(sc_total@meta.data[, c("HTO_maxID", "HTO_secondID")], 1, function(x) {
  paste(sort(x), collapse = "_")
})
sc_total$max.second.ID.library = paste0(sc_total$max.second.ID, "_", sc_total$library_id)

compare_categorical_continuous(sc_total@meta.data, "HTO_classification", "HTO_margin")
compare_categorical_continuous(sc_total@meta.data, "HTO_maxID", "HTO_margin")
compare_categorical_continuous(sc_total@meta.data, "HTO_classification.global", "HTO_margin")
compare_categorical_continuous(sc_total@meta.data, "library_id", "HTO_margin")
compare_categorical_continuous(sc_total@meta.data, "patient_id", "HTO_margin")
compare_categorical_continuous(sc_total@meta.data, "max.second.ID", "HTO_margin")
compare_categorical_continuous(sc_total@meta.data, "max.second.ID.library", "HTO_margin")

# Calculate average margin for each max.second.ID group
avg_margin <- sc_total@meta.data %>% 
  group_by(max.second.ID) %>% 
  summarise(avg_margin = mean(HTO_margin))

# Max scaling to 0-100 and rounding
range_min <- 0
range_max <- 100
avg_margin$scaled_margin <- round((avg_margin$avg_margin / max(avg_margin$avg_margin)) * (range_max - range_min) + range_min)
# Separate max.second.ID by "_" to get individual hashtags
avg_margin <- avg_margin %>% 
  separate(max.second.ID, c("Hashtag1", "Hashtag2"), sep = "_")

# Create upper triangle matrix
unique_hashtags <- unique(c(avg_margin$Hashtag1, avg_margin$Hashtag2))
num_hashtags <- length(unique_hashtags)
upper_triangle_matrix <- matrix(NA, nrow = num_hashtags, ncol = num_hashtags)
colnames(upper_triangle_matrix) <- unique_hashtags
rownames(upper_triangle_matrix) <- unique_hashtags

# Populate the upper triangle matrix
for (i in 1:(nrow(avg_margin))) {
  group1 <- avg_margin$Hashtag1[i]
  group2 <- avg_margin$Hashtag2[i]
  value <- avg_margin$scaled_margin[i]
  upper_triangle_matrix[group1, group2] <- value
}

# Convert to data frame for easier export
upper_triangle_df <- as.data.frame(upper_triangle_matrix)

# Export to CSV
write.csv(upper_triangle_df, "upper_triangle_avg_margin.csv", row.names = TRUE)
