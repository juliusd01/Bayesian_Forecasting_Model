# Libraries for visualization
library(ggplot2)
library(reshape2)


Data <- load("G:/Dokumente/UNI/Goethe/1. Semester/Bayesian Modelling/Homework6/Test_2.RData")

# Initialize an empty data frame
df <- data.frame()

# Loop through each element in lgtdata
for (i in seq_along(lgtdata)) {
  # Extract the price matrix (21 rows x 3 columns)
  prices <- lgtdata[[i]]$X
  # Ensure prices is a matrix with 21 rows and 3 columns
  if (!is.matrix(prices) || nrow(prices) != 21 || ncol(prices) != 3) {
    stop(paste("Unexpected price matrix dimensions in element", i))
  }
  # Reshape prices into 7 rows (1 row per choice, 3 brands per row)
  reshaped_prices <- matrix(prices[, 1], ncol = 3, byrow = TRUE)
  # Extract choices (vector of length 7)
  choices <- lgtdata[[i]]$y
  # Ensure the length of choices matches the number of rows
  if (length(choices) != 7) {
    stop(paste("Mismatch between price rows and choices in element", i))
  }
  # Create a temporary data frame for this consumer
  temp_df <- data.frame(
    Price_A = reshaped_prices[, 1],   # Prices for Brand A
    Price_B = reshaped_prices[, 2],   # Prices for Brand B
    Price_Outside = reshaped_prices[, 3], # Prices for the outside option
    Choice = choices                   # Consumer's choice
  )
  # Add an identifier for the choice set
  temp_df$Choice_Set <- i
  # Combine with the main data frame
  df <- rbind(df, temp_df)
}


# MARGINAL DISTRIBUTION
marginal_data <- melt(df[, c("Price_A", "Price_B", "Price_Outside")], 
                      variable.name = "Brand", value.name = "Price")
# Plot marginal price distributions
ggplot(marginal_data, aes(x = Price, fill = Brand)) +
  geom_density(alpha = 0.5) +
  labs(title = "Marginal Price Distributions",
       x = "Price",
       y = "Density") +
  theme_minimal()
# Summary statistics for marginal distributions
summary_marginal <- aggregate(Price ~ Brand, data = marginal_data, summary)
print(summary_marginal)




# Task 2:
choice_counts <- table(df$Chosen_Brand)
print(choice_counts)

# Calculate the proportions
choice_proportions <- prop.table(choice_counts)
print(choice_proportions)

# Prevalence of consumers that always choose same brand
library(dplyr)
# Group by Individual_ID and check if all choices are the same
result <- df %>%
  group_by(Choice_Set) %>%
  summarize(
    Always_Brand_A = n_distinct(Chosen_Brand) == 1 && Chosen_Brand=='Brand_A',
    Always_Brand_B = n_distinct(Chosen_Brand) == 1 && Chosen_Brand == 'Brand_B',
    Always_Outside_Option = n_distinct(Chosen_Brand) == 1 && Chosen_Brand=='Brand_Outside'
  )


# Count how many individuals always choose the same brand
total_always_Brand_A <- sum(result$Always_Brand_A)
total_always_Brand_B <- sum(result$Always_Brand_B)
total_always_Brand_Out <- sum(result$Always_Outside_Option)

# Print the result
cat("Number of individuals who always choose Brand A: ", total_always_Brand_A, "\n")
cat("Number of individuals who always choose Brand B: ", total_always_Brand_B, "\n")
cat("Number of individuals who always choose Outside Option: ", total_always_Brand_Out, "\n")



                