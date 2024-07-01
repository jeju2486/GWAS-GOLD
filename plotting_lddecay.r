library(ggplot2)
library(dplyr)
library(svglite)
library(optparse)

options(bitmapType = 'cairo')

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "input file directory", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL, 
              help = "output file directory", metavar = "character")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if input and output directories are provided
if (is.null(opt$input) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Both input and output directories must be supplied", call. = FALSE)
}

input_file <- opt$input
output_dir <- opt$output

# Load data
print("Starting data loading...")
ld_data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Calculate distances
maxBP <- max(ld_data$BP_B)
print("Calculating distances...")
ld_data$distance <- apply(ld_data, 1, function(x) {
  direct_distance <- abs(as.numeric(x['BP_B']) - as.numeric(x['BP_A']))
  circular_distance <- abs(maxBP - as.numeric(x['BP_B'])) + as.numeric(x['BP_A'])
  min(direct_distance, circular_distance)
})
print("Distance calculation completed.")

# Aggregate data
print("Aggregating data...")
interval_size_very_short <- 5
interval_size_short <- 10
interval_size_long <- 100

ld_data_summary <- ld_data %>%
  mutate(distance_interval = case_when(
    distance <= 20 ~ floor(distance / interval_size_very_short) * interval_size_very_short,
    distance <= 100 ~ floor(distance / interval_size_short) * interval_size_short,
    TRUE ~ floor(distance / interval_size_long) * interval_size_long
  )) %>%
  group_by(distance_interval) %>%
  summarise(average_R2 = mean(R2, na.rm = TRUE))
print("Aggregation completed.")

# Save aggregated data
print("Data after aggregation:")
subset_ld_data <- ld_data_summary %>%
  filter(distance_interval < 100000)
print(head(subset_ld_data))

write.csv(subset_ld_data, file.path(output_dir, "distance_zoom.tsv"), row.names = FALSE)

# Calculate average R2 after distance of 100,000
average_R2_after_100000 <- ld_data_summary %>%
  filter(distance_interval > 100000) %>%
  summarise(average_R2_2 = mean(average_R2, na.rm = TRUE)) %>%
  pull(average_R2_2)
print(paste("Average R2 after 100000 distance:", average_R2_after_100000))

# Remove rows with invalid values for the log model
cleaned_subset_ld_data <- subset_ld_data %>%
  filter(!is.na(average_R2) & average_R2 > 0 & !is.na(distance_interval) & distance_interval > 0)

# Plotting
print("Starting plotting...")

custom_breaks <- seq(0, 100000, by = 10000)
custom_log_breaks <- c(0, 10^(seq(0, 5, by = 1)))

# Cut first 20% for log trend line
num_rows <- nrow(cleaned_subset_ld_data)
trend_data <- cleaned_subset_ld_data[1:(num_rows / 10), ]

# Calculate the slope and y-intercept of the log-LD decaying trend-line
log_model <- lm(average_R2 ~ log10(distance_interval), data = trend_data)
slope <- coef(log_model)[2]
intercept <- coef(log_model)[1]
print(paste("Slope of the log-LD decaying trend-line:", slope))
print(paste("Y-intercept of the log-LD decaying trend-line:", intercept))

# Calculate recommended LD length
recommended_LD_length <- 10^((average_R2_after_100000 - intercept) / slope)
print(paste("Recommended LD length is:", recommended_LD_length))

# Write statistics to file
stats_file <- file.path(output_dir, "statistic.txt")
writeLines(c(
  paste("Average R2 after 100,000 distance:", average_R2_after_100000),
  paste("Slope of the log-LD decaying trend-line:", slope),
  paste("Y-intercept of the log-LD decaying trend-line:", intercept),
  paste("Recommended LD length is:", recommended_LD_length)
), con = stats_file)

# Write recommended LD length to a separate file
ld_length_file <- file.path(output_dir, "recommended_ld_length.txt")
writeLines(paste("Recommended LD length is:", recommended_LD_length), con = ld_length_file)

print("Statistics and recommended LD length saved successfully.")

p <- ggplot(cleaned_subset_ld_data, aes(x = distance_interval, y = average_R2)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", formula = y ~ x, color = "red") +
  scale_x_continuous(breaks = custom_breaks, labels = custom_breaks) +
  xlab("Distance (bp)") +
  ylab("Average R^2") +
  ggtitle("LD Decay") +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 1.0, by = 0.1))

p2 <- ggplot(cleaned_subset_ld_data, aes(x = distance_interval, y = average_R2)) +
  geom_point(size = 2) +
  geom_smooth(data = trend_data, method = "lm", formula = y ~ x, color = "red") +
  scale_x_continuous(trans = 'log10', breaks = custom_log_breaks, labels = custom_log_breaks) +
  xlab("Distance (bp)") +
  ylab("Average R^2") +
  ggtitle("LD Decay") +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 1.0, by = 0.1))

# Save the plots
ggsave(file.path(output_dir, "ld_decay.png"), plot = p, width = 10, height = 6, device = "png")
ggsave(file.path(output_dir, "ld_decay_log.png"), plot = p2, width = 10, height = 6, device = "png")
ggsave(file.path(output_dir, "ld_decay.svg"), plot = p, width = 10, height = 6, device = "svg")
ggsave(file.path(output_dir, "ld_decay_log.svg"), plot = p2, width = 10, height = 6, device = "svg")
write.csv(ld_data_summary, file.path(output_dir, "distance.tsv"), row.names = FALSE)

print("Plot saved successfully.")
