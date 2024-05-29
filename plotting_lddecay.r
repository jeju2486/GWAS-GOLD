library(ggplot2)
library(dplyr)
library(svglite)
library(optparse)

options(bitmapType='cairo')

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

print("Starting data loading...")
ld_data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

maxBP <- max(ld_data$BP_B)
print(paste("Maximum BP_B:", maxBP))

print("Calculating distances...")
# Calculate the minimum distance for each row using circular distance calculation
ld_data$distance <- apply(ld_data, 1, function(x) {
  direct_distance = abs(as.numeric(x['BP_B']) - as.numeric(x['BP_A']))
  circular_distance = abs(maxBP - as.numeric(x['BP_B'])) + as.numeric(x['BP_A'])
  min(direct_distance, circular_distance)
})
print("Distance calculation completed.")

print("Aggregating data...")
# Define interval sizes for aggregation
interval_size_very_short <- 5
interval_size_short <- 10       # For distances from 11bp up to 1000bp
interval_size_long <- 100       # For distances beyond 1000bp

# Adjusted aggregation with additional condition
ld_data_summary <- ld_data %>%
  mutate(distance_interval = case_when(
    distance <= 20 ~ floor(distance / interval_size_very_short) * interval_size_very_short,
    distance <= 100 ~ floor(distance / interval_size_short) * interval_size_short,
    TRUE ~ floor(distance / interval_size_long) * interval_size_long
  )) %>%
  group_by(distance_interval) %>%
  summarise(average_R2 = mean(R2, na.rm = TRUE))
print("Aggregation completed.")

print("Data after aggregation:")
subset_ld_data <- ld_data_summary[ld_data_summary$distance_interval < 40000, ]
print(head(subset_ld_data))

write.csv(subset_ld_data, file.path(output_dir, "distance_zoom.tsv"), row.names = FALSE)

print("Starting plotting...")

# Define custom breaks for the plot to accurately represent the mixed intervals
custom_breaks <- c(seq(0, 40000, by = 5000))
custom_log_breaks <- c(0, 10^(seq(0, 5, by = 1)))

#cut first 20% for log trend line
num_rows <- nrow(subset_ld_data)
trend_data <- subset_ld_data[1:(num_rows/5), ]

p <- ggplot(subset_ld_data, aes(x = distance_interval, y = average_R2)) +
  geom_point(size = 2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), color = "red") +  # Use method = "loess" for smooth curve
#  geom_smooth(method = "lm", formula=y ~ poly(x, 2), color = "red") +  
  scale_x_continuous(breaks = custom_breaks, labels = custom_breaks) +
  xlab("Distance (bp)") +
  ylab("Average R^2") +
  ggtitle("LD Decay") +
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank(), 
    legend.background = element_blank() 
  ) +
  scale_y_continuous(breaks = seq(0, 1.0 , by = 0.1))
  
p2 <- ggplot(subset_ld_data, aes(x = distance_interval, y = average_R2)) +
  geom_point(size = 2) +
  geom_smooth(data = trend_data, method = "lm", formula = y ~ x, color = "red") +  # Use method = "loess" for smooth curve  
  scale_x_continuous(trans='log10',breaks = custom_log_breaks, labels = custom_log_breaks) +
  xlab("Distance (bp)") +
  ylab("Average R^2") +
  ggtitle("LD Decay") +
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank(), 
    legend.background = element_blank() 
  ) +
  scale_y_continuous(breaks = seq(0, 1.0 , by = 0.1))

# Save the plots
ggsave(file.path(output_dir, "ld_decay_scc.png"), plot = p, width = 10, height = 6, device="png")
ggsave(file.path(output_dir, "ld_decay_log_scc.png"), plot = p2, width = 10, height = 6, device="png")
ggsave(file.path(output_dir, "ld_decay_scc.svg"), plot = p, width = 10, height = 6, device = "svg")
ggsave(file.path(output_dir, "ld_decay_log_scc.svg"), plot = p2, width = 10, height = 6, device = "svg")
write.csv(ld_data_summary, file.path(output_dir, "scc_distance.tsv"), row.names = FALSE)

print("Plot saved successfully.")
