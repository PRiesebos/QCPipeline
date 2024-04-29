library("ggplot2")
library("stringr")
library(scales)
library(dplyr)
library(reshape2)

# Getting Argsparser args.
args = commandArgs(trailingOnly=TRUE)
f_path = args[1]

# Creating DF.
raw_all_df <- read.csv(f_path, sep = "\t", header = TRUE)
raw_all_df$chr <- mapply(str_match, raw_all_df$var_id, pattern = ".*chr[0-9]{1,2}")
raw_all_df$call_rate <- (1 - raw_all_df$missing / raw_all_df$gt_total)

# Get unique values of var_status
var_status_levels <- unique(raw_all_df$var_status)

# Create labels for the legend
legend_labels <- paste0(var_status_levels, " (", table(raw_all_df$var_status)[var_status_levels], ")")

#Plotting
ggplot(raw_all_df, aes(x = call_rate, fill = var_status)) +
  geom_density(alpha= 0.5) +
  scale_y_continuous(labels = percent_format(accuracy = 1, scale = 1)) +
  expand_limits(x=0) +
  ggtitle("Density of Call Rates Grouped by Variant Regions\nAll Called Variants (No Filter)") +
  scale_fill_discrete(labels = legend_labels, name = "Region") +
  labs(x="Call Rate", y="Density Percentage") +
  theme(
    legend.position = c(.90, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

ggsave(
  "all_raw_call_rates.png",
  plot = last_plot(),
  scale = 1,
  width = 12,
  height = 8,
  dpi = 300
)

ggplot(raw_all_df, aes(x = alt_freq, fill = var_status)) +
  geom_density(alpha= 0.5) +
  scale_y_continuous(labels = percent_format(accuracy = 1, scale = 1)) +
  expand_limits(x=0) +
  ggtitle("Density of Alt.Freqs. Grouped by Variant Regions\nAll Called Variants (No Filter)") +
  scale_fill_discrete(labels = legend_labels, name = "Region") +
  labs(x="Alt. Freq", y="Density Percentage") +
  theme(
    legend.position = c(.90, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

ggsave(
  "all_raw_alt_freqs.png",
  plot = last_plot(),
  scale = 1,
  width = 12,
  height = 8,
  dpi = 300
)