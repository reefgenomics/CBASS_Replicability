# Simple R script for assessing replicability of CBASS runs as a function of sample set size
# Reads data, makes scatter plots, and computes agreement metrics
# The script will read all csv files in the folder and process them as individual datasets.

# Load libraries
library(pacman)
p_load(readr, dplyr, ggplot2, ggpubr, purrr, tidyr)

# --- Clear workspace ---
rm(list = ls())

# --- Parameters --- #
# Set working directory and input/output paths. Input and output folder are relative to the working directory. # Default reads and writes all in the working directory.
work_dir <- rstudioapi::selectDirectory()   # Select working directory. if working outside RStudio, set manually. e.g.; work_dir <- "/path/to/working/directory"
if (!exists("work_dir")) work_dir <- "."    # Set working directory if it fails to set in RStudio.
out_path <- "./"                            # Default reads and writes all in the working directory.
input_folder <- "./"                        # Needs trailing "/". # Default reads and writes all in the working directory.
# --- End Parameters --- #

# --- Functions --- #
# Function that reads all the CSV in a folder (for "multiple dataset" processing)
read_data_cbass_com <- function(path_to_csv_folder) {
    csv_files <- list.files(path_to_csv_folder, pattern = "\\.csv$", full.names = TRUE)
    data_list <- lapply(csv_files, function(f) {
        dat <- read_csv(f)
        colnames(dat)[1:3] <- c("ID", "ED50_R1", "ED50_R2")
        dat
    })
    names(data_list) <- tools::file_path_sans_ext(basename(csv_files))
    return(data_list)
}

# Function: get subsets, sampling to 1000 max
get_subsets <- function(items, k, max_subsets = 1000) {
    # Calculate the total number of possible subsets of size k
    total <- choose(length(items), k)
    # If the total number of subsets exceeds max_subsets, randomly sample max_subsets subsets
    if(total > max_subsets) {
        # Generate max_subsets random subsets of size k
        replicate(max_subsets, sample(items, k), simplify = FALSE)
    } else {
        # Otherwise, generate all possible subsets of size k
        combn(items, k, simplify = FALSE)
    }
}

# Assess replicability using midpoint threshold of range
# This function evaluates agreement by:
# - Iterating over decreasing subset sizes (from n to 2)
# - For each subset size k, generating up to 1000 random subsets of k samples
# - For each subset, computing a threshold (midpoint of min/max) for ED50_R1 and ED50_R2
# - Calculating the proportion of samples where the binary classification (above/below threshold) agrees between ED50_R1 and ED50_R2
# - Aggregating the mean and standard_error of these proportions for each k
# Returns a tibble with columns: samples (k), mean_prop, sd_prop
run_predictive <- function(df) {
    n <- nrow(df)
    results <- tibble()
    for(k in seq(n, 2, by = -1)){
        subs <- get_subsets(rownames(df), k)
        props <- map_dbl(subs, function(s) {
            subdat <- df[s, ]
            th1 <- (max(subdat$ED50_R1) + min(subdat$ED50_R1)) / 2
            th2 <- (max(subdat$ED50_R2) + min(subdat$ED50_R2)) / 2
            mean((subdat$ED50_R1 > th1) == (subdat$ED50_R2 > th2))
        })
        results <- bind_rows(results,
            tibble(samples = k,
                            mean_prop = mean(props),
                            se_prop = sd(props)/sqrt(length(props))
                        )
        )
    }
    return(results)
}

# --- End Functions --- #

# --- plotting theme  --- #
theme <- theme_minimal(base_size = 14) +
    theme(
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 12),
        plot.background = element_rect(fill = "white", color = NA), # White background
        panel.background = element_rect(fill = "white", color = NA), # White panel
        panel.grid.minor = element_blank(),  # Hide minor grid lines
        panel.grid.major = element_blank(),  # Hide major grid lines
        legend.key.size = unit(1.5, "lines"),  # Increase legend key size
        panel.spacing.y = unit(1.2, "lines"),  # Increase spacing between vertical panels
        panel.border = element_blank(),        # Remove full box border
        axis.line = element_line(color = "black") # Only x and y axis lines
    )
# --- end theme  --- #

# load data
data_sets <- read_data_cbass_com(paste0(work_dir, "/", input_folder))

# Scatter plots
for(name in names(data_sets)){
    dat <- data_sets[[name]]
    p <- ggplot(dat, aes(x = ED50_R1, y = ED50_R2)) +
        geom_point() +
        theme +
        geom_smooth(method = "lm", se = FALSE) +
        labs(title = paste(name, "samples"), x = "Replicate 1", y = "Replicate 2") +
        xlim(min(dat$ED50_R1,dat$ED50_R2, na.rm = TRUE) - 0.5, max(dat$ED50_R1,dat$ED50_R2, na.rm = TRUE)+ 0.5) +
        ylim(min(dat$ED50_R1,dat$ED50_R2, na.rm = TRUE) - 0.5, max(dat$ED50_R1,dat$ED50_R2, na.rm = TRUE) + 0.5)
        
    ggsave(paste0(work_dir , "/", out_path,"scatter_", name, ".png"), p, width = 5, height = 5)
}

# calculate and plot predictive agreement
# For each dataset, run the predictive function and plot the results
predictive_list <- list()
for (name in names(data_sets)) {
    res <- run_predictive(data_sets[[name]])
    p <- ggplot(res, aes(x = samples, y = mean_prop)) +
        geom_point() +
        geom_errorbar(aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop), width = 0.2) +
        scale_x_reverse(breaks = seq(min(res$samples), max(res$samples))) +
        scale_y_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(0, 1.2)) +
        theme +
        labs(title = paste0("Replicability-", name), x = "Number of samples", y = "Proportion correct")
    predictive_list[[name]] <- p
}

# Arrange all plots in predictive_list; if less than 4, all in one row, else 1 row and 2 columns
n_plots <- length(predictive_list)
if (n_plots == 1) {
    group_plot <- ggarrange(plotlist = predictive_list, nrow = 1)
        ggsave(
            filename = paste0(work_dir, "/", out_path, "ED50_replicability.pdf"),
            plot = group_plot,
            height = 21,
            width = 14.85,
            units = "cm",
            dpi = 300
        )
    } else if (n_plots < 4) {
    group_plot <- ggarrange(plotlist = predictive_list, nrow = 1)
        ggsave(
            filename = paste0(work_dir, "/", out_path, "ED50_replicability.pdf"),
            plot = group_plot,
            height = 21,
            width = 29.7,
            units = "cm",
            dpi = 300
        )
    } else if (n_plots >= 4 && n_plots <= 6) {
    group_plot <- ggarrange(plotlist = predictive_list, nrow = 2, ncol = 3)
        ggsave(
        filename = paste0(work_dir, "/", out_path, "ED50_replicability.pdf"),
        plot = group_plot,
        height = 21,
        width = 29.7,
        units = "cm",
        dpi = 300
    )
    } else {
    group_plot <- ggarrange(plotlist = predictive_list, ncol = 3, nrow = 2)
    ggexport(
        plotlist = group_plot,
        filename = paste0(work_dir, "/", out_path, "ED50_replicability.pdf"),
        height = 21,
        width = 29.7,
        units = "cm",
        res = 300
    )
}

# --- End of script ---