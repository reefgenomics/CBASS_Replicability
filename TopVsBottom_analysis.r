# Simple R script for analysis of CBASS ED50 data reproducibility
# Reads data, makes scatter plots, and computes predictive agreement metrics
# --- Setup --- #
# setwd("/Users/chris/Dropbox/job/ED50 ranks/")
library(pacman)
p_load(readr, Hmisc, data.table, ggplot2, ggpubr, tidyr)

# --- Parameters --- #
out_path <- "/home/colinl/Proj/202505_combinations_CV/plots/" # Output path for plots
input_folder <- "Input_csv_test/" # Needs trail "/" #TODO: make full path based on wd or out
take <- "NEW_TEST_LC"          # Was: "1 or 2: which dataset to process" now used for specifing a prefix. will work also empty ""
version <- "luigi"           # "ben" or "dan" for overlap metric calculation # luigi for Percentage based version. 

#TODO: join, add n_top and bottom as variable. cleanup names and layout.


# --- Data Loading and Preparation --- #
# Read data
# Function that reads all the CSV in a folder (for "multiple dataset" processing)
read_data_cbass_com <- function(path_to_csv_folder) {
    csv_files <- list.files(path_to_csv_folder, pattern = "\\.csv$", full.names = TRUE)
    data_list <- lapply(csv_files, function(f) {
        dat <- read_csv(f)
        colnames(dat)[1:3] <- c("ID", "ED50_R1", "ED50_R2")
        dat$meanED50 <- rowMeans(dat[, 2:3])
        dat <- dat[order(dat$ID), ]
        dat
    })
    names(data_list) <- tools::file_path_sans_ext(basename(csv_files))
    return(data_list)
}


data_sets <- read_data_cbass_com(input_folder)

results_list <- list()
for (dataset in names(data_sets)) {
    # --- Initialize Output Data Frame --- #
    ResOutputall <- data.frame(
        numsamples = numeric(),
        numoverlapping = numeric(),
        Ttestpvalue = numeric(),
        Wilcoxtestpvalue = numeric()
    )

    # --- Main Simulation Loop --- #
        # --- Random Sampling and Analysis --- #
    for (i in 10:nrow(data_sets[[dataset]])) {
        for (j in 1:1000) {
            # Randomly sample i rows
            x <- data_sets[[dataset]][sample(nrow(data_sets[[dataset]]), i), ]
            
            if (version == "dan") { # gets the top 5 and bottom 5.
            x$Rank1 <- rank(-x$ED50_R1)
            x$Rank2 <- rank(-x$ED50_R2)
            dt1 <- data.table(x, key = "Rank1")
            R1_top <- head(dt1, n = 5)
            R1_bot <- tail(dt1, n = 5)
            numoverlaptopbottom <- mean(c( # Calculates the mean of the two coutns below. If numoverlaptopbottom ≈ 0, the assays agree well. If numoverlaptopbottom ≈ 5, the assays disagree (genet well performing in one assay may be weak in the other).
                sum(R1_top$Rank2 > min(R1_bot$Rank2)), #  Counts how many ge top ED50_R1net have worse ED50_R2 rankings than the best ED50_R2 in the bottom ED50_R1 group.
                sum(R1_bot$Rank2 < max(R1_top$Rank2)) # Counts how many bottom ED50_R1 genet have better ED50_R2 rankings than the worst ED50_R2 in the top ED50_R1 group.
            ))
            } else if (version == "ben") { # gets the top 5 and bottom 5. adding up the number of overlapping genets between the two groups.
            R1_top <- x[order(x$ED50_R1, decreasing = TRUE)[1:5], ]
            R1_bot <- x[order(x$ED50_R1, decreasing = FALSE)[1:5], ]
            R1_top_genets <- R1_top$ID
            R1_bot_genets <- R1_bot$ID
            R2_top5_genets <- x$ID[order(x$ED50_R2, decreasing = TRUE)[1:5]]
            R2_bot5_genets <- x$ID[order(x$ED50_R2, decreasing = FALSE)[1:5]]
            numoverlaptopbottom <- mean(c(
                sum(R1_top_genets %in% R2_bot5_genets),
                sum(R1_bot_genets %in% R2_top5_genets)
            ))
        } else if (version == "luigi") {  ## Same as ben version but percentage based.
            n <- nrow(data_sets[[dataset]])
            n20 <- max(1, round((percent/100) * n))  # Ensure at least 1 sample
            R1_top <- x[order(x$ED50_R1, decreasing = TRUE)[1:n20], ]
            R1_bot <- x[order(x$ED50_R1, decreasing = FALSE)[1:n20], ]
            R1_top_genets <- R1_top$ID
            R1_bot_genets <- R1_bot$ID
            R2_top_genets <- x$ID[order(x$ED50_R2, decreasing = TRUE)[1:n20]]
            R2_bot_genets <- x$ID[order(x$ED50_R2, decreasing = FALSE)[1:n20]]
            numoverlaptopbottom <- mean(c(
                sum(R1_top_genets %in% R2_bot_genets),
                sum(R1_bot_genets %in% R2_top_genets)
            ))
        } else if (version == "luigi2") {  ## midrange based is inperfect as it dose can split the dataset very unevenly. #same as combinations_lc 
            th1 <- median(x$ED50_R1)
            th2 <- median(x$ED50_R2)
            # th1 <- (max(x$ED50_R1) + min(x$ED50_R1)) / 2
            # th2 <- (max(x$ED50_R2) + min(x$ED50_R2)) / 2
            R1_top <- x[(x$ED50_R1 > th1),]
            R1_bot <- x[(x$ED50_R1 < th1),]
            R1_top_genets <- R1_top$ID
            R1_bot_genets <- R1_bot$ID
            R2_top_genets <- x$ID[x$ED50_R1 > th2]
            R2_bot_genets <- x$ID[x$ED50_R1 < th2]
            numoverlaptopbottom <- mean(c(
                sum(R1_top_genets %in% R2_bot_genets),
                sum(R1_bot_genets %in% R2_top_genets)
            ))
        } 
        
        # Statistical tests
        TtestPvalue <- t.test(R1_top$meanED50, R1_bot$meanED50, alternative = "g")$p.value
        WilcoxtestPvalue <- wilcox.test(R1_top$meanED50, R1_bot$meanED50, alternative = "g")$p.value
        
        # Store results
        ResOutputall <- rbind(
            ResOutputall,
            data.frame(
                numsamples = i,
                numoverlapping = numoverlaptopbottom,
                Ttestpvalue = TtestPvalue,
                Wilcoxtestpvalue = WilcoxtestPvalue
            )
            )
        }
    }
    results_list[[dataset]] <- ResOutputall
}

adjusted_results_list <- list()
Output_Summary_all_list <- list()
lm_all_list <- list()
Spearman_all_list <- list()

for (dataset in names(data_sets)) {
    # --- Adjust P-values --- #
    AdjResOutputall <- data.frame()
    for (i in 10:nrow(data_sets[[dataset]])) {
        subset <- results_list[[dataset]][results_list[[dataset]]$numsamples == i, ]
        AdjResOutputall <- rbind(
            AdjResOutputall,
            cbind(
            subset,
            Ttestpadj = p.adjust(subset$Ttestpvalue, method = "fdr"),
            Wilcoxtestpadj = p.adjust(subset$Wilcoxtestpvalue, method = "fdr")
            )
        )
    }

    # --- Summarize Results --- #
    OutputSummaryall <- data.frame(
        numsamples = 10:nrow(data_sets[[dataset]]),
        Sigttestpadj = aggregate(Ttestpadj ~ numsamples, function(x) sum(x < 0.05), data = AdjResOutputall)[[2]],
        SigWilcoxpadj = aggregate(Wilcoxtestpadj ~ numsamples, function(x) sum(x < 0.05), data = AdjResOutputall)[[2]],
        meanoverlap = aggregate(numoverlapping ~ numsamples, mean, data = AdjResOutputall)[[2]],
        minoverlap = aggregate(numoverlapping ~ numsamples, min, data = AdjResOutputall)[[2]],
        maxoverlap = aggregate(numoverlapping ~ numsamples, max, data = AdjResOutputall)[[2]],
        medianoverlap = aggregate(numoverlapping ~ numsamples, median, data = AdjResOutputall)[[2]],
        stddevoverlap = aggregate(numoverlapping ~ numsamples, sd, data = AdjResOutputall)[[2]]
    )
# --- Correlation and Regression --- #
Spearman.all <- cor.test(data_sets[[dataset]]$ED50_R1, data_sets[[dataset]]$ED50_R2, method = "spearman")
lm.all <- lm(ED50_R2 ~ ED50_R1, data = data_sets[[dataset]])

# --- Save Results to file list --- #
adjusted_results_list[[dataset]] <- AdjResOutputall
Output_Summary_all_list[[dataset]] <- OutputSummaryall
lm_all_list[[dataset]] <- lm.all
Spearman_all_list[[dataset]] <- Spearman.all
}

# --- plotting theme  --- #
        theme <- theme_minimal() +
            theme(
                strip.text.y = element_blank(),
                strip.text.x = element_text(size = 12),
                panel.grid.minor = element_blank(),  # Hide minor grid lines
                panel.grid.major = element_blank(),  # Hide minor grid lines
                legend.key.size = unit(1.5, "lines"),  # Increase legend key size
                panel.spacing.y = unit(1.2, "lines"),  # Increase spacing between vertical panels
                panel.border = element_blank(),        # Remove full box border
                axis.line = element_line(color = "black") # Only x and y axis lines
            )


# --- Plotting LC --- #
for (dataset in names(data_sets)) {
    # 1. Scatter plot of ED50_R1 vs ED50_R2 with regression line and annotations
    p1 <- ggplot(data_sets[[dataset]], aes(x = ED50_R1, y = ED50_R2)) +
        geom_point() +
        geom_smooth(method = "lm", color = "red", se = FALSE, linewidth = 1.2) +
        theme + 
        labs(
            x = "ED50_run1 °C",
            y = "ED50_run2 °C",
            title = paste0(
                "Take ", take, ": Correlation R1 vs. R2 ED50s by ID",
                if (take == 2 && remove.outlier) ". 1 outlier removed!" else ""
            )
        ) +
        xlim(min(data_sets[[dataset]]$ED50_R1) - 0.2, max(data_sets[[dataset]]$ED50_R1) + 0.2) +
        ylim(min(data_sets[[dataset]]$ED50_R2) - 0.2, max(data_sets[[dataset]]$ED50_R2) + 0.2) +
        annotate(
            "text",
            x = max(data_sets[[dataset]]$ED50_R1),
            y = min(data_sets[[dataset]]$ED50_R2) + 0.3,
            hjust = 1,
            label = paste("R-squared:", round(summary(lm.all)$r.squared, 2))
        ) +
        annotate(
            "text",
            x = max(data_sets[[dataset]]$ED50_R1),
            y = min(data_sets[[dataset]]$ED50_R2) + 0.6,
            hjust = 1,
            label = paste("Spearman's Rho:", round(Spearman.all$estimate, 2))
        )

    # 2. Line plot of number of significant comparisons vs. number of samples
    p2 <- ggplot(Output_Summary_all_list[[dataset]], aes(x = numsamples, y = Sigttestpadj)) +
        geom_line() +
        geom_point() +
        theme + 
        labs(
            x = "Number of samples",
            y = "Number of significant comparisons",
            title = ifelse(version == "luigi", paste0("Num Sig comparisons top ", percent,"% vs. bottom ", percent,"%  ED50s"), "Num Sig comparisons top5 vs. bottom5 ED50s")
        )

    # 3. Errorbar plot for mean overlap ± 1 stddev
    p3 <- ggplot(Output_Summary_all_list[[dataset]], aes(x = numsamples, y = meanoverlap)) +
        geom_line() +
        geom_point() +
        geom_errorbar(
            aes(
                ymin = meanoverlap - stddevoverlap,
                ymax = meanoverlap + stddevoverlap
            ),
            width = 0.2
        ) +
        theme +
        labs(
            x = "Number of samples",
            y = "Mean overlap ± 1 stddev",
            title = ifelse(version == "luigi", paste0("Num Sig comparisons top ", percent,"% vs. bottom ", percent,"%  ED50s ", version, " method"), paste0("Mean overlap top5 vs. bottom5 ED50s ", version, " method"))
        ) +
        if (version == "dan") geom_hline(yintercept = 3, linetype = "dashed", color = "blue") else NULL

    # 4. Save plots
    pdf_name <- if (take == 2 && remove.outlier) {
        paste0(out_path, take, "_", dataset,"_reproducibility_1000reps_all_", version, "_version_outlier_removed.pdf")
        } else {
        paste0(out_path, take, "_", dataset,"_reproducibility_1000reps_all_", version, "_version.pdf")
        }
    
    # # Save plot togheter
    # group_plot <- ggarrange(p1, p2, p3, ncol = 1)
    # ggsave(pdf_name, group_plot, width = 21, height = 39.7, units = "cm", dpi = 300)
    
    # Save each plot separately
    pdf(pdf_name, width = 12, height = 6)
    print(p1)
    print(p2)
    print(p3)
    dev.off()
}
# --- End of script --- #
quit() #Added to ensure clearing of session in between "version changes"
