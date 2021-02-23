#!/usr/bin/Rscript 
library(ggplot2)
library(patchwork)
library(rlist)
theme_set(
  theme_bw() +
    theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white")) +
    theme(legend.position = "top")
)

# R version 3.5.2

#########################################################################################################################################
# This script gets poly T tables and generates box plots representing the distribution of the percents of each poly T out of all poly Ts.
#########################################################################################################################################

# A function for decimal representation
scaleFUN <- function(x) sprintf("%.2f", x)

table_path <- "poly_T/tables/counting"
file_name <- "genes_seq_counting.csv"
figure_path <- "figures"
poly_Ts <<- c("T_1", "T_2", "T_3", "T_4", "T_5", "T_6", "T_7", "T_8", "T_9", "T_10")

adjust_p_value <- function(poly_T_data, alternate) {
  "This function gets data of 2 populations and compare them by using Wilcoxon test according to the givven alternative hypothesis"
  customSymnum <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05) / (nrow(poly_T_data) / 2), 1), symbols = c("**", "****", "***", "*", "ns")) # Define the p values symbols based on Bonferroni correction.

  transcribed_data <- subset(poly_T_data, Strand == "transcribed")$Value # Define the group of the transcribed strand.
  non_transcribed_data <- subset(poly_T_data, Strand == "non_transcribed")$Value # "" non transcribed strand.
  wilcoxon_stats <- wilcox.test(transcribed_data, non_transcribed_data, paired = T)
  p_value <- wilcoxon_stats$p.value
  # Scan the customSymnum list in order to adjust the matching symbol to the calculated p value.
  for (i in 1:(length(customSymnum$cutpoints) - 1)) {
    if (customSymnum$cutpoints[i] <= p_value && p_value < customSymnum$cutpoints[i + 1]) {
      symbol <- customSymnum$symbols[i]
    }
  }
  return(symbol)
}


make_table <- function(table, nuc) {
  "This function gets counting table and a nucleotide and generate a table with all the relative frequencies of the poly nucleotide."

  counting_table <- table[, 1:2] # Extract the Id and the length of the sequence.
  tmp_table <- table[, grepl(nuc, names(table))] # Intialize a table with the values of the given nucleotides.
  tmp_table <- tmp_table[, -ncol(tmp_table)] # Remove the "N_R" column.
  counting_table <- cbind(counting_table, table[[paste0(nuc, "_NR")]], tmp_table)
  names(counting_table) <- c("Id", "length", "Ts_num", poly_Ts)
  counting_table[poly_Ts] <- counting_table[poly_Ts] / counting_table$Ts_num * 100 # Divide each value in the total number of Ts of each interval and multi[le by 100 (in order to get percentage).

  return(counting_table)
}

prepare_data <- function(table) {
  "This function prepare the given table for plotting using ggplot2 package."

  ggplot_table <- data.frame(matrix(ncol = 2, nrow = 0)) # Initalize a data frame for the further ggplot use.

  # Scan all the counting values in the table and "tag" them by their column name.
  for (col in colnames(table[-c(1, 2, 3)])) {
    tmp <- cbind(rep(col, length(table[[col]])), as.data.frame(table[[col]])) # Create a column of poly T values and their category (i.e. T_1, T2 etc.)
    ggplot_table <- rbind(ggplot_table, tmp) # Concatenate the new columns into one data frame.
  }
  names <- c("poly_T", "value")
  colnames(ggplot_table) <- names # Adjust the table names.

  return(ggplot_table)
}

prepare_data_ggplot_list <- function(transcribed_table, non_transcribed_table) {
  "This function preprare data for ggplot of multiple box plots."

  poly_Ts_list <- replicate(n = 10, expr = {
    data.frame(matrix(ncol = 2, nrow = 0))
  }, simplify = F) # Initialize a list of 10 data frames.
  table_names <- c("Value", "Strand")

  # Scan the given tables and merge the transcribed and non values of each poly T into the same data frame.
  for (i in 1:length(poly_Ts))
  {
    tmp_trans <- as.data.frame(subset(transcribed_table["value"], transcribed_table["poly_T"] == poly_Ts[i])) # Extract the matching values of the given poly T.
    tmp_trans <- cbind(tmp_trans, rep("transcribed", nrow(tmp_trans))) # Tage the values as the transcribed strand.
    names(tmp_trans) <- table_names
    tmp_non_trans <- as.data.frame(subset(non_transcribed_table["value"], non_transcribed_table["poly_T"] == poly_Ts[i])) # Extract the matching values of the given poly T.
    tmp_non_trans <- cbind(tmp_non_trans, rep("non_transcribed", nrow(tmp_non_trans))) # Tage the values as the non-transcribed strand.
    names(tmp_non_trans) <- table_names
    poly_Ts_list[[i]] <- rbind(tmp_trans, tmp_non_trans)
  }
  return(poly_Ts_list)
}

make_box_plot_facets <- function(poly_Ts_list, figure_name) {
  "This function gets a list of data frames and generate box plots for each data frame."
  ggplots <- list() # Intialize an empty list for ggplot elements.

  yLimBreaks <- list()
  # Scan the list of the data frames and build a ggplot element for every data frame.
  for (i in 1:length(poly_Ts_list)) {
    symbol <- adjust_p_value(poly_Ts_list[[i]], "less")
    ranges <- boxplot.stats(poly_Ts_list[[i]]$Value)$stats[c(1, 5)]
    # Set the y limities based on the values without the outliers
    yLim <- ranges + c(-ranges[1] * 0.2, ranges[2] * 0.2)
    # Set the axis break values for 4 breaks.
    breaksY <- c(yLim[1], (yLim[2] - yLim[1]) / 3 + yLim[1], ((yLim[2] - yLim[1]) / 3) * 2 + yLim[1], yLim[2])
    # breaksY <- round(breaksY, 2)
    yLimBreaks <- list.append(yLimBreaks, breaksY)
    p <- ggplot(poly_Ts_list[[i]], aes(x = Strand, y = Value, group = Strand, fill = Strand)) +
      geom_boxplot(colour = "black", outlier.shape = NA) + theme(legend.position = "none") +
      xlab("") + ylab("") +
      ggtitle(bquote((T)[.(i)])) + # Extract the "poly number".
      scale_y_continuous(limits = yLim, breaks = breaksY, labels = scaleFUN) +
      stat_summary(fun = mean, geom = "point", shape = 18) +
      scale_fill_manual(values = c("goldenrod", "darkgreen")) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20)) +
      annotate("text", label = symbol, x = 1.5, y = Inf, vjust = 2, size = 6)
    name <- paste("p", i, sep = "_") # Define the name of the element.
    ggplots[[name]] <- p # Add the element to the list.
  }
  pdf(paste(figure_path, paste0(figure_name, "_poly_T_violin_plots_facets.pdf"), sep = "/"), height = 8, width = 12)
  print((ggplots$p_1 | ggplots$p_2 | ggplots$p_3 | ggplots$p_4 | ggplots$p_5) / (ggplots$p_6 | ggplots$p_7 | ggplots$p_8 | ggplots$p_9 | ggplots$p_10))
  dev.off()
}



counting_table <- read.csv(paste(table_path, file_name, sep = "/"))

figure_name <- gsub("_counting.csv", "", file_name)
# Generate tables
transcribed_counting_table <- make_table(counting_table, "A")
non_transcribed_counting_table <- make_table(counting_table, "T")
# Prepare tables for ggplot
transcribed_table_plot <- prepare_data(transcribed_counting_table)
non_transcribed_table_plot <- prepare_data(non_transcribed_counting_table)
# Prepare data for ggplot facets plot.
ggplot_data_list <- prepare_data_ggplot_list(transcribed_table_plot, non_transcribed_table_plot)
# Generate box plots of comparison all poly Ts on the transcribed and non-transcribed strands.
make_box_plot_facets(ggplot_data_list, figure_name)

