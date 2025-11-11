### Scale Rows Function 

#' Scale Matrix Rows
#' @title Scale Matrix Rows
#' @description
#' Scales rows of a matrix by centering and/or scaling. Similar to base::scale()
#' but operates on rows instead of columns.
#' @param x A numeric matrix
#' @param center Logical indicating whether to center the rows (default TRUE)
#' @param scale Logical indicating whether to scale the rows (default TRUE)
#' @param add_attr Logical indicating whether to add attributes for center and scale (default TRUE)
#' @param rows Optional vector of row indices to subset
#' @param cols Optional vector of column indices to subset
#' @return A matrix with scaled rows
#' @export
#' @examples
#' # mat <- matrix(rnorm(100), nrow=10)
#' # scaled_mat <- scale_matrix_rows(mat)
scale_matrix_rows <- function(x,
                             center = TRUE,
                             scale = TRUE,
                             add_attr = TRUE,
                             rows = NULL,
                             cols = NULL) {
  
  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }
  
  # Get the row means
  cm <- rowMeans(x, na.rm = TRUE)
  
  # Get the row sd
  if (scale) {
    csd <- matrixStats::rowSds(x, center = cm, na.rm = TRUE)
  } else {
    # just divide by 1 if not
    csd <- rep(1, length = length(cm))
  }
  
  if (!center) {
    # just subtract 0
    cm <- rep(0, length = length(cm))
  }
  
  x <- (x - cm) / csd
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  
  return(x)
}

# Function: performGO
# Takes a dataframe with a gene column and performs gene ontology on the
# specified ontology
#
# Args:
#   - df: dataframe containing a column labelled gene with the hgnc gene symbol
#   - ontology = "BP", other options are "CC" for cellular component and
#                "MF" for molecular function
# Returns:
#   - dataframe containing GO results
#' Perform Gene Ontology Analysis
#' @title Gene Ontology Analysis
#' @description
#' Performs gene ontology enrichment analysis on a set of genes using clusterProfiler.
#' @param df A dataframe containing a 'Gene' column with gene symbols
#' @param ontology The ontology to use: "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
#' @param database The organism database to use (default: org.Rn.eg.db for rat)
#' @param keytype The type of gene identifier (default: "SYMBOL")
#' @return A dataframe containing GO enrichment results
#' @export
#' @examples
#' # df <- data.frame(Gene = c("ACTB", "GAPDH", "TP53"))
#' # go_results <- geneOntology(df, ontology = "BP")
geneOntology <- function(df, ontology = "BP", database = org.Rn.eg.db, keytype = "SYMBOL") {
  genes <- df$Gene
  go_analysis <- as.data.frame(clusterProfiler::enrichGO(gene = genes, OrgDb = database,
                                        keyType = keytype, ont = ontology))
  return(go_analysis)
}

# Function: performGO
# Takes a list of dataframes and performs gene ontology on each dataframe and returns
# a list of gene ontology results
#
# Args:
#   - data: list of dataframes containing each dataframe to perform GO Analysis on
#
# Returns:
#   - goresult: list containing dataframes with the go_results
#' Perform Gene Ontology on Multiple Dataframes
#' @title Perform GO Analysis on List
#' @description
#' Takes a list of dataframes and performs gene ontology analysis on each,
#' across all three ontologies (BP, MF, CC).
#' @param data A list of dataframes, each containing a 'Gene' column
#' @return A list of lists containing GO results for each dataframe and ontology
#' @export
#' @examples
#' # df1 <- data.frame(Gene = c("ACTB", "GAPDH"))
#' # df2 <- data.frame(Gene = c("TP53", "MYC"))
#' # results <- performGO(list(df1, df2))
performGO <- function(data) {
  ontologies <- c("BP", "MF", "CC")
  goresult <- list() # Empty list to store the final results
  for (i in seq_along(data)) {
    gene_ontology_data <- data[[i]]
    gene_ont <- list() # Empty list to store result for each go before merging with goresult
    n <- 1
    for (onto in ontologies) {
      gene_ont[[n]] <- geneOntology(gene_ontology_data, ontology = onto, database = org.Rn.eg.db)
      n <- n + 1    # move on to the next ontology
    }
    goresult[[i]] <- gene_ont
    }
  return(goresult)
}

# Function: performTTest
# Takes two dataframes and performs a t-test on each row of the dataframes
#
# Args:
# @param df1: dataframe containing the first cohort
# @param df2: dataframe containing the second cohort
# @param equal_var: logical, whether the t-test assumes equal variance
#
# Returns: 
# @return: dataframe containing the results of the t-test
#' Perform T-Test on Two Dataframes
#' @title Row-wise T-Test
#' @description
#' Performs a t-test on each row between two dataframes and calculates
#' fold changes and adjusted p-values.
#' @param df1 Dataframe containing the first cohort
#' @param df2 Dataframe containing the second cohort
#' @return A dataframe with the original data plus pvalue, padj, foldchange, nlog10p, and log2fc columns
#' @export
#' @examples
#' # df1 <- data.frame(s1=rnorm(10), s2=rnorm(10), s3=rnorm(10))
#' # df2 <- data.frame(s4=rnorm(10), s5=rnorm(10), s6=rnorm(10))
#' # results <- performTTest(df1, df2)
performTTest <- function(df1, df2) {
  # Ensure dataframes have the same number of rows
  if (nrow(df1) != nrow(df2)) {
    stop("df1 and df2 must have the same number of rows")
  }
  
  # Initialize result vectors
  pvalues <- numeric(nrow(df1))
  foldchange <- numeric(nrow(df1))
  
  for (i in 1:nrow(df1)) {
    # Extract values for row i (exclude NAs)
    values_df1 <- unlist(df1[i, ], use.names = FALSE)
    values_df1 <- values_df1[!is.na(values_df1)]
    
    values_df2 <- unlist(df2[i, ], use.names = FALSE)
    values_df2 <- values_df2[!is.na(values_df2)]
    
    # Check valid sample size for t-test
    if (length(values_df1) >= 3 && length(values_df2) >= 3) {
      # Compute t-test
      test_result <- tryCatch(
        t.test(values_df1, values_df2),
        error = function(e) list(p.value = 1)
      )
      
      pvalues[i] <- test_result$p.value
      foldchange[i] <- mean(values_df2, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
    } else {
      pvalues[i] <- 1
      foldchange[i] <- 0
    }
  }
  
  # Adjust p-values and compute logs
  padj <- p.adjust(pvalues, method = "BH")
  nlog10p <- -log10(pvalues)
  log2fc <- foldchange
  
  # Combine results with original data
  c_df <- cbind(df1, df2)
  c_df$pvalue <- pvalues
  c_df$padj <- padj
  c_df$foldchange <- foldchange
  c_df$nlog10p <- nlog10p
  c_df$log2fc <- log2fc
  
  return(c_df)
}

# Function: PerformMWTest
#
# Args:
# @param df1: dataframe containing the first cohort
# @param df2: dataframe containing the second cohort
#
# Returns:
# @return: dataframe containing the results of the Mann-Whitney test
#' Perform Mann-Whitney Test on Two Dataframes
#' @title Row-wise Mann-Whitney Test
#' @description
#' Performs a Mann-Whitney U test on each row between two dataframes.
#' @param df1 Dataframe containing the first cohort
#' @param df2 Dataframe containing the second cohort
#' @return A dataframe with the original data plus pvalue, padj, foldchange, and nlog10p columns
#' @export
#' @examples
#' # df1 <- data.frame(s1=rnorm(10), s2=rnorm(10), s3=rnorm(10))
#' # df2 <- data.frame(s4=rnorm(10), s5=rnorm(10), s6=rnorm(10))
#' # results <- performMWTest(df1, df2)
performMWTest <- function(df1, df2) {
  pvalues <- numeric(nrow(df1))
  foldchange <- numeric(nrow(df1))
  
  for (i in 1:nrow(df1)) {
    values_df1 <- as.numeric(df1[i, ])
    values_df2 <- as.numeric(df2[i, ])
    
    # Check if both rows have at least 4 non-NA values
    if (sum(!is.na(values_df1)) >= 3 && sum(!is.na(values_df2)) >= 3) {
      pvalues[i] <- wilcox.test(values_df1, values_df2)$p.value
      foldchange[i] <- mean(values_df2, na.rm = TRUE) / mean(values_df1, na.rm = TRUE)
    } else {
      pvalues[i] <- 1
      foldchange[i] <- 1
    }
  }
  padj <- p.adjust(pvalues, method = "BH")
  nlog10p <- -log10(pvalues)
  
  c_df <- cbind(df1, df2)
  c_df <- c_df %>%
    mutate(
      pvalue = pvalues, padj = padj, foldchange = foldchange, nlog10p = nlog10p
    )
  
  return(c_df)
}

# Function: categorizeData
# Given a dataframe and column name of Identifier columns, categorizes the data
# according to the provided classifications and returns a list of dataframes
# 
# Args: 
# @param df: numerical dataframe containing the data with proteins in rows
# and samples in columns
# @param identifier: vector, column names for identifier columns
# @param cohorts: vector, cohort identifiers to identify the columns
#
# Returns:
# @return: list of dataframes containing the categorized data
#' Categorize Data by Cohorts
#' @title Categorize Data
#' @description
#' Splits a dataframe into multiple dataframes based on cohort identifiers in column names.
#' @param df A dataframe with samples in columns
#' @param cohorts A character vector of cohort identifiers
#' @return A list of dataframes: first element contains character columns, subsequent elements contain numeric columns for each cohort
#' @export
#' @examples
#' # df <- data.frame(Gene=c("A", "B"), Control_1=c(1,2), Control_2=c(3,4), Treated_1=c(5,6))
#' # categorized <- categorizeData(df, c("Control", "Treated"))
categorizeData <- function(df, cohorts) {
  data_list <- list()
  df_char <- df %>% select(where(is.character))
  df_num <- df %>% select(where(is.numeric))
  data_list[[1]] <- df_char
  for (i in cohorts) {
    data_list[[length(data_list) + 1]] <- df_num %>% select(contains(i))
  }
  return(data_list)
}


# Function: compareCohorts
# Given a list of dataframes, performs the specified test on each pair of dataframes
#
# Args:
# @param data: list of dataframes
# @param test: character, the test to perform. Options are "ttest" or "mwtest"
#
# Returns:
# @return: list of dataframes containing the results of the test
#' Compare Cohorts with Statistical Tests
#' @title Compare Cohorts
#' @description
#' Performs pairwise statistical comparisons between all cohorts in a list.
#' @param data_list A list of dataframes containing categorized data
#' @param test The statistical test to use: "mwtest" (Mann-Whitney) or "ttest" (t-test)
#' @return A list of dataframes containing comparison results
#' @export
#' @examples
#' # data_list <- list(df_char, df_cohort1, df_cohort2)
#' # results <- compareCohorts(data_list, test = "mwtest")
compareCohorts <- function(data_list, test = "mwtest") {
  results <- list()
  
  num_data_list <- data_list[-1]
  if (test == "mwtest") {
    for (i in 1:(length(num_data_list) - 1)) {
      for (j in (i + 1):length(num_data_list)) {
        results[[length(results) + 1]] <- performMWTest(num_data_list[[i]], num_data_list[[j]])
      }
    }
  } else if (test == "ttest") {
    for (i in 1:(length(num_data_list) - 1)) {
      for (j in (i + 1):length(num_data_list)) {
        results[[length(results) + 1]] <- performTTest(num_data_list[[i]], num_data_list[[j]])
      }
    }
  } else {
    stop("Test not supported")
  }
  return(cbind(data_list[[1]], results))
}
# Function: performPCA
# Given a numeric dataframe, transpose the dataframe and perform PCA
#
# Args:
# @param df: dataframe containing numeric values
# @param sample: A vector containing the sample names
# @param scale: logical, whether to scale the data
# @param center: logical, whether to center the data
# 
# Returns:
# @return: list containing the PCA results
#' Perform Principal Component Analysis
#' @title Perform PCA
#' @description
#' Performs PCA on a dataframe and returns results with variance explained.
#' @param df A dataframe containing numeric values
#' @param sample A vector of sample names
#' @param scale Logical indicating whether to scale the data (default TRUE)
#' @param center Logical indicating whether to center the data (default TRUE)
#' @return A list containing: PCA results dataframe, x-axis label, y-axis label
#' @export
#' @examples
#' # df <- data.frame(s1=rnorm(10), s2=rnorm(10), s3=rnorm(10))
#' # pca_results <- performPCA(df, sample=c("A", "B", "C"))
performPCA <- function(df, sample, scale = TRUE, center = TRUE) {
  pca_results <- list()
  df <- df %>% select(where(is.numeric))
  pca_df <- df %>% t()
  pca_matrix <- as.matrix(pca_df)
  pca <- prcomp(pca_matrix, scale. = scale, center = center)
  pca_result_df <- as.data.frame(pca$x)
  pca_result_df$Sample <- sample
  variance_explained <- pca$sdev^2 / sum(pca$sdev^2)
  xlabel <- paste0("Principal Component 1", " (", round(variance_explained[1] * 100, 2), "%)")
  ylabel <- paste0("Principal Component 2", " (", round(variance_explained[2] * 100, 2), "%)")
  pca_results[[1]] <- pca_result_df
  pca_results[[2]] <- xlabel
  pca_results[[3]] <- ylabel
  return(pca_results)
}

# Function: plotPCA
# Given a list containing the PCA results, plot the PCA
#
# Args:
# @param pca_results: list containing the PCA results
# @param title: character, the title of the plot
#
# Returns:
# @return: ggplot2 object containing the PCA plot
#' Plot PCA Results
#' @title Plot PCA
#' @description
#' Creates a PCA plot from PCA analysis results.
#' @param pca_list A list containing PCA results from performPCA()
#' @param title The title for the plot (default: "Principal Component Analysis")
#' @return A ggplot2 object containing the PCA plot
#' @export
#' @examples
#' # pca_results <- performPCA(df, sample=c("A", "B", "C"))
#' # plot <- plotPCA(pca_results)
plotPCA <- function(pca_list,
                    title = "Principal Component Analysis") {
  pca_df <- pca_list[[1]]
  xlabel <- pca_list[[2]]
  ylabel <- pca_list[[3]]
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = Sample, color = Sample)) +
    geom_point(shape = 21, color = "black", size = 4) +
    ggprism::theme_prism(base_size = 18) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 22, face = "bold"),
          axis.text = element_text(size = 22, face = "bold"),
          axis.title.x = element_text(size = 22, face = "bold"),
          axis.minor.ticks.length = rel(0.5),
          plot.title = element_text(size = 24, face = "bold")) +
    labs(title = title, x = xlabel, y = ylabel)
  return(pca_plot)
}

# Function: makeClustermap
# Given a numeric dataframe, perform a heatmap
#
# Args:
# @param df: dataframe containing numeric values
# @param show_rownames: logical, whether to show the rownames
# @param rownames: vector, the rownames of the dataframe
# @param scale: logical, whether to scale the data
# @param col: vector containing the colors for the heatmap
# @param width: unit, the width of the heatmap
# @param height: unit, the height of the heatmap
#
# Returns:
# @return: A ComplexHeatmap object containing the heatmap
#' Create a Heatmap with Clustering
#' @title Create Clustermap
#' @description
#' Creates a heatmap using ComplexHeatmap with optional row scaling.
#' @param df A dataframe containing numeric values
#' @param scale Logical indicating whether to scale the rows (default TRUE)
#' @param show_rownames Logical indicating whether to show row names (default FALSE)
#' @param rownames Optional vector of row names to display
#' @param col Color scale for the heatmap
#' @param width Width of the heatmap
#' @param height Height of the heatmap
#' @return A ComplexHeatmap object
#' @export
#' @examples
#' # df <- data.frame(s1=rnorm(10), s2=rnorm(10), s3=rnorm(10))
#' # hm <- makeClustermap(df)
makeClustermap <- function(df, scale = TRUE, show_rownames = FALSE, rownames = NULL,
                           col = circlize::colorRamp2(c(-2, 0, 2), c("green", "black", "red")),
                           width = NULL, height = NULL) {
  df <- df %>% dplyr::select(where(is.numeric))
  # require(ComplexHeatmap) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  colnames <- colnames(df)
  if (is.null(rownames)) {
    rownames(df) <- rownames
  }
  if (scale) {
    df <- as.data.frame(scale_matrix_rows(as.matrix(df), center = TRUE, scale = TRUE))
    colnames(df) <- colnames
    df <- replace(df, is.na(df), 0)
  }
  if (show_rownames) {
    hm <- ComplexHeatmap::Heatmap(df, cluster_columns = FALSE, col = col,
                                  show_row_names = TRUE, row_labels = rownames,
                                  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                  width = width, height = height,
                                  rect_gp = gpar(col = "black"), border_gp = gpar(col = "black"),
                                  heatmap_legend_param = list(
                                    title = "Z-Score", at = c(-2, 0, 2),
                                    labels = c("-2", "0", "2"),
                                    direction = "horizontal",
                                    legend_width = width,
                                    border = TRUE,
                                    font_size = 8
                                  ))
  } else {
    hm <- ComplexHeatmap::Heatmap(df, cluster_columns = FALSE, col = col,
                                  show_row_names = FALSE,
                                  column_names_gp = gpar(fontsize = 12),
                                  width = width, height = height, border_gp = gpar(col = "black"),
                                  heatmap_legend_param = list(
                                    title = "Z-Score", at = c(-2, 0, 2),
                                    labels = c("-2", "0", "2"),
                                    direction = "horizontal",
                                    legend_width = width,
                                    border = TRUE,
                                    labels_gp = gpar(fontsize = 12, fontface = "bold")
                                  ))
  }
  return(hm)
}

# Function: makeVolcano
# Takes a dataframe with 2 columns labelled nlog10p and log2fc. It creates a new column
# in the dataframe labelled status containing information for that row (upregulated,
# downregulated, or outside statistical parameters). Then it creates and plots a volcano plot.
#
# Args:
# @param df: dataframe containing at least two columns. columns 1 must be nlog10p and
#            column 2 must be log2fc.
# @param fc_cutoff: numeric, the cutoff fold change, default is 0
# @param p_cutoff: numeric, the cutoff p-value for significance, default is 0.05

# Returns:
# @return: ggplot2 object containing the volcano plot
#' Create Volcano Plot
#' @title Make Volcano Plot
#' @description
#' Creates a volcano plot from differential expression results.
#' @param df A dataframe containing 'log2fc' and 'nlog10p' columns
#' @param fc_cutoff Fold change cutoff (default 0)
#' @param p_cutoff P-value cutoff (default 0.05)
#' @return A ggplot2 object containing the volcano plot
#' @export
#' @examples
#' # df <- data.frame(log2fc=rnorm(100), nlog10p=runif(100, 0, 5), pvalue=runif(100))
#' # plot <- makeVolcano(df, fc_cutoff=1, p_cutoff=0.05)
makeVolcano <- function(df, fc_cutoff = 0, p_cutoff = 0.05) {
  # require(ggplot2) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggprism) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  if (!("log2fc" %in% names(df))) {
    stop("Must have a columns named log2fc")
  }
  if (!("nlog10p" %in% names(df))) {
    stop("Must have a column named nlog10p")
  }
  df <- df %>% 
    mutate(
      pvalue = ifelse(is.na(pvalue), 1, pvalue),
      log2fc = ifelse(is.na(log2fc), 0, log2fc)
    )
  ymax <- ceiling(max(df$nlog10p, na.rm = TRUE) * 1.3)
  xmax <- ceiling(max(abs(df$log2fc), na.rm = TRUE) * 1.3)
  xmin <- -xmax
  x_lims <- max(xmax, abs(xmin))
  df <- df %>% 
    mutate(status = case_when(
      (nlog10p > -log10(p_cutoff)) & (log2fc >= fc_cutoff) ~ "Upregulated",
      (nlog10p > -log10(p_cutoff)) & (log2fc <= -fc_cutoff) ~ "Downregulated",
      TRUE ~ "Not Significant"
    ))
  volcanoPlot <- ggplot() +
    geom_point(data = df %>% filter(status != "Significant"),
               aes(x = log2fc, y = nlog10p, color = factor(status)),
               alpha = 0.6, shape = 19, size = 2) +
    geom_point(data = df %>% filter(status == "Significant"),
               aes(x = log2fc, y = nlog10p, color = factor(status)),
               alpha = 0.8, shape = 19, size = 3) +
    ggprism::theme_prism(base_size = 18) +
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 24, face = "bold"),
          axis.text = element_text(size = 18, face = "bold"),
          axis.title.x = element_text(size = 24, face = "bold"),
          axis.minor.ticks.length = rel(0.5),
          plot.title = element_text(size = 24, face = "bold")) +
    scale_color_manual(values = c("Upregulated" = "red2",
                                 "Downregulated" = "forestgreen",
                                 "Not Significant" = "darkgrey")) +
    labs(x = expression(log[2]*foldchange), y = expression(-log[10]*p)) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black",
               size = 0.75) 
  num_upregulated <- sum(df$status == "Upregulated")
  num_downregulated <- sum(df$status == "Downregulated")
  volcanoPlot <- volcanoPlot +
    geom_label(aes(x = (x_lims * 0.3), y = (ymax/2)), label = paste(num_upregulated),
              size = 5, fontface = "bold") +
    geom_label(aes(x = (-x_lims * 0.3), y = (ymax/2)), label = paste(num_downregulated),
              size = 5, fontface = "bold")
  if (fc_cutoff != 0) {
    volcanoPlot <- volcanoPlot + geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
                                            linetype = "dashed", color = "black",
                                            size = 0.75)
  } else {
    volcanoPlot <- volcanoPlot
  }
  return(volcanoPlot)
}

#' performANOVA4
#' @description
#' Takes in 4 dataframes, each belonging to a cohort, and performs an ANOVA on each row.
#' The function returns a dataframe with the p-value from the ANOVA, Tukey's HSD p-values for pairwise
#' comparisons, and fold changes for each pairwise comparison. Only meant to be used with compareCohorts.
#'
#'
#' @param df1 dataframe containing cohort 1
#' @param df2 dataframe containing cohort 2
#' @param df3 dataframe containing cohort 3
#' @param df4 dataframe containing cohort 4
#'
#' @return A dataframe with the results of the ANOVA
#' @export
#'
#' @examples
#' test <- performANOVA4(df1, df2, df3, df4)
performANOVA4 <- function(df1, df2, df3, df4) {
  # require(dplyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  pvalues <- numeric(nrow(df1))         # Store p-values from ANOVA
  tukey_pvalue12 <- numeric(nrow(df1))  # Store Tukey's p-values for Group1 vs Group2
  tukey_pvalue13 <- numeric(nrow(df1))  # Store Tukey's p-values for Group1 vs Group3
  tukey_pvalue14 <- numeric(nrow(df1))  # Store Tukey's p-values for Group1 vs Group4
  tukey_pvalue23 <- numeric(nrow(df1))  # Store Tukey's p-values for Group2 vs Group3
  tukey_pvalue24 <- numeric(nrow(df1))  # Store Tukey's p-values for Group2 vs Group4
  tukey_pvalue34 <- numeric(nrow(df1))  # Store Tukey's p-values for Group3 vs Group4
  foldchange12 <- numeric(nrow(df1))    # Store fold change for Group1 vs Group2
  foldchange13 <- numeric(nrow(df1))    # Store fold change for Group1 vs Group3
  foldchange14 <- numeric(nrow(df1))    # Store fold change for Group1 vs Group4
  foldchange23 <- numeric(nrow(df1))    # Store fold change for Group2 vs Group3
  foldchange24 <- numeric(nrow(df1))    # Store fold change for Group2 vs Group4
  foldchange34 <- numeric(nrow(df1))    # Store fold change for Group3 vs Group4
  tryCatch({
  for (i in 1:nrow(df1)) {
    values_df1 <- as.numeric(df1[i, ])
    values_df2 <- as.numeric(df2[i, ])
    values_df3 <- as.numeric(df3[i, ])
    values_df4 <- as.numeric(df4[i, ])

    # Combine all values from the 4 dataframes
    values <- c(values_df1, values_df2, values_df3, values_df4)
    groups <- factor(c(rep("Group1", length(values_df1)),
                       rep("Group2", length(values_df2)),
                       rep("Group3", length(values_df3)),
                       rep("Group4", length(values_df4))))

    # Ensure there are enough non-NA values in each group for ANOVA
    if (sum(!is.na(values_df1)) >= 3 && sum(!is.na(values_df2)) >= 3 && sum(!is.na(values_df3)) >= 3 && sum(!is.na(values_df4)) >= 3) {
      # Perform ANOVA
      anova_result <- aov(values ~ groups)
      pvalues[i] <- summary(anova_result)[[1]][["Pr(>F)"]][1]  # Extract p-value from ANOVA

      # Perform Tukey's HSD test if ANOVA is significant
      if (pvalues[i] <= 1) {
        tukey_result <- TukeyHSD(anova_result)

        # Extract adjusted p-values for pairwise comparisons from Tukey's test
        tukey_comparisons <- tukey_result$groups

        # Extract Tukey p-values for specific pairwise comparisons
        tukey_pvalue12[i] <- tukey_comparisons["Group2-Group1", "p adj"]
        tukey_pvalue13[i] <- tukey_comparisons["Group3-Group1", "p adj"]
        tukey_pvalue14[i] <- tukey_comparisons["Group4-Group1", "p adj"]
        tukey_pvalue23[i] <- tukey_comparisons["Group3-Group2", "p adj"]
        tukey_pvalue24[i] <- tukey_comparisons["Group4-Group2", "p adj"]
        tukey_pvalue34[i] <- tukey_comparisons["Group4-Group3", "p adj"]

        # Calculate fold changes for pairwise comparisons
        foldchange12[i] <- mean(values_df2, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
        foldchange13[i] <- mean(values_df3, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
        foldchange14[i] <- mean(values_df4, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
        foldchange23[i] <- mean(values_df3, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
        foldchange24[i] <- mean(values_df4, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
        foldchange34[i] <- mean(values_df4, na.rm = TRUE) - mean(values_df3, na.rm = TRUE)
      }
    } else {
      # Set p-values to 1 if insufficient non-NA values
      pvalues[i] <- 1
      tukey_pvalue12[i] <- 1
      tukey_pvalue13[i] <- 1
      tukey_pvalue14[i] <- 1
      tukey_pvalue23[i] <- 1
      tukey_pvalue24[i] <- 1
      tukey_pvalue34[i] <- 1
      foldchange12[i] <- 0
      foldchange13[i] <- 0
      foldchange14[i] <- 0
      foldchange23[i] <- 0
      foldchange24[i] <- 0
      foldchange34[i] <- 0
    }
  }
  }, error = function(e) {
    # Set p-values to 1 if error occurs
    pvalues[i] <- 1
    tukey_pvalue12[i] <- 1
    tukey_pvalue13[i] <- 1
    tukey_pvalue14[i] <- 1
    tukey_pvalue23[i] <- 1
    tukey_pvalue24[i] <- 1
    tukey_pvalue34[i] <- 1
    foldchange12[i] <- 0
    foldchange13[i] <- 0
    foldchange14[i] <- 0
    foldchange23[i] <- 0
    foldchange24[i] <- 0
    foldchange34[i] <- 0
  })

  # Create final dataframe including ANOVA p-values and Tukey's HSD results
  c_df <- dplyr::bind_cols(df1, df2, df3, df4)
  c_df <- c_df %>%
    dplyr::mutate(
      anova_pvalue = pvalues,
      tukey_pvalue12 = tukey_pvalue12,   # P-value for Group1 vs Group2
      tukey_pvalue13 <- tukey_pvalue13,   # P-value for Group1 vs Group3
      tukey_pvalue14 <- tukey_pvalue14,   # P-value for Group1 vs Group4
      tukey_pvalue23 <- tukey_pvalue23,   # P-value for Group2 vs Group3
      tukey_pvalue24 <- tukey_pvalue24,   # P-value for Group2 vs Group4
      tukey_pvalue34 <- tukey_pvalue34,   # P-value for Group3 vs Group4
      foldchange12 <- foldchange12,       # Fold change for Group1 vs Group2
      foldchange13 <- foldchange13,       # Fold change for Group1 vs Group3
      foldchange14 <- foldchange14,       # Fold change for Group1 vs Group4
      foldchange23 <- foldchange23,       # Fold change for Group2 vs Group3
      foldchange24 <- foldchange24,       # Fold change for Group2 vs Group4
      foldchange34 <- foldchange34        # Fold change for Group3 vs Group4
    )

  return(c_df)
}

# Function: makeBarplot
# Takes 2 dataframes, one long dataframe to be used to make the dataframe and a
# wide dataframe containing the pvalues
#
# Args:
# @param df: dataframe containing the data
# @param pvalue_df: dataframe containing the pvalues
# @param samples: vector containing sample labels (must be in order you want the barplot)
#
# Returns:
# @return: ggplot2 object containing the barplot
#' Create Barplot with Statistical Annotations
#' @title Make Barplot
#' @description
#' Creates a barplot for a specific gene with Tukey HSD p-values from ANOVA.
#' @param df A long-format dataframe with 'Treatment' and 'Abundance' columns
#' @param pvalue_df A wide-format dataframe containing Tukey p-values
#' @param cohort_labels A vector of cohort labels
#' @param gene The gene name to plot
#' @return A ggplot2 object containing the barplot
#' @export
#' @examples
#' # df_long <- data.frame(Treatment=rep(c("A","B","C","D"),each=5), Abundance=rnorm(20))
#' # pvalue_df <- data.frame(Gene="ACTB", tukey_pvalue12=0.01, tukey_pvalue13=0.05, tukey_pvalue14=0.001)
#' # plot <- makeBarplot(df_long, pvalue_df, c("A","B","C","D"), "ACTB")
makeBarplot <- function(df, pvalue_df, cohort_labels, gene) {
  # require(ggplot2) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggprism) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(dplyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(tidyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggpubr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggprism) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  
  ymax <- (ceiling(max(df$Abundance, na.rm = TRUE)/5) * 5) + 5
  ymin <- (floor(min(df$Abundance, na.rm = TRUE)/5) * 5)
  breaks <- seq(ymin, ymax, by = 5)
  
  pvalue_list <- list()
  
  # Create p-value dataframes for comparisons
  pvalue_df_1 <- data.frame(
    group1 = cohort_labels[1],
    group2 = cohort_labels[2],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.05, "*", ""))),
    y.position = ymax - 5.8
  )
  if (pvalue_df_1$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_1))
  
  pvalue_df_2 <- data.frame(
    group1 = cohort_labels[1],
    group2 = cohort_labels[3],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.05, "*", ""))),
    y.position = ymax - 3.5
  )
  if (pvalue_df_2$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_2))
  
  pvalue_df_3 <- data.frame(
    group1 = cohort_labels[1],
    group2 = cohort_labels[4],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue14),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue14) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue14) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue14) < 0.05, "*", ""))),
    y.position = ymax - 1.5
  )
  if (pvalue_df_3$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_3))
  
  # Create the plot
  p <- ggplot(df, aes(x = Treatment, y = Abundance)) +
    geom_bar(aes(fill = Treatment), stat = "summary", fun = "mean", position = "dodge", width = 0.6, alpha = 0.6, color = 'black') +
    geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.6), width = 0.2) +
    geom_jitter(aes(fill = Treatment), shape = 21, color = "black", size = 1.5, width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    ggprism::theme_prism(base_size = 8) +
    labs(title = gene, x = NULL, y = "Log2Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    scale_y_continuous(breaks = breaks, expand = c(0, 0)) +
    theme(legend.position = "none")
  
  # Add p-values to the plot
  if (length(pvalue_list) > 0) {
    for (pval_df in pvalue_list) {
      p <- p + ggprism::add_pvalue(pval_df, label = "p.adj.signif", y.position = pval_df$y.position,
                                   tip.length = 0.04, bracket.size = 0.3, label.size = 5)
    }
  }
  
  return(p)
}



#' makeBarplot3
#'
#' @param df long dataframe containing the values
#' @param pvalue_df the main wide dataframe containing the pvalues from the ANOVA result
#' @param cohort_labels vector containing the cohort labels
#' @param gene the gene to be plotted
#'
#' @return ggplot2 object containing the barplot
#' @export
#'
#' @examples
#' \dontrun{
#'   df_long <- data.frame(
#'     Treatment = rep(c("A", "B", "C"), each = 3),
#'     Abundance = rnorm(9)
#'   )
#'   tukey_df <- data.frame(
#'     Gene = "Gene1",
#'     tukey_pvalue12 = 0.03,
#'     tukey_pvalue13 = 0.2,
#'     tukey_pvalue23 = 0.04
#'   )
#'   makeBarplot3(df_long, tukey_df, c("A", "B", "C"), "Gene1")
#' }
makeBarplot3 <- function(df, pvalue_df, cohort_labels, gene) {
  # require(ggplot2) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggprism) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(dplyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(tidyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggpubr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  
  ymax <- (ceiling(max(df$Abundance, na.rm = TRUE)/5) * 5) + 5
  ymin <- (floor(min(df$Abundance, na.rm = TRUE)/5) * 5)
  breaks <- seq(ymin, ymax, by = 5)
  
  pvalue_list <- list()
  
  # Create p-value dataframes for comparisons (for 3 groups)
  pvalue_df_1 <- data.frame(
    group1 = cohort_labels[1],
    group2 = cohort_labels[2],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.05, "*", ""))),
    y.position = ymax - 5.8
  )
  if (pvalue_df_1$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_1))
  
  pvalue_df_2 <- data.frame(
    group1 = cohort_labels[1],
    group2 = cohort_labels[3],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.05, "*", ""))),
    y.position = ymax - 3.5
  )
  if (pvalue_df_2$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_2))
  
  pvalue_df_3 <- data.frame(
    group1 = cohort_labels[2],
    group2 = cohort_labels[3],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue23),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue23) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue23) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue23) < 0.05, "*", ""))),
    y.position = ymax - 1.5
  )
  if (pvalue_df_3$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_3))
  
  # Create the plot
  p <- ggplot(df, aes(x = Treatment, y = Abundance)) +
    geom_bar(aes(fill = Treatment), stat = "summary", fun = "mean", position = "dodge", width = 0.6, alpha = 0.6, color = 'black') +
    geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.6), width = 0.2) +
    geom_jitter(aes(fill = Treatment), shape = 21, color = "black", size = 1.5, width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    ggprism::theme_prism(base_size = 8) +
    labs(title = gene, x = NULL, y = "Log2Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    scale_y_continuous(breaks = breaks, expand = c(0, 0)) +
    theme(legend.position = "none")
  
  # Add p-values to the plot
  if (length(pvalue_list) > 0) {
    for (pval_df in pvalue_list) {
      p <- p + ggprism::add_pvalue(pval_df, label = "p.adj.signif", y.position = pval_df$y.position,
                                   tip.length = 0.04, bracket.size = 0.3, label.size = 5)
    }
  }
  
  return(p)
}

  

  
#' Create Simple Barplot
#' @title Make Simple Barplot
#' @description
#' Creates a simple barplot without statistical annotations.
#' @param df A long-format dataframe with 'Treatment' and 'Abundance' columns
#' @param cohort_labels A vector of cohort labels
#' @param gene The gene name to plot
#' @param num Additional space to add to y-axis maximum (default 0)
#' @return A ggplot2 object containing the barplot
#' @export
#' @examples
#' # df_long <- data.frame(Treatment=rep(c("A","B"),each=5), Abundance=rnorm(10))
#' # plot <- makeBarplot_simp(df_long, c("A","B"), "ACTB")
makeBarplot_simp <- function(df, cohort_labels, gene, num = 0) {
  # require(ggplot2) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggprism) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  
  ymax <- (ceiling(max(df$Abundance, na.rm = TRUE)/4) * 4)  + num
  ymin <- floor(min(df$Abundance, na.rm = TRUE))
  breaks <- seq(ymin, ymax, by = 5)
  
  # Create the plot
  p <- ggplot(df, aes(x = Treatment, y = Abundance)) +
    geom_bar(aes(fill = Treatment), stat = "summary", fun = "mean", position = "dodge", width = 0.6, alpha = 0.6, color = 'black') +
    geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.6), width = 0.2) +
    geom_jitter(aes(fill = Treatment), shape = 21, color = "black", size = 1.5, width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    ggprism::theme_prism(base_size = 8) +
    labs(title = gene, x = NULL, y = "Log2Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none")
  
  return(p)
}

#' calculateROC
#'
#' @description
#' Given a dataframe of predictors and a target vector,
#' this function calculates ROC curves, AUC values, confidence intervals,
#' and associated metrics for each specified predictor.
#'
#' @param df A dataframe containing predictor columns.
#' @param identifiers A character vector of column names in \code{df}
#'   for which the ROC analysis should be performed.
#' @param target_vector A numeric or factor vector containing the target values
#'   (commonly 0/1 for binary classification).
#' @param ci_level A numeric value between 0 and 1 specifying the confidence
#'   level for the AUC CI (default is 0.95).
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item \strong{identifier}: Name of the predictor column
#'     \item \strong{auc}: The AUC for that predictor
#'     \item \strong{ci_lower}: Lower bound of the AUC confidence interval
#'     \item \strong{ci_upper}: Upper bound of the AUC confidence interval
#'     \item \strong{sensitivity}: Sensitivity at various thresholds
#'     \item \strong{specificity}: 1 - Specificity at those thresholds
#'   }
#'
#' @importFrom pROC roc ci.auc
#' @importFrom dplyr bind_rows
#'
#' @examples
#' \dontrun{
#' library(pROC)
#' library(dplyr)
#'
#' # Sample data
#' set.seed(123)
#' df_example <- data.frame(
#'   Predictor1 = runif(50, min = 0, max = 1),
#'   Predictor2 = runif(50, min = 0, max = 1)
#' )
#'
#' # Binary outcome (0 or 1)
#' target_vec <- rbinom(50, size = 1, prob = 0.3)
#'
#' # Calculate ROC for both columns
#' roc_results <- calculateROC(
#'   df = df_example,
#'   identifiers = c("Predictor1", "Predictor2"),
#'   target_vector = target_vec
#' )
#'
#' head(roc_results)
#' }
#'
#' @export
calculateROC <- function(df, identifiers, target_vector, ci_level = 0.95) {
  
  # Ensure pROC is installed
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required. Please install it with install.packages('pROC').")
  }
  
  # Check for missing columns
  if (!all(identifiers %in% colnames(df))) {
    missing_cols <- setdiff(identifiers, colnames(df))
    warning(
      "The following identifiers are not in 'df': ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # Initialize an empty data frame for the combined results
  roc_df <- data.frame()
  
  for (identifier in identifiers) {
    
    # Skip if column does not exist
    if (!identifier %in% colnames(df)) {
      warning(paste("Column", identifier, "not found in dataframe. Skipping."))
      next
    }
    
    predictor <- df[[identifier]]
    
    # Check if predictor is numeric
    if (!is.numeric(predictor)) {
      warning(paste("Column", identifier, "is not numeric. Skipping ROC calculation."))
      next
    }
    
    # Calculate ROC curve
    roc_obj <- tryCatch(
      pROC::roc(response = target_vector, predictor = predictor),
      error = function(e) {
        warning(paste("Error calculating ROC for", identifier, ":", e$message))
        return(NULL)
      }
    )
    
    # If ROC calculation succeeded, extract metrics
    if (!is.null(roc_obj)) {
      auc_value <- roc_obj$auc
      ci_vals   <- pROC::ci.auc(roc_obj, conf.level = ci_level)
      
      # pROC's ci.auc typically returns c(lower, median, upper).
      ci_lower  <- ci_vals[1]
      ci_upper  <- ci_vals[3]
      std_error <- (ci_upper - ci_lower) / 3.92
      
      # Build a data frame of results for each threshold
      roc_data <- data.frame(
        identifier   = identifier,
        auc          = as.numeric(auc_value),
        ci_lower     = as.numeric(ci_lower),
        ci_upper     = as.numeric(ci_upper),
        std_error    = as.numeric(std_error),
        sensitivity  = roc_obj$sensitivities,
        specificity  = 1 - roc_obj$specificities
      )
      
      # Append to the master results
      roc_df <- dplyr::bind_rows(roc_df, roc_data)
    }
  }
  
  return(roc_df)
}


#' calculateROC_SVM_Poly
#'
#' @description
#' Fits a polynomial-kernel SVM model (via \code{e1071::svm()}) to a binary
#' target using all columns in \code{df} as predictors, then calculates ROC
#' curves, AUC, confidence intervals, and metrics similar to \code{calculateROC}.
#'
#' @param df A dataframe containing predictor columns (no target column).
#' @param target_vector A numeric or factor vector containing the binary target
#'   values (commonly 0/1).
#' @param ci_level A numeric value between 0 and 1 specifying the confidence
#'   level for the AUC CI (default is 0.95).
#'
#' @return A data frame similar in structure to \code{calculateROC}, containing:
#'   \itemize{
#'     \item \strong{identifier}: A character label (e.g., "SVM_Poly")
#'     \item \strong{auc}: The AUC for that SVM model
#'     \item \strong{ci_lower}: Lower bound of the AUC confidence interval
#'     \item \strong{ci_upper}: Upper bound of the AUC confidence interval
#'     \item \strong{std_error}: An approximate standard error of the AUC
#'     \item \strong{sensitivity}: Sensitivity at various thresholds
#'     \item \strong{specificity}: 1 - Specificity at those thresholds
#'   }
#'
#' @examples
#' \dontrun{
#' library(e1071)
#' library(pROC)
#'
#' # Example data
#' set.seed(123)
#' df_example <- data.frame(
#'   x1 = runif(50),
#'   x2 = runif(50)
#' )
#'
#' # Binary outcome (0 or 1)
#' target_vec <- rbinom(50, size = 1, prob = 0.3)
#'
#' # Fit polynomial-kernel SVM & compute ROC
#' svm_roc <- calculateROC_SVM_Poly(
#'   df = df_example,
#'   target_vector = target_vec
#' )
#'
#' head(svm_roc)
#' }
#'
#' @importFrom pROC roc ci.auc
#' @importFrom e1071 svm
#' @export
calculateSVM <- function(df, target_vector, ci_level = 0.95) {
  # Check for e1071
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("Package 'e1071' is required. Install with install.packages('e1071').")
  }
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required. Install with install.packages('pROC').")
  }
  
  # Ensure binary target
  target_factor <- as.factor(target_vector)
  if (length(levels(target_factor)) != 2) {
    stop("Target vector must have exactly 2 levels for binary classification.")
  }
  positive_class <- levels(target_factor)[2]  # We'll treat the 2nd level as "positive"
  
  # Fit polynomial-kernel SVM (all columns of df as predictors)
  svm_model <- e1071::svm(
    x = df,
    y = target_factor,
    kernel = "polynomial",
    probability = TRUE
  )
  
  # Predict probabilities
  pred_class <- predict(svm_model, df, probability = TRUE)
  # Extract the probability of the positive class
  prob_matrix <- attr(pred_class, "probabilities")
  prob_positive <- prob_matrix[, positive_class]
  
  # Calculate ROC
  roc_obj <- tryCatch(
    pROC::roc(response = target_vector, predictor = prob_positive),
    error = function(e) {
      stop("Error calculating ROC: ", e$message)
    }
  )
  
  # Extract AUC and CI
  auc_value <- roc_obj$auc
  ci_vals   <- pROC::ci.auc(roc_obj, conf.level = ci_level)
  ci_lower  <- ci_vals[1]
  ci_upper  <- ci_vals[3]
  
  # Approximate SE from the 95% CI width
  # 3.92 ~ 2 * 1.96 for a 95% CI normal approximation
  std_error <- (ci_upper - ci_lower) / 3.92
  
  # Build a data frame for each threshold in the ROC
  roc_df <- data.frame(
    identifier   = "Combined (SVM)",  # or customize as needed
    auc          = as.numeric(auc_value),
    ci_lower     = as.numeric(ci_lower),
    ci_upper     = as.numeric(ci_upper),
    std_error    = as.numeric(std_error),
    sensitivity  = roc_obj$sensitivities,
    specificity  = 1 - roc_obj$specificities
  )
  
  return(roc_df)
}

#' performANOVA6
#' @description
#' Takes in 6 data frames, each belonging to a cohort, and performs an ANOVA on each row.
#' The function returns a data frame with the p-value from the ANOVA, Tukey's HSD p-values
#' for pairwise comparisons, and fold changes for each pairwise comparison.
#'
#' @param df1 Data frame containing cohort 1
#' @param df2 Data frame containing cohort 2
#' @param df3 Data frame containing cohort 3
#' @param df4 Data frame containing cohort 4
#' @param df5 Data frame containing cohort 5
#' @param df6 Data frame containing cohort 6
#'
#' @return A data frame with the results of the ANOVA, Tukey’s HSD p-values, and fold changes
#' @export
#'
#' @examples
#' test <- performANOVA6(df1, df2, df3, df4, df5, df6)
performANOVA6 <- function(df1, df2, df3, df4, df5, df6) {
  # require(dplyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(stats) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  
  # Number of rows (assuming df1,...,df6 have the same number of rows)
  n_rows <- nrow(df1)
  
  # Store overall ANOVA p-values
  pvalues <- numeric(n_rows)
  
  # Store Tukey's p-values for all 15 pairwise comparisons
  tukey_pvalue12 <- numeric(n_rows)
  tukey_pvalue13 <- numeric(n_rows)
  tukey_pvalue14 <- numeric(n_rows)
  tukey_pvalue15 <- numeric(n_rows)
  tukey_pvalue16 <- numeric(n_rows)
  tukey_pvalue23 <- numeric(n_rows)
  tukey_pvalue24 <- numeric(n_rows)
  tukey_pvalue25 <- numeric(n_rows)
  tukey_pvalue26 <- numeric(n_rows)
  tukey_pvalue34 <- numeric(n_rows)
  tukey_pvalue35 <- numeric(n_rows)
  tukey_pvalue36 <- numeric(n_rows)
  tukey_pvalue45 <- numeric(n_rows)
  tukey_pvalue46 <- numeric(n_rows)
  tukey_pvalue56 <- numeric(n_rows)
  
  # Store fold changes for all 15 pairwise comparisons
  foldchange12 <- numeric(n_rows)
  foldchange13 <- numeric(n_rows)
  foldchange14 <- numeric(n_rows)
  foldchange15 <- numeric(n_rows)
  foldchange16 <- numeric(n_rows)
  foldchange23 <- numeric(n_rows)
  foldchange24 <- numeric(n_rows)
  foldchange25 <- numeric(n_rows)
  foldchange26 <- numeric(n_rows)
  foldchange34 <- numeric(n_rows)
  foldchange35 <- numeric(n_rows)
  foldchange36 <- numeric(n_rows)
  foldchange45 <- numeric(n_rows)
  foldchange46 <- numeric(n_rows)
  foldchange56 <- numeric(n_rows)
  
  # Loop over each row
  tryCatch({
    for (i in seq_len(n_rows)) {
      values_df1 <- as.numeric(df1[i, ])
      values_df2 <- as.numeric(df2[i, ])
      values_df3 <- as.numeric(df3[i, ])
      values_df4 <- as.numeric(df4[i, ])
      values_df5 <- as.numeric(df5[i, ])
      values_df6 <- as.numeric(df6[i, ])
      
      # Combine all values
      values <- c(values_df1, values_df2, values_df3, values_df4, values_df5, values_df6)
      groups <- factor(c(
        rep("Group1", length(values_df1)),
        rep("Group2", length(values_df2)),
        rep("Group3", length(values_df3)),
        rep("Group4", length(values_df4)),
        rep("Group5", length(values_df5)),
        rep("Group6", length(values_df6))
      ))
      
      # Check for >=3 non-NA values in each group to run ANOVA
      enough_data <- (
        sum(!is.na(values_df1)) >= 3 &&
          sum(!is.na(values_df2)) >= 3 &&
          sum(!is.na(values_df3)) >= 3 &&
          sum(!is.na(values_df4)) >= 3 &&
          sum(!is.na(values_df5)) >= 3 &&
          sum(!is.na(values_df6)) >= 3
      )
      
      if (enough_data) {
        # Perform ANOVA
        anova_result <- aov(values ~ groups)
        pvalues[i] <- summary(anova_result)[[1]][["Pr(>F)"]][1]
        
        # Perform Tukey's HSD test
        # (Checking "if (pvalues[i] <= 1)" always succeeds unless p-value is NA;
        #  presumably you'd check p < 0.05 if you only want post-hoc for significant ANOVA)
        if (!is.na(pvalues[i]) && pvalues[i] <= 1) {
          tukey_result <- TukeyHSD(anova_result)
          # Access the matrix of pairwise comparisons
          comp <- tukey_result$groups
          
          # Extract Tukey p-values for pairwise comparisons
          tukey_pvalue12[i] <- comp["Group2-Group1", "p adj"]
          tukey_pvalue13[i] <- comp["Group3-Group1", "p adj"]
          tukey_pvalue14[i] <- comp["Group4-Group1", "p adj"]
          tukey_pvalue15[i] <- comp["Group5-Group1", "p adj"]
          tukey_pvalue16[i] <- comp["Group6-Group1", "p adj"]
          tukey_pvalue23[i] <- comp["Group3-Group2", "p adj"]
          tukey_pvalue24[i] <- comp["Group4-Group2", "p adj"]
          tukey_pvalue25[i] <- comp["Group5-Group2", "p adj"]
          tukey_pvalue26[i] <- comp["Group6-Group2", "p adj"]
          tukey_pvalue34[i] <- comp["Group4-Group3", "p adj"]
          tukey_pvalue35[i] <- comp["Group5-Group3", "p adj"]
          tukey_pvalue36[i] <- comp["Group6-Group3", "p adj"]
          tukey_pvalue45[i] <- comp["Group5-Group4", "p adj"]
          tukey_pvalue46[i] <- comp["Group6-Group4", "p adj"]
          tukey_pvalue56[i] <- comp["Group6-Group5", "p adj"]
          
          # Calculate fold changes = mean(GroupX) - mean(GroupY)
          foldchange12[i] <- mean(values_df2, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange13[i] <- mean(values_df3, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange14[i] <- mean(values_df4, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange15[i] <- mean(values_df5, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange16[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange23[i] <- mean(values_df3, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
          foldchange24[i] <- mean(values_df4, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
          foldchange25[i] <- mean(values_df5, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
          foldchange26[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
          foldchange34[i] <- mean(values_df4, na.rm = TRUE) - mean(values_df3, na.rm = TRUE)
          foldchange35[i] <- mean(values_df5, na.rm = TRUE) - mean(values_df3, na.rm = TRUE)
          foldchange36[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df3, na.rm = TRUE)
          foldchange45[i] <- mean(values_df5, na.rm = TRUE) - mean(values_df4, na.rm = TRUE)
          foldchange46[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df4, na.rm = TRUE)
          foldchange56[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df5, na.rm = TRUE)
        }
      } else {
        # If not enough data for ANOVA, set p-values = 1, fold changes = 0
        pvalues[i]       <- 1
        tukey_pvalue12[i] <- 1
        tukey_pvalue13[i] <- 1
        tukey_pvalue14[i] <- 1
        tukey_pvalue15[i] <- 1
        tukey_pvalue16[i] <- 1
        tukey_pvalue23[i] <- 1
        tukey_pvalue24[i] <- 1
        tukey_pvalue25[i] <- 1
        tukey_pvalue26[i] <- 1
        tukey_pvalue34[i] <- 1
        tukey_pvalue35[i] <- 1
        tukey_pvalue36[i] <- 1
        tukey_pvalue45[i] <- 1
        tukey_pvalue46[i] <- 1
        tukey_pvalue56[i] <- 1
        
        foldchange12[i]  <- 0
        foldchange13[i]  <- 0
        foldchange14[i]  <- 0
        foldchange15[i]  <- 0
        foldchange16[i]  <- 0
        foldchange23[i]  <- 0
        foldchange24[i]  <- 0
        foldchange25[i]  <- 0
        foldchange26[i]  <- 0
        foldchange34[i]  <- 0
        foldchange35[i]  <- 0
        foldchange36[i]  <- 0
        foldchange45[i]  <- 0
        foldchange46[i]  <- 0
        foldchange56[i]  <- 0
      }
    }
  }, error = function(e) {
    # If an error occurs, fill the current row with safe defaults
    # (In a robust design, you'd handle or log the error more gracefully.)
    pvalues[i]       <- 1
    tukey_pvalue12[i] <- 1
    tukey_pvalue13[i] <- 1
    tukey_pvalue14[i] <- 1
    tukey_pvalue15[i] <- 1
    tukey_pvalue16[i] <- 1
    tukey_pvalue23[i] <- 1
    tukey_pvalue24[i] <- 1
    tukey_pvalue25[i] <- 1
    tukey_pvalue26[i] <- 1
    tukey_pvalue34[i] <- 1
    tukey_pvalue35[i] <- 1
    tukey_pvalue36[i] <- 1
    tukey_pvalue45[i] <- 1
    tukey_pvalue46[i] <- 1
    tukey_pvalue56[i] <- 1
    
    foldchange12[i]  <- 0
    foldchange13[i]  <- 0
    foldchange14[i]  <- 0
    foldchange15[i]  <- 0
    foldchange16[i]  <- 0
    foldchange23[i]  <- 0
    foldchange24[i]  <- 0
    foldchange25[i]  <- 0
    foldchange26[i]  <- 0
    foldchange34[i]  <- 0
    foldchange35[i]  <- 0
    foldchange36[i]  <- 0
    foldchange45[i]  <- 0
    foldchange46[i]  <- 0
    foldchange56[i]  <- 0
  })
  
  # Combine original data frames into one and append computed columns
  c_df <- dplyr::bind_cols(df1, df2, df3, df4, df5, df6) %>%
    dplyr::mutate(
      anova_pvalue   = pvalues,
      tukey_pvalue12 = tukey_pvalue12,  # Group1 vs Group2
      tukey_pvalue13 <- tukey_pvalue13,  # Group1 vs Group3
      tukey_pvalue14 <- tukey_pvalue14,  # Group1 vs Group4
      tukey_pvalue15 <- tukey_pvalue15,  # Group1 vs Group5
      tukey_pvalue16 <- tukey_pvalue16,  # Group1 vs Group6
      tukey_pvalue23 <- tukey_pvalue23,  # Group2 vs Group3
      tukey_pvalue24 <- tukey_pvalue24,  # Group2 vs Group4
      tukey_pvalue25 <- tukey_pvalue25,  # Group2 vs Group5
      tukey_pvalue26 <- tukey_pvalue26,  # Group2 vs Group6
      tukey_pvalue34 <- tukey_pvalue34,  # Group3 vs Group4
      tukey_pvalue35 <- tukey_pvalue35,  # Group3 vs Group5
      tukey_pvalue36 <- tukey_pvalue36,  # Group3 vs Group6
      tukey_pvalue45 <- tukey_pvalue45,  # Group4 vs Group5
      tukey_pvalue46 <- tukey_pvalue46,  # Group4 vs Group6
      tukey_pvalue56 <- tukey_pvalue56,  # Group5 vs Group6
      
      foldchange12 <- foldchange12,
      foldchange13   = foldchange13,
      foldchange14   = foldchange14,
      foldchange15   = foldchange15,
      foldchange16   = foldchange16,
      foldchange23   = foldchange23,
      foldchange24   = foldchange24,
      foldchange25   = foldchange25,
      foldchange26   = foldchange26,
      foldchange34   = foldchange34,
      foldchange35   = foldchange35,
      foldchange36   = foldchange36,
      foldchange45   = foldchange45,
      foldchange46   = foldchange46,
      foldchange56   = foldchange56
    )
  
  return(c_df)
}

#' Create Grouped Boxplot with Significance
#' @title Create Grouped Boxplot
#' @description
#' Creates a boxplot with beeswarm points and statistical significance annotations.
#' @param data A dataframe with 'group' and 'value' columns
#' @param alpha Significance threshold (default 0.05)
#' @return A ggplot2 object containing the boxplot
#' @export
#' @examples
#' # df <- data.frame(group=rep(c("A","B","C"),each=10), value=rnorm(30))
#' # plot <- create_grouped_boxplot(df, alpha=0.05)
create_grouped_boxplot <- function(data, alpha = 0.05) {
  # Check required packages
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
    stop("Please install ggbeeswarm: install.packages('ggbeeswarm')")
  }
  
  # Check required columns
  if (!all(c("group", "value") %in% colnames(data))) {
    stop("Data must contain 'group' and 'value' columns")
  }
  
  # Convert group to factor
  data$group <- factor(data$group)
  
  # Generate all possible group pairs
  all_pairs <- combn(levels(data$group), 2, simplify = FALSE)
  
  # Compute p-values for each pair
  pairwise_p <- purrr::map_dbl(all_pairs, ~{
    t.test(data$value[data$group == .x[1]], 
           data$value[data$group == .x[2]])$p.value
  })
  
  # Create pairwise data
  pairwise_data <- dplyr::tibble(
    pair = all_pairs,
    p_value = pairwise_p
  ) %>%
    dplyr::filter(p_value < alpha)
  
  # Return basic plot if no significant pairs
  if (nrow(pairwise_data) == 0) {
    warning("No significant pairs found at alpha = ", alpha)
    return(
      ggplot2::ggplot(data, ggplot2::aes(x = group, y = value, fill = group)) +
        ggplot2::geom_boxplot() +
        ggbeeswarm::geom_beeswarm() +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       legend.position = "none")
    )
  }
  
  # Create annotations
  pairwise_data <- pairwise_data %>%
    dplyr::mutate(
      annotation = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*"
      )
    )
  
  # Calculate y positions
  y_range <- range(data$value, na.rm = TRUE)
  y_start <- max(data$value, na.rm = TRUE) + 0.2
  step <- (y_range[2] - y_range[1]) * 0.1
  y_positions <- y_start + (0:(nrow(pairwise_data) - 1)) * step
  
  # Create plot
  ggplot2::ggplot(data, ggplot2::aes(x = group, y = value, fill = group)) +
    ggplot2::geom_boxplot(staplewidth = 0.5, outliers = FALSE) +
    ggbeeswarm::geom_quasirandom(aes(fill = group), shape = 21, color = "black", size = 2.2) +
    ggsignif::geom_signif(
      comparisons = pairwise_data$pair,
      annotations = pairwise_data$annotation,
      y_position = y_positions,
      tip_length = 0.01,
      map_signif_level = FALSE,
      textsize = 10,
      vjust = 0.6
    ) +
    ggplot2::ylim(c(min(data$value), max(y_positions) + step)) +
    ggprism::theme_prism(base_size = 18) +
    labs(x = NULL, y = "Log2Abundance") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      legend.position = "none"
    )
}

#' Perform Row-wise ANOVA, Extract Specific Tukey HSD Comparisons, and Calculate Log2 Fold Change
#'
#' This function iterates through rows of a combined data matrix, performs
#' ANOVA across groups defined by a grouping factor, runs Tukey's HSD post-hoc
#' test, calculates the log2 fold change (difference of means, assuming log-transformed input)
#' for user-specified comparisons, and extracts adjusted p-values.
#'
#' @param data_list A named list of numeric data frames or matrices. Each element
#'   represents a cohort (e.g., list(Q = Q_df, G = G_df, QZ = QZ_df, ...)).
#'   All data frames must have the same number of rows, corresponding to the
#'   same features (proteins) in the same order. Columns are replicates.
#'   Input data is assumed to be log-transformed (e.g., log2).
#' @param comparisons_list A character vector listing the desired comparisons
#'   in the format "Group1-Group2" (e.g., c("Q-G", "Q-QZ", "QZ-QAE")).
#'   The fold change will be calculated as mean(Group1) - mean(Group2).
#' @param feature_ids An optional vector of feature identifiers (e.g., protein IDs,
#'   gene names) corresponding to the rows of the data frames in `data_list`.
#'   If provided, these will be included in the output data frame. Must have
#'   the same length as the number of rows in the data frames.
#'
#' @return A data frame containing the results. If `feature_ids` were provided,
#'   the first column will be `FeatureID`. Subsequent columns correspond to the
#'   requested comparisons, with pairs of columns for log2 fold change (`_log2FC`)
#'   and adjusted p-value (`_p.adj`). Rows where analysis failed will have NA.
#' @export
#'
#' @examples
#' # --- Create Sample Data (simulate user's dataframes on log2 scale) ---
#' set.seed(42)
#' n_proteins <- 50
#' Q <- as.data.frame(matrix(rnorm(n_proteins * 4, mean = 10, sd = 1), nrow = n_proteins))
#' G <- as.data.frame(matrix(rnorm(n_proteins * 5, mean = 11, sd = 1), nrow = n_proteins)) # G higher mean (log2FC ~1)
#' QZ <- as.data.frame(matrix(rnorm(n_proteins * 4, mean = 10, sd = 1), nrow = n_proteins))
#' QAE <- as.data.frame(matrix(rnorm(n_proteins * 3, mean = 9, sd = 1), nrow = n_proteins)) # QAE lower mean than QZ (log2FC ~ -1)
#' GX <- as.data.frame(matrix(rnorm(n_proteins * 5, mean = 11, sd = 1), nrow = n_proteins)) # Same as G
#' GYM <- as.data.frame(matrix(rnorm(n_proteins * 4, mean = 12, sd = 1), nrow = n_proteins)) # GYM higher mean than GX (log2FC ~ -1 for GX-GYM)
#'
#' # Assign meaningful column names (optional but good practice)
#' colnames(Q) <- paste0("Q_", 1:ncol(Q))
#' colnames(G) <- paste0("G_", 1:ncol(G))
#' colnames(QZ) <- paste0("QZ_", 1:ncol(QZ))
#' colnames(QAE) <- paste0("QAE_", 1:ncol(QAE))
#' colnames(GX) <- paste0("GX_", 1:ncol(GX))
#' colnames(GYM) <- paste0("GYM_", 1:ncol(GYM))
#'
#' # Create the input list
#' data_for_analysis <- list(
#'   Q = Q, G = G, QZ = QZ, QAE = QAE, GX = GX, GYM = GYM
#' )
#'
#' # Define desired comparisons
#' my_comparisons <- c("Q-G", "Q-QZ", "QZ-QAE", "G-GX", "GX-GYM")
#'
#' # Generate dummy protein IDs
#' protein_ids <- paste0("Prot_", 1:n_proteins)
#'
#' # --- Run the function ---
#' analysis_results <- perform_rowwise_anova_tukey_fc(
#'   data_list = data_for_analysis,
#'   comparisons_list = my_comparisons,
#'   feature_ids = protein_ids
#' )
#'
#' # --- View Results ---
#' print(head(analysis_results))
#' # Note columns like Q-G_log2FC and Q-G_p.adj
#' summary(analysis_results)
#'
#' # Example: Filter for significant results with absolute log2FC > 1
#' # library(dplyr)
#' # significant_fc_Q_G <- analysis_results %>%
#' #   filter(`Q-G_p.adj` < 0.05 & abs(`Q-G_log2FC`) > 1)
#' # print(head(significant_fc_Q_G))

perform_rowwise_anova_tukey_fc <- function(data_list, comparisons_list, feature_ids = NULL) {
  
  # --- Input Validation ---
  if (!is.list(data_list) || is.data.frame(data_list)) {
    stop("Error: 'data_list' must be a named list of data frames or matrices.")
  }
  if (length(names(data_list)) != length(data_list) || any(names(data_list) == "")) {
    stop("Error: 'data_list' must be a *named* list.")
  }
  n_rows <- unique(sapply(data_list, nrow))
  if (length(n_rows) != 1) {
    stop("Error: All data frames in 'data_list' must have the same number of rows.")
  }
  # Check if data appears numeric (basic check on first data frame)
  if(!all(sapply(data_list[[1]], is.numeric))) {
    warning("Warning: Input data in the first element of 'data_list' does not appear fully numeric. Ensure all data frames contain only numeric values.")
  }
  if (!is.null(feature_ids) && length(feature_ids) != n_rows) {
    stop("Error: 'feature_ids' length must match the number of rows in the data frames.")
  }
  if(!is.character(comparisons_list)) {
    stop("Error: 'comparisons_list' must be a character vector.")
  }
  # Check comparison format
  if(!all(grepl("^[^-]+-[^-]+$", comparisons_list))) {
    stop("Error: 'comparisons_list' elements must be in the format 'Group1-Group2'.")
  }
  
  
  # --- Combine Data and Create Grouping Factor ---
  group_names <- names(data_list)
  combined_data <- do.call(cbind, data_list)
  # Ensure unique column names, necessary if original dfs had overlapping names
  colnames(combined_data) <- make.unique(unlist(sapply(data_list, colnames, USE.NAMES = FALSE)))
  
  group_vector <- factor(rep(group_names, times = sapply(data_list, ncol)))
  
  # --- Helper function to standardize comparison names (alphabetical order) ---
  # Used only for matching Tukey output, not for FC calculation direction
  standardize_comparison <- function(comp_str) {
    parts <- sort(strsplit(comp_str, "-", fixed = TRUE)[[1]])
    return(paste(parts, collapse = "-"))
  }
  
  # --- Function to process a single row ---
  process_row <- function(row_values, grouping_factor, target_comparisons_original) {
    
    # Create data frame for aov
    df_row <- data.frame(intensity = as.numeric(row_values), group = grouping_factor)
    df_row <- na.omit(df_row) # Remove NAs for this specific row
    
    # Check if enough data remains for ANOVA (at least 2 groups with >1 data point)
    group_counts <- table(df_row$group)
    valid_groups <- names(group_counts[group_counts >= 1]) # Need at least 1 point per group for mean calc
    valid_groups_anova <- names(group_counts[group_counts > 1]) # Need >1 point per group for ANOVA variance
    
    # Initialize result vector (interleaved p-val, fc)
    n_comp <- length(target_comparisons_original)
    result_vector <- rep(NA_real_, n_comp * 2)
    names(result_vector) <- character(n_comp * 2)
    
    # Calculate group means if possible
    group_means <- tapply(df_row$intensity, df_row$group, mean) # Will have NA for groups not present
    
    # --- Calculate Fold Changes ---
    fc_idx <- 2 # Start index for FC results in the interleaved vector
    for (comp in target_comparisons_original) {
      groups_to_compare <- strsplit(comp, "-", fixed = TRUE)[[1]]
      group1 <- groups_to_compare[1]
      group2 <- groups_to_compare[2]
      
      mean1 <- group_means[group1] # Returns NA if group1 not present in df_row
      mean2 <- group_means[group2] # Returns NA if group2 not present in df_row
      
      # Calculate FC only if both means are available
      if (!is.na(mean1) && !is.na(mean2)) {
        result_vector[fc_idx] <- mean1 - mean2
      } # else remains NA
      
      # Set names for FC column
      names(result_vector)[fc_idx] <- paste0(comp, "_log2FC")
      fc_idx <- fc_idx + 2 # Move to next FC slot
    }
    
    
    # --- Perform ANOVA and Tukey HSD (if possible) ---
    p_val_idx <- 1 # Start index for p-value results
    if (length(valid_groups_anova) < 2) {
      # Cannot perform ANOVA/Tukey, fill p-values with NA
      for (comp in target_comparisons_original) {
        names(result_vector)[p_val_idx] <- paste0(comp, "_p.adj")
        # result_vector[p_val_idx] is already NA
        p_val_idx <- p_val_idx + 2
      }
      return(result_vector) # Return vector with NAs for p-values and calculated FCs
    }
    
    # Try ANOVA/Tukey
    tukey_results <- tryCatch({
      aov_res <- aov(intensity ~ group, data = df_row)
      TukeyHSD(aov_res)
    }, error = function(e) {
      NULL # Return NULL if error occurs
    })
    
    # Extract Tukey p-values if successful
    if (!is.null(tukey_results) && "group" %in% names(tukey_results)) {
      tukey_table <- as.data.frame(tukey_results$group)
      tukey_comparisons_std <- sapply(rownames(tukey_table), standardize_comparison, USE.NAMES = FALSE)
      
      for (comp in target_comparisons_original) {
        req_std_comp <- standardize_comparison(comp) # Standardize requested comparison for matching
        match_idx <- which(tukey_comparisons_std == req_std_comp)
        
        if (length(match_idx) == 1) {
          result_vector[p_val_idx] <- tukey_table[match_idx, "p adj"] # Extract adjusted p-value
        } # else remains NA
        
        names(result_vector)[p_val_idx] <- paste0(comp, "_p.adj")
        p_val_idx <- p_val_idx + 2
      }
    } else {
      # ANOVA/Tukey failed, fill p-values with NA
      for (comp in target_comparisons_original) {
        names(result_vector)[p_val_idx] <- paste0(comp, "_p.adj")
        # result_vector[p_val_idx] is already NA
        p_val_idx <- p_val_idx + 2
      }
    }
    
    return(result_vector) # Return interleaved vector
    
  } # End process_row function
  
  # --- Apply process_row to each row ---
  results_matrix <- t(apply(combined_data, 1, process_row,
                            grouping_factor = group_vector,
                            target_comparisons_original = comparisons_list))
  results_df <- as.data.frame(results_matrix)
  
  # --- Add Feature IDs if provided ---
  if (!is.null(feature_ids)) {
    results_df <- cbind(FeatureID = feature_ids, results_df)
  }
  
  return(results_df)
  
} # End perform_rowwise_anova_tukey_fc function


