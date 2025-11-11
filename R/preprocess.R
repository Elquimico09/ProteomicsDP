#' Import Data from File
#' @title Import Data
#' @description
#' Imports data from various file formats (.csv, .xlsx, .xls).
#' @param filepath Path to the file to import
#' @return A dataframe containing the imported data
#' @export
#' @examples
#' # data <- importData("data.csv")
importData <- function(filepath) {
  if (grepl(".csv", filepath)) {
    data <- readr::read_csv(filepath)
  } else if (grepl(".xlsx", filepath)) {
    data <- readxl::read_excel(filepath)
  } else if (grepl(".xls", filepath)) {
    data <- readxl::read_excel(filepath)
  } else {
    stop("File type not supported")
  }
    return(data)
}

#' Filter Rows with Missing Values
#' @title Filter Missing Values
#' @description
#' Removes rows with more than a specified threshold of missing values.
#' @param df A dataframe containing the data
#' @param threshold The minimum fraction of non-missing values required (default 0.7)
#' @return A filtered dataframe
#' @export
#' @examples
#' # df <- data.frame(Gene=c("A", "B"), s1=c(1, NA), s2=c(NA, NA), s3=c(2, 3))
#' # filtered <- filterMissingValues(df, threshold=0.5)
filterMissingValues <- function(df, threshold = 0.7) {
  df_numeric <- df %>% select(where(is.numeric))
  df_character <- df %>% select(where(is.character))
  df_numeric <- df_numeric %>% mutate(coverage = rowSums(!is.na(.)) / ncol(df_numeric))
  df_merged <- cbind(df_character, df_numeric) %>% filter(coverage > threshold) %>% 
    select(-coverage)
  return(df_merged)
}

#' Filter High Abundance Proteins
#' @title Filter High Abundance
#' @description
#' Removes common high abundance proteins (albumin, immunoglobulins, etc.) from the dataset.
#' @param df A dataframe containing an 'Accession' column
#' @return A filtered dataframe with high abundance proteins removed
#' @export
#' @examples
#' # df <- data.frame(Accession=c("P02768", "Q12345"), value=c(100, 50))
#' # filtered <- filterHighAbundance(df)
filterHighAbundance <- function(df){
  high_abundant_accession <- c("P02768", "P0DOX5", "P02671",
                               "P02675", "P02679", "P02647",
                               "P02763", "P02787", "P01024",
                               "P01009", "P02766", "P01023")
  df_filtered <- df %>% dplyr::filter(!Accession %in% high_abundant_accession)
  return(df_filtered)
}

#' Filter Keratin Proteins
#' @title Filter Keratin
#' @description
#' Removes keratin proteins (common contaminants) from the dataset.
#' @param df A dataframe containing a 'Description' column
#' @return A filtered dataframe with keratin proteins removed
#' @export
#' @examples
#' # df <- data.frame(Description=c("Keratin type I", "Actin"), value=c(100, 50))
#' # filtered <- filterKeratin(df)
filterKeratin <- function(df){
  df_filtered <- df %>% dplyr::filter(!str_detect(Description, "Keratin"))
  return(df_filtered)
}

#' Extract Gene Name and Protein Name from Description column
#' @title Extract Gene Name and Protein Name
#' @description
#' Given a dataframe of PD Result with a 'Description' column, this function extracts the
#' Gene Name and Protein Name
#' @param df Data frame containing a 'Description' column
#' @return Data frame with added 'Gene' and 'Protein' columns
#' @export 
#' @examples
#' df <- data.frame(Description = c("sp|P12345|PROT1_HUMAN Protein Name 1 GN=GENE1",
#'                                  "sp|P67890|PROT2_HUMAN Protein Name 2 GN=GENE2"))
#' df_extracted <- extract_gene_protein(df)
extractGeneName <- function(df) {
  if (!"Description" %in% colnames(df)) {
    stop("Data frame must contain a 'Description' column")
  }
  
  # Extract Gene Name
  df$Gene <- sub(".* GN=([^ ]+).*", "\\1", df$Description)
  
  # Extract Protein Name by splitting on OS= and taking the first part
  df$Protein <- sub(" OS=.*", "", df$Description)
  
  return(df)
}


#' Normalize Proteomics Data
#' @title Normalize Data
#' @description
#' Normalizes data using various methods: log2, quantile, or combined approaches.
#' @param df A dataframe containing the data to normalize
#' @param method Normalization method: "log2quantile", "log2", "quantile", or "relative"
#' @return A normalized dataframe
#' @export
#' @examples
#' # df <- data.frame(Gene=c("A", "B"), s1=c(100, 200), s2=c(150, 250))
#' # normalized <- normalizeData(df, method="log2quantile")
normalizeData <- function(df, method = "log2quantile") {
  original_colnames <- colnames(df)
  df_chr <- df %>% select(where(is.character))
  df_num <- df %>% select(where(is.numeric))
  num_colnames <- colnames(df_num)
  if (method == "log2quantile") {
    df_quantile_matrix <- normalize.quantiles(as.matrix(df_num))
    df_quantile <- as.data.frame(df_quantile_matrix)
    colnames(df_quantile) <- num_colnames 
    df_log <- log2(df_quantile)
    result <- cbind(df_chr, df_log)
  } else if (method == "log2") {
    df_log <- log2(df_num + 1)
    result <- cbind(df_chr, df_log) 
  } else if (method == "quantile") {
    df_quantile_matrix <- normalize.quantiles(as.matrix(df_num))
    df_quantile <- as.data.frame(df_quantile_matrix)
    colnames(df_quantile) <- num_colnames 
    result <- cbind(df_chr, df_quantile)
  } else if (method == "relative") {
    df_relative <- as.data.frame(
      apply(df_num, 2, function(x) x / sum(x, na.rm = TRUE) * 100)
    )
    colnames(df_relative) <- num_colnames 
    result <- cbind(df_chr, df_relative)
  } else {
    stop("Method not supported")
  }
  result <- result[, original_colnames] 
  return(result)
}

#' Impute Missing Values using a Shifted Normal Distribution (MinProb)
#'
#' Imputes missing values (NA) in a numeric matrix using random draws
#' from a normal distribution tailored for each row (protein). The distribution's
#' mean is shifted down from the row's observed mean, and its standard deviation
#' is scaled down from the row's observed standard deviation. This method is
#' often used for proteomics data assuming missingness is primarily due to
#' low abundance (MNAR / left-censored).
#'
#' @param data_matrix A numeric matrix where rows represent features (e.g., proteins)
#'   and columns represent samples (replicates). Missing values should be NA.
#'   It is highly recommended to use log-transformed data as input.
#' @param shift Numeric scalar. Controls how many standard deviations the mean
#'   of the imputation distribution is shifted *down* from the observed mean
#'   for that row. Default is 1.8 (a common value used in Perseus).
#' @param scale Numeric scalar. Controls the standard deviation of the imputation
#'   distribution as a fraction of the observed standard deviation for that row.
#'   Default is 0.3 (a common value used in Perseus).
#' @param warn_rows_with_few_values Logical. If TRUE (default), prints a warning
#'   for rows where imputation could not be performed due to having fewer than 2
#'   observed values (mean/sd cannot be reliably calculated). NAs will remain
#'   in these rows.
#'
#' @return A numeric matrix with the same dimensions as `data_matrix`, where
#'   NA values have been imputed based on the described method for rows with
#'   sufficient observed data. Rows with insufficient data will retain NAs.
#'
#' @examples
#' # --- Create Sample Data (log2 scale recommended) ---
#' set.seed(123) # for reproducibility
#' mat <- matrix(rnorm(50, mean = 8, sd = 1.5), nrow = 10, ncol = 5)
#' colnames(mat) <- paste0("Sample_", 1:5)
#' rownames(mat) <- paste0("Protein_", 1:10)
#'
#' # Introduce some NAs, especially for lower abundance proteins (simulate MNAR)
#' mat[1, 1:3] <- NA  # Protein 1 low in samples 1-3
#' mat[5, 4:5] <- NA  # Protein 5 low in samples 4-5
#' mat[8, 2] <- NA   # Sporadic NA
#' mat[9, ] <- NA    # Protein 9 completely missing
#' mat[10, 1] <- rnorm(1, mean=4, sd=0.5) # Add one lower value protein
#' mat[10, 2:5] <- NA
#'
#' print("Original Matrix:")
#' print(mat)
#'
#' # --- Perform Imputation ---
#' imputed_matrix <- impute_minprob(mat, shift = 1.8, scale = 0.3)
#'
#' print("Imputed Matrix:")
#' print(imputed_matrix)
#'
#' # Check if NAs remain (should be Protein_9 and Protein_10)
#' print("Remaining NAs after imputation:")
#' print(which(is.na(imputed_matrix), arr.ind = TRUE))
#'
#' # Check imputed values are lower than observed for Protein 1
#' print("Protein 1 original observed:")
#' print(mat[1, !is.na(mat[1,])])
#' print("Protein 1 imputed values:")
#' print(imputed_matrix[1, is.na(mat[1,])])
#'
#' @export

imputeMinProb <- function(data_matrix, shift = 1.8, scale = 0.3, warn_rows_with_few_values = TRUE) {
  
  # --- Input Validation ---
  if (!is.matrix(data_matrix) || !is.numeric(data_matrix)) {
    stop("Error: 'data_matrix' must be a numeric matrix.")
  }
  if (!is.numeric(shift) || length(shift) != 1 || shift <= 0) {
    stop("Error: 'shift' must be a single positive numeric value.")
  }
  if (!is.numeric(scale) || length(scale) != 1 || scale <= 0) {
    stop("Error: 'scale' must be a single positive numeric value.")
  }
  
  # Create a copy to store results
  imputed_matrix <- data_matrix
  n_rows <- nrow(data_matrix)
  skipped_rows <- c()
  skipped_row_names <- c() # Store names for better warning message
  
  # --- Iterate through each row (protein) ---
  for (i in 1:n_rows) {
    row_data <- data_matrix[i, ]
    na_indices <- which(is.na(row_data))
    observed_values <- row_data[!is.na(row_data)]
    n_na <- length(na_indices)
    
    # Proceed only if there are missing values and enough observed values
    if (n_na > 0) {
      if (length(observed_values) >= 2) {
        # Calculate observed mean and sd
        obs_mean <- mean(observed_values)
        obs_sd <- sd(observed_values)
        
        # Handle case where sd is zero (all observed values are identical)
        # Use a very small SD instead of zero to allow rnorm to work
        if (obs_sd == 0) {
          # Use a small fraction of the mean if mean is not zero, else a tiny absolute value
          obs_sd <- if (obs_mean != 0) abs(obs_mean * 1e-6) else 1e-6
        }
        
        # Calculate parameters for the imputation distribution
        impute_mean <- obs_mean - (shift * obs_sd)
        impute_sd <- scale * obs_sd
        
        # Ensure imputation sd is positive
        if (impute_sd <= 0) {
          # Fallback to a small fraction of observed sd or tiny absolute value
          impute_sd <- if(obs_sd > 1e-6) obs_sd * 1e-3 else 1e-6
        }
        
        # Generate random numbers from the shifted distribution
        imputed_values <- rnorm(n = n_na, mean = impute_mean, sd = impute_sd)
        
        # Replace NAs in the result matrix
        imputed_matrix[i, na_indices] <- imputed_values
        
      } else {
        # Not enough observed values to calculate mean/sd reliably
        skipped_rows <- c(skipped_rows, i)
        # Store row name if available, otherwise store row index
        row_name <- rownames(data_matrix)[i]
        skipped_row_names <- c(skipped_row_names, ifelse(is.null(row_name), as.character(i), row_name))
      }
    } # End if (n_na > 0)
  } # End for loop
  
  # --- Warning for skipped rows ---
  if (length(skipped_rows) > 0 && warn_rows_with_few_values) {
    warning("Imputation skipped for rows (insufficient observed data < 2): ",
            paste(skipped_row_names, collapse = ", "))
  }
  
  return(imputed_matrix)
}