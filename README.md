# proteomicsUtils

A collection of utility functions for proteomics data analysis.

## Description

This package provides a comprehensive set of functions for proteomics data analysis, including:

- Data import and preprocessing
- Normalization and filtering
- Statistical testing (t-tests, Mann-Whitney, ANOVA)
- Data visualization (volcano plots, heatmaps, PCA plots, barplots)
- Gene ontology analysis
- Machine learning utilities (ROC curves, SVM)
- Missing value imputation

## Installation

### From GitHub (Recommended)

Once this package is on GitHub, you can install it directly:

```r
# Install remotes if you don't have it
install.packages("remotes")

# Install from GitHub (replace 'yourusername' with your GitHub username)
remotes::install_github("yourusername/proteomicsUtils")
```

Or using devtools:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install from GitHub
devtools::install_github("yourusername/proteomicsUtils")
```

The package will automatically install most dependencies. However, you may need to manually install Bioconductor packages first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "ComplexHeatmap"))
```

### From Local Source

If you have the source files locally:

```r
install.packages("proteomicsUtils_0.1.0.tar.gz", repos = NULL, type = "source")
```

Or from terminal:

```bash
R CMD INSTALL proteomicsUtils_0.1.0.tar.gz
```

## Usage

After installation, load the package in your R session:

```r
library(proteomicsUtils)
```

Now you can use any of the 29 functions directly without sourcing the utilities.R file!

### Example Workflow

```r
library(proteomicsUtils)

# Import data
data <- importData("proteomics_data.csv")

# Filter and normalize
filtered <- filterMissingValues(data, threshold = 0.5)
normalized <- normalizeData(filtered)

# Impute missing values
imputed <- imputeValues(normalized)

# Statistical testing
results <- performTTest(imputed, group1_cols, group2_cols)

# Visualization
volcano_plot <- makeVolcano(results, p_cutoff = 1.3, fc_cutoff = 1)
pca_results <- performPCA(imputed)
pca_plot <- plotPCA(pca_results)
```

### Example Functions

- `importData()` - Import proteomics data from CSV or Excel files
- `normalizeData()` - Normalize your data
- `filterMissingValues()` - Filter rows with too many missing values
- `imputeValues()` - Impute missing values
- `performTTest()` - Perform t-tests
- `performANOVA4()` - Perform ANOVA with 4 cohorts
- `makeVolcano()` - Create volcano plots
- `performPCA()` - Perform PCA analysis
- `plotPCA()` - Create PCA plots
- `makeClustermap()` - Create heatmaps
- `geneOntology()` - Perform gene ontology analysis
- `calculateROC()` - Calculate and plot ROC curves
- `create_grouped_boxplot()` - Create grouped boxplots

## Available Functions

The package includes 29 functions:

1. makeVolcano_deprecated
2. scale_matrix_rows
3. extractGeneName
4. geneOntology
5. performGO
6. performTTest
7. performMWTest
8. importData
9. categorizeData
10. compareCohorts
11. performPCA
12. plotPCA
13. makeClustermap
14. filterHighAbundance
15. filterKeratin
16. normalizeData
17. filterMissingValues
18. makeVolcano
19. performANOVA4
20. makeBarplot
21. makeBarplot3
22. makeBarplot_simp
23. imputeValues
24. calculateROC
25. calculateSVM
26. performANOVA6
27. create_grouped_boxplot
28. impute_minprob
29. perform_rowwise_anova_tukey_fc

## License

MIT License

## Author

Vishal Sandilya  
vishal.sandilya@ttu.edu