# Video Games Sales Analysis

This repository provides a statistical and data analysis of video game sales using R. The analysis examines sales data for games with over 100,000 copies sold, ranked historically up to 2016. This project includes an in-depth report, R code, and a complete R Markdown file for full reproducibility and insights into the analysis process.

## Project Overview

The project includes three main files:

1. **Report_vg_sales.pdf** - A comprehensive report detailing the statistical and data analysis of the video game sales dataset. This report includes both univariate and multivariate analyses, as well as clustering methods to understand patterns within the data.
2. **vg_sales_r.R** - The R code used for the analysis, written in RStudio. This file includes all steps for data loading, cleaning, visualization, and analysis.
3. **vg_sales_r.Rmd** - The full R Markdown file of the report, containing both the R code and Markdown text to recreate the entire report, including code, outputs, and explanations.

## Dataset Summary

The dataset used in this analysis is sourced from Kaggle ([https://www.kaggle.com/](https://www.kaggle.com/)) and was originally scraped from [VGChartz](http://www.vgchartz.com/). It contains data on video games with over 100,000 copies sold and provides insights into various sales regions worldwide until 2016.

### Key Fields in the Dataset

- **Name** - The name of the video game
- **Platform** - Platform of the game's release (e.g., PC, PS4, etc.)
- **Year** - Year of the game's release
- **Genre** - Genre of the game
- **Publisher** - Publisher of the game
- **NA_Sales** - Sales in North America (millions)
- **EU_Sales** - Sales in Europe (millions)
- **JP_Sales** - Sales in Japan (millions)
- **Other_Sales** - Sales in the rest of the world (millions)
- **Global_Sales** - Total worldwide sales

The full dataset originally included 16,598 observations, but for optimization purposes, only the top 200 ranked games were analyzed.

## Report Highlights

The analysis in the report covers:

- **Univariate Analysis**: Descriptive statistics for each field, including attributes like `Platform`, `Genre`, `Publisher`, and regional sales.
- **Multivariate Analysis**: Relationships between fields, including Principal Component Analysis (PCA) for dimensionality reduction and insights on variance.
- **Clustering Analysis**: Both hierarchical and partitioning clustering methods, using Euclidean and Manhattan distances, to identify groups and patterns among games based on sales data. Cluster validity was also evaluated.

## Code and Reproducibility

The analysis was conducted in R, and the **vg_sales_r.R** file contains the full code to reproduce each step. The **vg_sales_r.Rmd** file provides the R Markdown version of the report, allowing you to view both code and outputs in one file. Key methods include data loading, cleaning, statistical summaries, PCA, and clustering techniques. Users can run the R Markdown file in RStudio to recreate the report or explore the analysis further.

### Requirements
Before running `vg_sales_r.R` or `vg_sales_r.Rmd`, make sure:
1. The dataset file `vg_sales.xlt` is saved in the same directory.
2. Required R packages are installed, including `readxl`.

The initial code in `vg_sales_r.R` loads the dataset as follows:

```r
# Load the necessary library
library(readxl)

# Read the Excel file, limiting to the first 200 records for analysis
vg_sales <- read_excel("vg_sales.xlt", 
    col_types = c("numeric", "text", "text", 
        "numeric", "text", "text", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric"), n_max = 200)

# Inspect the structure of the dataset
str(vg_sales)

## Usage

To explore or reproduce this analysis:

1. **Download or clone the repository.**
2. **Place `vg_sales.xlt` in the same directory** as `vg_sales_r.R` and `vg_sales_r.Rmd`.
3. **Open `vg_sales_r.R` or `vg_sales_r.Rmd` in RStudio** and execute the code blocks to recreate the analyses.
4. **Review `Report_vg_sales.pdf`** for an in-depth explanation of the findings and visualizations.

## Dataset Access

The dataset can be found on the [Kaggle Video Game Sales Dataset](https://www.kaggle.com/) page.

