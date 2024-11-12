# Video Games Sales Analysis

This repository provides a statistical and data analysis of a dataset of video game sales using R. The analysis examines sales data for games with over 100,000 copies sold, ranked historically up to 2016. This work includes both an in-depth report and the R code used to generate it.

## Project Overview

The project includes two main files:

1. **Report_vg_sales.pdf** - A comprehensive report covering the statistical and data analysis of the video game sales dataset. The report delves into univariate and multivariate analyses, as well as clustering methods to explore patterns and trends within the data.
2. **vg_sales_r.Rmd** - The complete R code used for the analysis, created in RStudio. This file includes all steps for data cleaning, visualization, and analysis, facilitating full reproducibility.

## Dataset Summary

The dataset used in this analysis is sourced from Kaggle (https://www.kaggle.com/) and was originally scraped from [VGChartz](http://www.vgchartz.com/). It contains data on video games with over 100,000 copies sold and provides insights into various sales regions worldwide until 2016.

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

- **Univariate Analysis**: Descriptive statistics for each field, with details on key attributes such as `Platform`, `Genre`, `Publisher`, and sales in different regions.
- **Multivariate Analysis**: Examines relationships between fields, including Principal Component Analysis (PCA) for dimensionality reduction and insights on variance.
- **Clustering Analysis**: Includes both hierarchical and partitioning methods, utilizing Euclidean and Manhattan distances, to identify groups and patterns among games based on sales data. Evaluation of cluster validity is also performed.

## Code and Reproducibility

The analysis was conducted in R, and the **vg_sales_r.Rmd** file contains the complete code to reproduce each step. Key methods include data cleaning, statistical summaries, PCA, and various clustering techniques. Users can run this code in RStudio to replicate the study or adapt it for further exploration.

## Usage

To explore or reproduce this analysis:
1. Download or clone the repository.
2. Open **vg_sales_r.Rmd** in RStudio and execute the code blocks to recreate the analyses.
3. Review **Report_vg_sales.pdf** for an in-depth explanation of the findings and visualizations.

## Dataset Access

The dataset can be found on [Kaggle's Video Game Sales Dataset](https://www.kaggle.com/) page.
