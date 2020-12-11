Replication code for "Brazil's Amazon Soy Moratorium reduced deforestation"

## Authors
Analysis was designed by Robert Heilmayr, Holly Gibbs, Lisa Rausch 
and Jacob Munger. Code was written by Robert Heilmayr.

## Data
Replication data is available at https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/LE42B1. 
To run these scripts, the primary data file (wide.csv) should be downloaded and placed in a new '\data\' directory located inside the root directory of this repository (e.g. '..\data\wide.csv')

## Dependencies
Code was written in a combination of Python 3.6.8, Stata 14 and R 4.0.2.
Required python packages include pandas (0.24.1), numpy (1.16.2), matplotlib (3.0.3) and seaborn (0.9.0).
Required R packages include bife (0.7) and tidyverse (1.3.0).

## Summary of scripts
Code should be run in the following order:
1) data_prep.py: Converts wide dataset into long formats for statistical analysis.
2) regressions.do: Runs statistical analyses and generates latex tables
3) asm_bife.R: Adds robustness results from the binomial fixed effects model.
4) results.ipynb: Generates figures and primary results from the paper.