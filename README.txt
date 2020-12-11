Replication code for "Brazil's Amazon Soy Moratorium reduced deforestation"

Authors: Analysis was designed by Robert Heilmayr, Holly Gibbs, Lisa Rausch 
and Jacob Munger. Code was written by Robert Heilmayr.

Data: Replication data is available at . The primary data file (wide.csv) should be downloaded and placed in a new '..\data\wide.csv' folder within the root directory of this repository.

Dependencies: Code was written in a combination of Python 3, Stata 14 and R (4.0.2).
Required python packages include pandas (v), numpy (v), matplotlib (v) and seaborn (v).
Required R packages include bife (v) and tidyverse (v).

Summary of scripts: Code should be run in the following order:
1) data_prep.py: Converts wide dataset into long formats for statistical analysis.
2) regressions.do: Runs statistical analyses and generates latex tables
3) asm_bife.R: Adds robustness results from the binomial fixed effects model.
4) results.ipynb: Generates figures and primary results from the paper.