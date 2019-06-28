# Replication code for "The impacts of the Amazon Soy Moratorium on deforestation"

Authors: Analysis was designed by Robert Heilmayr, Holly Gibbs, Lisa Rausch 
and Jacob Munger. Code was written by Robert Heilmayr.

Data: Replication data is available at ...

Depdendencies: Code was written in a combination of Python 3 and Stata 14.
Required python packages include pandas, numpy, matplotlib and seaborn.

Summary of scripts: Code should be run in the following order:
1) data_prep.py: Converts wide dataset into long formats for statistical analysis.
2) regressions.do: Runs statistical analyses and generates latex tables
3) results.ipynb: Generates figures and descriptive statistics.