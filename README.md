Metabolites Correlation Analysis Tool for Galaxy
================================================

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)  [![Build Status](https://travis-ci.org/workflow4metabolomics/correlation_analysis.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/correlation_analysis)

Our project
-----------
The [Workflow4Metabolomics](http://workflow4metabolomics.org), W4M in short, is a French infrastructure offering software tool processing, analyzing and annotating metabolomics data. It is based on the Galaxy platform.


Metabolites Correlation Analysis
--------------------------------
This tool takes as inputs, tabular table files from the metabolomic workflow.(variableMetadata,dataMatrix and SampleMetadata) and executes two possible functions ("sorting" and "corrdel").

The second option takes as input a table of your own and execute the "corr_matrix" function.

The first function "sorting":

1- It sorts the dataframe by rtmed column.
2- It computes the mean operation of all the signal values of the metabolites by sample, and put the results in a new column "signal_moy".
3- It finally creates a tabular output "sorted_table.tsv".
    
The second function "corrdel" for each pcgroup of the previous sorted tabular file "sorted_table.tsv", does the following things:

1- The first metabolite which has the highest value of mean signal intensity is selected.
2- Computes a global correlation matrix. 
3- Select the metabolites which are not correlated to others from the same pcgroup based on a threshold value "Correlation threshold".
4- It creates two tabular files: "correlation_matrix.tsv" (correlation matrix of all the metabolites) and "selected_metabolites_transpo.tsv" (metabolites tagged as deleted are removed and the dataframe is transposed)

The "corr_matrix" function computes a correlation matrix named "correlation_matrix.tsv" and creates a sif file named "sif_table.tsv" (for visualization in CytoScape).v


Galaxy
------
Galaxy is an open, web-based platform for data intensive biomedical research. Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses. 

Homepage: [https://galaxyproject.org/](https://galaxyproject.org/)

Dependencies using Conda
------------------------
[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)  

[Conda](http://conda.pydata.org/) is package manager that among many other things can be used to manage Python packages.


```
#To install miniconda2
#http://conda.pydata.org/miniconda.html
#To install the dependencies libraries using conda:
conda install r-batch r-reshape r-mass
#To set an environment:
conda create -n w4m_correlation r-batch r-reshape r-mass`
#To activate the environment:
. activate w4m_correlation
```

Travis
------
[![Build Status](https://travis-ci.org/workflow4metabolomics/correlation_analysis.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/correlation_analysis)

Test and Deploy with Confidence. Easily sync your GitHub projects with Travis CI and you'll be testing your code in minutes!

Historic contributors
---------------------
 - Antoine Gravot - Protocole conception - Réseau Corsaire - antoine.gravot $ univ-rennes1.fr - France
 - Misharl Monsoor @mmonsoor - for galaxy wrapper and R script - [ABiMS](http://abims.sb-roscoff.fr/) / [IFB](http://www.france-bioinformatique.fr/) - [CNRS](www.cnrs.fr)/[UPMC](www.upmc.fr) - [Station Biologique de Roscoff](http://www.sb-roscoff.fr/) - France
