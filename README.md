# Overview

The R code in this repository reproduces the results in the paper:

Emmanuelle A. Dankwa, Martyn Plummer, Daniel Chapman, Rima Jeske, Julia Butt,
Michael Hill, Tim Waterboer, Iona Y. Millwood, Ling Yang, and Christiana Kartsonaki. Calibrating multiplex serology for *Helicobacter pylori*.

Analyses were run in R version 4.0.5. Package versions and references are provided in manuscript.

# Components and description 

* `analysis.Rmd`: main analysis file containing R script to reproduce results in manuscript.
  
* `data`: contains subfolder `new-data` which contains

    +  `training-data.RDS`: standardized training data on a log scale. This is a 498 (rows) x 13 (columns) data frame with rows corresponding to individual observations and columnns corresponding to the immunoblot outcome for each observation (`helicoblot_ckb_overall`) and the 12 predictor antigens: CagA, VacA, GroEL, Omp (renamed to HP1564 for specificity), HpaA, HyuA-C, Cad, UreA, NapA, HcPc, Catalase and HP0305. 

    + `test-data.RDS`: standardized test data on a log scale. This is a 922 (rows) x 13 (columns) data frame with rows corresponding to individual observations and columnns corresponding to the immunoblot outcome for each observation (`helicoblot_ckb_overall`) and the 12 predictor antigens.
 
* `results`: contains the outputs of `analysis.Rmd`. Figures are separately stored in the subfolder `figs`.
  
* `monbart_0.3`: R package to run the monotone BART algorithm, based on the C++ implementation by Chipman et al. (2022), available at [https://bitbucket.org/remcc/mbart/src/master/](https://bitbucket.org/remcc/mbart/src/master/). We created a 
custom R interface to the C++ code for the purposes of this study.


# To reproduce results in paper

* Install the monbart_0.3 package. You may use one of the following ways.
   +  Packages > Install from > Package Archive File (.tgz; .tar.gz) > monbart_0.3/monbart_0.3.tar.gz. 
   +  Run `R CMD INSTALL monbart_0.3/monbart_0.3.tar.gz` at the terminal. 
* Run the script `analysis.Rmd`.

# Questions 

Please contact Emmanuelle at [edankwa@hsph.harvard.edu](edankwa@hsph.harvard.edu). 

# Reference

Hugh A. Chipman. Edward I. George. Robert E. McCulloch. Thomas S. Shively. "mBART: Multidimensional Monotone BART." Bayesian Anal. 17 (2) 515 - 544, June 2022.


