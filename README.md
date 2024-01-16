# multiplex-serology-calibration

The code in this repository reproduces the results in the paper:

Emmanuelle A. Dankwa, Martyn Plummer, Daniel Chapman, Rima Jeske, Julia Butt,
Michael Hill, Tim Waterboer, Iona Y. Millwood, Ling Yang, and Christiana Kartsonaki. Calibrating multiplex serology for *Helicobacter pylori*.

# Components and description. 

* `analysis.Rmd`: main analysis file containing R script to reproduce results in manuscript.
* `data`: contains subfolder `new-data` which contains

    +  `training-data.RDS`: standardized training data on a log scale. This is a 498 (rows) x 13 (columns) data frame with rows corresponding to individual observations and columnns corresponding to the immunoblot outcome for each observation (`helicoblot_ckb_overall`) and the 12 predictor antigens: CagA, VacA, GroEL, HP1564, HpaA, HyuA-C, Cad, UreA, NapA, HcPc, Catalase and HP0305. 

    + `training-data.RDS`: standardized training data on a log scale. This is a 498 (rows) x 13 (columns) data frame with rows corresponding to individual observations and columnns corresponding to the immunoblot outcome for each observation (`helicoblot_ckb_overall`) and the 12 predictor antigens.
 
* Displaying surveillance data: `surveillance` folder
* Comparison of relative risks: `risk_analysis` folder
* Hepatitis A virus (HAV) transmission model: `model` folder. 

These have been structured such that each analysis can be run independently of the others. The contents of each folder and guidance on reproducing the results are outlined below. 
