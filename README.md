# MSDM
Overprediction correction approaches in species distribution models


## Install MSDM package

```r
require(devtools)  
install_github("sjevelazco/MSDM")  
require(MSDM)
```


## Description

MSDM provides tools to correct overprediction of species distribution models. There is two main functions `MSDM_Priori()` and `MSDM_Posteriori()`. 

*.* **`MSDM_Priori`**: offers four methods, named *XY*, *MIN*,  *CML*, and *KER*. These creates spatial predictor variables which have to be used with environmental variable to construct species distribution models. These approaches constrain the suitability predicted by SDMs to species occurrences.

*.* **`MSDM_Posteriori`**: provides four methods, named *OBR*, *PRES*, *LQ*, *MCP*, and *BMCP*. All theme correct overprediction of species distribution models based on occurrences and suitability patterns of species. 


