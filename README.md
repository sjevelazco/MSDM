[![DOI](https://zenodo.org/badge/DOI/10.1016/j.ecolmodel.2020.109180.svg)](https://doi.org/10.1016/j.ecolmodel.2020.109180)
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-lightgrey.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![Lifecycle:
manturing](https://img.shields.io/badge/lifecycle-manturing-blue.svg)](https://www.tidyverse.org/lifecycle/#manturing)
[![R-CMD-check](https://github.com/sjevelazco/MSDM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sjevelazco/MSDM/actions/workflows/R-CMD-check.yaml)


# MSDM
Methods to deal with overprediction of species distribution models


## Installation

```r
require(devtools)  
install_github("sjevelazco/MSDM")  
require(MSDM)

# See more details of funcions and examples
?MSDM_Priori
?MSDM_Posteriori
```

## Description

MSDM provides tools to correct or reduce overprediction of species distribution models. There are two groups of methods compiled in two main functions `MSDM_Priori()` and `MSDM_Posteriori()`. The main difference between both methods is that MSDM_Priori generates predictive spatial variables that will be used (along with the environmental variables) to construct species distribution models, whereas MSDM_Posteriori works on the already fitted models. For more information on MSDM see [Mendes *et al.* (2020)](https://doi.org/10.1016/j.ecolmodel.2020.109180).

*.* **`MSDM_Priori`**: offers four methods, named *XY*, *MIN*,  *CML*, and *KER*. All these methods consist of creating spatial predictor variables to be used with environmental variable to construct species distribution models. These approaches constrain the suitability predicted by species distribution models to species occurrences.

*.* **`MSDM_Posteriori`**: provides four methods, named *OBR*, *PRES*, *LQ*, *MCP*, and *BMCP*. These methods reduce overprediction of species distribution models already fitted based on the occurrences and suitability patterns of species. 

 
 <a href='https://github.com/sjevelazco/MSDM'><img src="https://raw.githubusercontent.com/sjevelazco/MSDM/master/man/Figure/readme_figure.svg" align="centre" height="500"/></a>

  ( a) Example of suitability pattern of an SDM based on a generalized linear model (GLM) constructed without spatial restriction and (b) different A priori and A posteriori MSDM methods </figcaption>

### CITATION:
**Mendes, P., Velazco, S. J. E., Andrade, A. F. A. de, & De Marco, P. (2020). Dealing with overprediction in species distribution models: How adding distance constraints can improve model accuracy. *Ecological Modelling*, 431, 109180.** https://doi.org/10.1016/j.ecolmodel.2020.109180

### NEWS
MSDM functions are now available in the new R package [flexsdm](https://sjevelazco.github.io/flexsdm/) 

## Package website
See [MSDM](https://sjevelazco.github.io/MSDM) package website (https://sjevelazco.github.io/MSDM) for further details of [functions](https://sjevelazco.github.io/MSDM/reference/index.html) and examples 




> Test the package and give us feedback [here](https://github.com/sjevelazco/MSDM/issues) or send an e-mail to sjevelazco@gmail.com.

> MSDM package is integrated to [ENMTML](https://github.com/andrefaa/ENMTML) R package

