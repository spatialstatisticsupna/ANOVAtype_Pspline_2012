## **R CODE**

This folder contains the R code to fit the P-spline ANOVA-type interaction model described in the paper *A P-spline ANOVA type model in space-time disease mapping. Stochastic Environmental Research and Risk Assessment, 26(6), 835-845.* To fit this model, the user should open the **MainFile.R** file

In the **steps 1-to-6** of this file, the data set **ProstH.txt** is loaded, the B-splines bases and penalty matrices are created and afterwards, the matrices of fixed and random effects of the mixed model are generated. These elements are automatically saved to the Dumpsdata folder.

In the **step 7**, the ANOVA type model is fitted using the PQL estimation technique. In this step **PQLs.R** is used. This will be automatically done running `source("PQLs.r")` code. This step could consume some computational time.

In the **step 8,** the mean squared error estimation is provided. To do so, the code `source("msespline.r")` should be run. Finally, the values of AIC and BIC are computed. The process ends saving the results as **ANOVA_Results.RDATA**.

More about the **ProstH.txt** data set. This file includes information on the following variables:

-   **region**: Takes the values 1 to 50 representing different regions of Spain

-   **Year**: Takes values 1975 to 2008 representing a calendar year

-   **cases**: Observed number of deaths

-   **pop**: Population at risk in each region/year.

-   **esp**: Espected number of deaths

-   **long**: Longitud of the centroid of each region

-   **lat**: Latitude of the centroid of each region

## R and R packages. Version Info

``` r
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 17763)

Matrix products: default

Random number generation:
 RNG:     Mersenne-Twister 
 Normal:  Inversion 
 Sample:  Rounding 
 
locale:
[1] LC_COLLATE=Spanish_Spain.1252  LC_CTYPE=Spanish_Spain.1252   
[3] LC_MONETARY=Spanish_Spain.1252 LC_NUMERIC=C                  
[5] LC_TIME=Spanish_Spain.1252    

attached base packages:
[1] splines   stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] MASS_7.3-57

loaded via a namespace (and not attached):
 [1] compiler_4.2.1  fastmap_1.1.0   cli_3.3.0       htmltools_0.5.3 tools_4.2.1    
 [6] rstudioapi_0.14 yaml_2.3.5      rmarkdown_2.16  knitr_1.40      xfun_0.32      
[11] digest_0.6.29   rlang_1.0.5     evaluate_0.16 
```
