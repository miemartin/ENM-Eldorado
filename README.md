# ENM-Eldorado
This repository contains the scripts that we have developed to build Ecological Niche Models and potential geographical distribution of Aedes aegypti and Aedes albopictus in Eldorado city, Misiones, Argentina.

All pre- and post-processing of remote sensing data was previously carried out in QGIS 3.18.3 software (https://www.qgis.org/, accessed in 2023).

The order in which the scripts should be used is as follows:

1. Occurrence data: This script imports the occurrence data and does all the processing for it. It eliminates duplicate data, assigns only one data per pixel and finally creates a .csv file where the values of the environmental variables are obtained for each location.

2. VIF analysis: this is the script where the selection of environmental variables is carried out based on the VIF values and Pearson correlation.

3. Modeling candidate Ae. aegypti: this is the script to perform the calibration of candidate models from the combination of different feature classes and regularization multipliers.

4. Best model Ae. aegypti: this is the script that runs the best selected model. Applies 25-fold cross validation. The model metrics, the average predictive map, and the response curves of the environmental variables used are obtained. Then the model is transferred to the new study area (in this case, Northeast Argentina) and finally the binary map is obtained by using the “p10” threshold rule.

NOTE: this entire procedure was repeated for each season of the year and for Ae. albopictus.

5. Average maps: This is the script, the average predictive maps and binary maps are obtained by season of the year (average of, for example, spring-summer 2016 and spring-summer 2017), and the area (in percentages) suitable for the species is calculated from the previously obtained binary maps.
