# freeforall
**************************************************************************************

README

**Replication package for the paper titled "Free-for-all: Does crowding impact outcomes because hospital emergency departments do not prioritise effectively?"**

Author: Igor Francetic

**************************************************************************************

Original version: 08/08/2023

This version: 28/02/2024 (R&R version)

NOTE ON DATA AVAILABILITY:

The data that support the findings of this study are restricted and available from NHS Digital. This repository includes the Stata code used for the analyses and a comprehensive data dictionary. Note: do-files that clean raw NHS data and prepare datasets ready for analsysis can be provided upon request.

The package includes a folder structure and two do-files, that should be run in the following order:

1 - analysis.do: Runs the regression models described in the paper

2 - outputs.do:  Generates the tables and figures used in the paper

3 - fdr_sharpened_qvalues.do: Generate inference corrected for MHT, the code is adapted directly from Anderson (2008): https://are.berkeley.edu/~mlanderson/ARE_Website/Research.html 

For any question or feedback, please email the corresponding author. Contact: igor[dot]francetic[at]manchester[dot]ac[dot]uk
