# freeforall
**************************************************************************************

README

**Replication package for the paper titled "Free-for-all: Does crowding impact outcomes because hospital emergency departments do not prioritise effectively?"**

Author: Igor Francetic

**************************************************************************************

Original version: 08/08/2023

This version: 28/02/2024 (version after R&R used for accepted version of manuscript)

NOTE ON DATA AVAILABILITY:

The data that support the findings of this study are restricted and available from NHS Digital. This repository includes the Stata code used for the analyses and a comprehensive data dictionary. 

The package includes a folder structure wiht 4 do-files. The main do-file to open is named "analysis.do". Running this will, in order:

1 - Call the do-file "preliminary.do", which cleans raw NHS data and prepare datasets ready for analsysis. This is not included and can be provided upon request.

2 - Calls the do-file "avoidable.do", which identifies avoidable/non-avoidable attendances and computes daily (and hourly) volumes of these patients

3 - Run the regression models described in the paper

4 - Call the do-file "fdr_sharpened_qvalues.do" to obtain inference corrected for MHT. Note: the code is adapted directly from Anderson (2008), to which goes all credit. See here https://are.berkeley.edu/~mlanderson/ARE_Website/Research.html.

5 - Call the do-file "outputs.do" to generates the tables and figures used in the paper

For any question or feedback, please email the corresponding author. Contact: igor[dot]francetic[at]manchester[dot]ac[dot]uk
