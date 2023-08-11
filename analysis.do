////////////////////////////////////////////////////////////////////////////////
// Do-file running the analyses for the paper
// Version: 08/08/2023
///////////////////////////////////////////////////////////////////////////////

capture log close
clear all
set more off, perm
set matsize 11000
set maxvar 100000

/* Set paths */
if c(username)=="Add your username" {
	global path "Add your working directory"
}

cd "$path"
*do preliminary // Uncomment to run preliminary.do saved in the same folder
*use finaldata, clear // The data that support the findings of this study are restricted and only available upon request from NHS Digital

////////////////////////////////////////////////////////////////////////////////
*********************************** ANALYSIS ***********************************
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// CREATE GLOBALS
////////////////////////////////////////////////////////////////////////////////

// Covariates

// Pooling elixhauser comorbidities together
global elixyear nchf2yr narrhy2yr nvd2yr npcd2yr npvd2yr nhptn_uc2yr nhptn_c2yr npara2yr nothnd2yr ncopd2yr ndiab_uc2yr ndiab_c2yr nhypothy2yr nrf2yr nld2yr npud_nb2yr nhiv2yr nlymp2yr nmets2yr ntumor2yr nrheum_a2yr ncoag2yr nobesity2yr nwl2yr nfluid2yr nbla2yr nda2yr nalcohol2yr ndrug2yr npsycho2yr ndep2yr

// Full set of covariates
global x1 i.agesex ib0.eth i.imd2015_cat ib41.prim_diag ib0.nattsyr ib0.nemadyr ${elixyear} ib0.nelixyr ib2.arrmode i.attcat ib80.patgrp ib1.refsource ib10.incloc ib1.bedocc_yday_q

// Subset of covariates to be used in the coefficient stability test
global x2 i.agesex ib0.eth i.imd2015_cat ib41.prim_diag ib0.nattsyr ib0.nemadyr ib2.arrmode i.attcat ib80.patgrp ib1.refsource ib10.incloc ib1.bedocc_yday_q

// Outcomes 
global cont_out aecost ninv1 ntrt1 totdur initdur durinvtret
global bin_out tot4hr admit lwt disch_nof disch_gp refclin refoth
global bin_outalt tot4hr // Global exlcuding discharge destinations from binary outcomes for stratified analyses conditioning on discharge
global after_out reatt7_all reatt30_all dead7_ons dead30_ons 

////////////////////////////////////////////////////////////////////////////////
//  MAIN ESTIMATES
////////////////////////////////////////////////////////////////////////////////
// Estimating Linear/Poisson models with HDFEs scaling attendance volumes by  SD

** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_out $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/main/`y'_allhdfe", replace
}
** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/main/`y'_poisson_allhdfe", replace
}
     
////////////////////////////////////////////////////////////////////////////////
// Heterogeneity 1: Use deciles of unexpected volumes instead of main exposure
// Goal: is there a threshold beyond which detrimental effects are observed?
////////////////////////////////////////////////////////////////////////////////

** Linear HDFE model for binary outcomes 
eststo clear
foreach y in $bin_out $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
xtile dailyvol_avoid_q = unexp_avoid_na if e(sample), nq(10)
xtile dailyvol_nonavoid_q = unexp_nonavoid_na if e(sample), nq(10)
reghdfe `y' ib5.dailyvol_avoid_q ib5.dailyvol_nonavoid_q ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop dailyvol_avoid_q dailyvol_nonavoid_q
eststo `y'_decil_nonavoid
estimates save "output/main/`y'_allhdfe_dec", replace
}

** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
xtile dailyvol_avoid_q = unexp_avoid_na if e(sample), nq(10)
xtile dailyvol_nonavoid_q = unexp_nonavoid_na if e(sample), nq(10)
ppmlhdfe `y' ib5.dailyvol_avoid_q ib5.dailyvol_nonavoid_q  ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins rb5.dailyvol_avoid_q rb5.dailyvol_nonavoid_q, post
drop dailyvol_avoid_q dailyvol_nonavoid_q 
estimates save "output/main/`y'_poisson_dec", replace
}

////////////////////////////////////////////////////////////////////////////////
// Heterogeneity 2: Estimates model across tertiles of NAV/AV ratio
// Goal: Are the results consistent with the model in relation to higher NA %?
////////////////////////////////////////////////////////////////////////////////
gen rationavav=dailyvol_nonavoid/dailyvol_total
xtile tert_ratio=rationavav, nq(3)

forvalues i=1/3 {
preserve
keep if tert_ratio==`i'
// Estimating Linear/Poisson models with HDFEs and transforming both attendances
// to the SD scale dividing by SD of unexpected volumes
** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_out $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/main/`y'_ratio`i'", replace
}
** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/main/`y'_ratio`i'", replace
}

restore
}
  
////////////////////////////////////////////////////////////////////////////////
// ROBUSTNESS CHECKS
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Supply side 1: Reverse causality in defining avoidable/non-avoidable patients
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Determinants of NHS avoidable state across distribution of unexpected volumes
////////////////////////////////////////////////////////////////////////////////
eststo clear
bysort unexp_avoid_na_q: eststo: quietly estpost summarize nhs_avoidable nhstreat nhsinvest nhsdisch nhsarrmode reatt7_all reatt30_all reatt7_unplan reatt30_unplan
esttab using "output/robustness/qvolume_avoid.rtf", cells("mean(fmt(3)) sd(fmt(3))") title("Quintile unexp avoidable") replace
eststo clear
bysort unexp_avoid_na_q: eststo: quietly estpost summarize nhs_avoidable nhstreat nhsinvest nhsdisch nhsarrmode reatt7_all reatt30_all reatt7_unplan reatt30_unplan
esttab using "output/robustness/qvolume_total.rtf", cells("mean sd") title("Quintile unexp total") replace
eststo clear
bysort unexp_avoid_na_q: eststo: quietly estpost summarize $count $bin $after
esttab using "output/robustness/qvolume_outcomes.rtf", cells("mean sd") title("Outcomes across quintiles of avoidable att (full sample)") replace  
   
////////////////////////////////////////////////////////////////////////////////
// Tammes definition instead of NHS 
// Transforming attendance volume to the SD scale
////////////////////////////////////////////////////////////////////////////////

bysort prov_id arrivaldate: egen dailyvol_tamav=total(tammes_avoidable)
bysort prov_id arrivaldate: egen dailyvol_tamnon=total(tammes_nonavoidable)
areg dailyvol_tamav if tammes_avoidable==0, absorb(prov_dow_fe)
predict exp_tamav, xbd
predict unexp_tamav, residual
areg dailyvol_tamnon if tammes_avoidable==0, absorb(prov_dow_fe)
predict exp_tamnon, xbd
predict unexp_tamnon, residual

** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_out $after_out {
qui: areg `y' dailyvol_tamav dailyvol_tamnon ${x1} if tammes_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_tamav if e(sample)
gen est_dailyvol_avoid=dailyvol_tamav/`r(sd)' if e(sample)
qui: sum unexp_tamnon if e(sample)
gen est_dailyvol_nonavoid=dailyvol_tamnon/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if tammes_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_tammes
estimates save "output/robustness/`y'_tammes", replace
}

** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_tamav dailyvol_tamnon ${x1} if tammes_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_tamav if e(sample)
gen est_dailyvol_avoid=dailyvol_tamav/`r(sd)' if e(sample)
qui: sum unexp_tamnon if e(sample)
gen est_dailyvol_nonavoid=dailyvol_tamnon/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if tammes_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_tammes
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_poisson_tammes", replace
}
       
//////////////////////////////////////////////////////////////////////////////// 
// Ambulance arrival as definition of non-avoidable instead of NHS
// Transforming attendance volume to the SD scale
////////////////////////////////////////////////////////////////////////////////

bysort prov_id arrivaldate: egen dailyvol_nonavoidamb=total(ambarrival)
gen dailyvol_avoidamb=dailyvol_total-dailyvol_nonavoidamb
areg dailyvol_avoidamb if ambarrival==1, absorb(prov_dow_fe)
predict exp_avoidamb, xbd
predict unexp_nonavoidamb, residual
areg dailyvol_nonavoidamb if ambarrival==1, absorb(prov_dow_fe)
predict exp_avoidamb, xbd
predict unexp_avoidamb, residual

** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_out $after_out {
qui: areg `y' dailyvol_avoidamb dailyvol_nonavoidamb ${x1} if ambarrival==1, absorb(prov_dow_fe)
qui: sum unexp_avoidamb if e(sample)
gen est_dailyvol_avoidamb=unexp_avoidamb/`r(sd)' if e(sample)
qui: sum unexp_nonavoidamb if e(sample)
gen est_dailyvol_nonavoidamb=unexp_nonavoidamb/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoidamb est_dailyvol_nonavoidamb ${x1} if ambarrival==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoidamb est_dailyvol_nonavoidamb
eststo `y'_ambdef
estimates save "output/robustness/`y'_ambdef", replace
}

** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoidamb dailyvol_nonavoidamb ${x1} if ambarrival==1, absorb(prov_dow_fe)
qui: sum dailyvol_avoidamb if e(sample)
gen est_dailyvol_avoidamb=dailyvol_avoidamb/`r(sd)' if e(sample)
qui: sum dailyvol_nonavoidamb if e(sample)
gen est_dailyvol_nonavoidamb=dailyvol_nonavoidamb/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoidamb est_dailyvol_nonavoidamb ${x1} if ambarrival==1, d absorb(prov_dow_fe) vce(cluster prov_id)
margins, dydx(est_dailyvol_avoidamb est_dailyvol_nonavoidamb) post
drop est_dailyvol_avoidamb est_dailyvol_nonavoidamb
estimates save "output/robustness/`y'_poisson_ambdef", replace
}

////////////////////////////////////////////////////////////////////////////////
// Supply side 2: Heterogeneity across hospital sizes
// Explore heterogeneity across hospitals of different size by dividing the
// the sample into 3 groups based on total volume of ED attendances
////////////////////////////////////////////////////////////////////////////////

xtile size=dailyvol_total, nq(3)
forvalues i=1/3 {

preserve
keep if size==`i'
** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_out $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/robustness/`y'_size`i'", replace
}

** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_size`i'", replace
}

restore	
	
}
       
////////////////////////////////////////////////////////////////////////////////
// Supply side 3: Heterogeneity between times of day (in-hours/out-of-hours)
// May reveal different input resource constraints (more stringent out-of-hours)
////////////////////////////////////////////////////////////////////////////////       

// Define out of hours
recode hour (0/5=1) (18/23=1) (6/17=0), gen(night)
recode dow (1/5=0) (6/7=1), gen(weekend)
gen ooh=0
replace ooh=1 if (weekend==1)|(weekend==0&night==1)
replace ooh=. if weekend==.|night==.

*** Out-of-hours ***
preserve
keep if ooh==1
** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_out $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/robustness/`y'_ooh", replace
}

** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_ooh", replace
}
restore	


*** In-hours ***
preserve
keep if ooh==0
** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_out $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/robustness/`y'_ih", replace
}

** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_ih", replace
}
restore	

////////////////////////////////////////////////////////////////////////////////
// Supply side 4: Heterogeneity across quintiles of hospital bed occupation
// May reveal broader hospital capacity constraints limiting proba of admission
////////////////////////////////////////////////////////////////////////////////

forvalues i=1/5 {

preserve
keep if bedocc_yday_q==`i'

** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_out $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/robustness/`y'_bed`i'", replace
}

** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_bed`i'", replace
}

restore	
	
}

////////////////////////////////////////////////////////////////////////////////
// Supply side 5: Exclude central London where ambulances may decide where to go
////////////////////////////////////////////////////////////////////////////////

// Flag London CCGs
gen london=0
replace london=1 if inlist(ccg_treatment,"07L","07M","07N","07P","07Q","07R","09A","07T")
replace london=1 if inlist(ccg_treatment, "07V","07W","07X","08A","08C","08D","08E","08F")
replace london=1 if inlist(ccg_treatment, "08G","07Y","08H","08J","08K","08L","08R","08M")
replace london=1 if inlist(ccg_treatment, "08N","08P","08Q","08T","08V","08W","08X","08Y")

// Estimating Linear/Poisson models with HDFEs and scaling attendances volumes
// by the SD of unexpected volumes, and excluding EDs in London
preserve
drop if london==1
** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_out $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/robustness/`y'_nolondon", replace
}
** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_nolondon", replace
}
restore   

////////////////////////////////////////////////////////////////////////////////
// Demand side 1: Heterogeneity across age groups
// Reveal age-related prioritisation patterns (e.g. young children)
////////////////////////////////////////////////////////////////////////////////

gen agegr=1 if arrivalage>0&arrivalage<10
replace agegr=2 if arrivalage>10&arrivalage<20
replace agegr=3 if arrivalage>19&arrivalage<30
replace agegr=4 if arrivalage>29&arrivalage<60
replace agegr=5 if arrivalage>59&arrivalage<75
replace agegr=6 if arrivalage>74&arrivalage!=.

forvalues i=1/6 {

preserve
keep if agegr==`i'
** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_out $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/robustness/`y'_agegr`i'", replace
}

** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_agegr`i'", replace
}

restore	
	
}

////////////////////////////////////////////////////////////////////////////////
// Demand side 2: Results across discharge destinations
// Ruling out that effects on outcomes are influenced by discharge destination
////////////////////////////////////////////////////////////////////////////////            
    
*(a) leaving without treatment
preserve
keep if lwt==1
** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_outalt $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/robustness/`y'_lwt", replace
}
** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_lwt", replace
}
restore	

*(b) discharged without follow-up
preserve
keep if disch_nof==1
** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_outalt $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/robustness/`y'_nof", replace
}

** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_nof", replace
}
restore	
       
*(c) discharged with GP follow-up
preserve
keep if disch_gp==1
** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_outalt $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/robustness/`y'_gpfu", replace
}

** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_gpfu", replace
}
restore	

*(d) referred to other provider
preserve
keep if refclin==1|refoth==1
** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_outalt $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/robustness/`y'_refprov", replace
}

** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_refprov", replace
}
restore	
       
////////////////////////////////////////////////////////////////////////////////
// Demand side 3: Coefficient stability to adding detailed severity controls 
// Goal: Check coefficient stability to inclusion of better risk-adjustment vars
// on the sub-sample of patients admitted to inpatient care for which we have
// more accurate ICD-10 codes
////////////////////////////////////////////////////////////////////////////////
       
////////////////////////////////////////////////////////////////////////////////
// Set up primary diagnosis data for admitted patients
////////////////////////////////////////////////////////////////////////////////

preserve
clear all
use aekey admidate diag_01 using "$rawdatapath/apc2017.dta" // Opens inpatient data for 2017
keep if aekey!=.
rename diag_01 inpat_primdiag
rename admidate arrivaldate
gen date=date(arrivaldate, "YMD")
replace date=date-1
keep if arrivaldate=="2017-04-01"
save inpatdiag17, replace

clear all
use aekey admidate diag_01 admisorc using "$rawdatapath/apc2016.dta" // Opens inpatient data for 2016
keep if aekey!=.
rename diag_01 inpat_primdiag
rename admidate arrivaldate
gen date=date(arrivaldate, "YMD")
replace date=date-1
save inpatdiag16, replace
append using inpatdiag17
gen icd10=substr(inpat_primdiag,1,3)
save inpatdiag, replace
restore

preserve
keep if admit==1
gen date=date(arrivaldate, "YMD")
merge 1:1 aekey arrivaldate using inpatdiag, keepusing(icd10)
drop if _merge==2
drop _merge
merge 1:1 aekey date using inpatdiag, keepusing(icd10)
drop if _merge==2
drop _merge
replace icd10="Unknown" if icd10==""
egen recode_icd10=group(icd10)
save admitted, replace
restore

preserve
clear all
use admitted, clear

** WITHOUT ELIXHAUSERS COMORBIDITIES

** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_outalt $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x2} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x2} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/robustness/`y'_stab1", replace
}
** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x2} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x2} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_stab1", replace
}

** WITH ELIXHAUSERS COMORBIDITIES

** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_outalt $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/robustness/`y'_stab2", replace
}
** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_stab2", replace
}
     
** WITH ELIXHAUSERS COMORBIDITIES AND INPATIENT ICD10 CODE

** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_outalt $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} i.recode_icd10 if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe recode_icd10) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
estimates save "output/robustness/`y'_stab3", replace
}
** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} i.recode_icd10 if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe recode_icd10) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid) post
drop est_dailyvol_avoid est_dailyvol_nonavoid
estimates save "output/robustness/`y'_stab3", replace
}
restore

////////////////////////////////////////////////////////////////////////////////
// Specification check: Run main models including ineraction term
////////////////////////////////////////////////////////////////////////////////

// Estimating Linear/Poisson models with HDFEs scaling attendance volumes by  SD

** Linear HDFE model for binary outcomes
eststo clear
foreach y in $bin_out $after_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
gen vol_inter=est_dailyvol_nonavoid*est_dailyvol_avoid
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid vol_inter ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid vol_inter
eststo `y'_linear
estimates save "output/interacted/`y'_allhdfe", replace
}
** Poisson HDFE model for continous outcomes
eststo clear
foreach y in $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
gen vol_inter=est_dailyvol_nonavoid*est_dailyvol_avoid
ppmlhdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid vol_inter ${x1} if nhs_nonavoidable==1, d absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_poisson_nonavoid
margins, dydx(est_dailyvol_avoid est_dailyvol_nonavoid vol_inter) post
drop est_dailyvol_avoid est_dailyvol_nonavoid vol_inter
estimates save "output/interacted/`y'_poisson_allhdfe", replace
}