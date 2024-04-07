////////////////////////////////////////////////////////////////////////////////
// Do-file running the analyses for the paper
// Version: 28/02/2024
///////////////////////////////////////////////////////////////////////////////

capture log close
clear all
set more off, perm
set matsize 11000
set maxvar 100000
set processors 12

/* Set paths */
if c(username)=="Add your username" {
	global path "Add your working directory"
}

cd "$path"
*do preliminary // Uncomment to run preliminary.do saved in the same folder, which prepares the dataset from raw NHS data files
use finaldata, clear // The data that support the findings of this study are restricted and only available upon request from NHS Digital

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
global bin_out tot4hr admit lwt disch_nof disch_gp refclin refoth 
global after_out reatt7_all dead30_ons 
global cont_out aecost ninv1 ntrt1 totdur initdur durinvtret
egen day=group(attdate)

compress

////////////////////////////////////////////////////////////////////////////////
// MAIN ESTIMATES (Table 3)
// Estimating Linear models with HDFEs and transforming both attendances
// to the SD scale dividing by SD of unexpected volumes
////////////////////////////////////////////////////////////////////////////////

** Linear HDFE model
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_main
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/main/`y'_main", replace
}

// Test differences between avoidable and non-avoidable on volume scale (not SD scale)
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
scalar sd1=`r(sd)' 
qui: sum unexp_nonavoid_na if e(sample)
scalar sd2=`r(sd)'
display `"`y'"', `"`:var label `y''"', `"`:val label `y''"'
estimates use "Output/Main/`y'_main.ster"
eststo `y': test sd1*est_dailyvol_avoid=sd2*est_dailyvol_nonavoid
estadd scalar wald=`r(p)'
}
esttab $bin_out $after_out $cont_out using "output/main/wald.rtf", keep(*cons*) scalar(wald) p(%9.4f) label replace

////////////////////////////////////////////////////////////////////////////////
// Main models estimating discrete outcomes with Poisson models (Table A.5)
////////////////////////////////////////////////////////////////////////////////

** Poisson HDFE model 
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
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_poisson_main", replace
}

////////////////////////////////////////////////////////////////////////////////
// Sharpened Q for Multiple Hypothesis Testing (multiple outcomes, Table A.4)
** Source: Anderson (2008), "Multiple Inference and Gender Differences in the Effects of Early Intervention: 
** A Reevaluation of the Abecedarian, Perry Preschool, and Early Training Projects", JASA, 103(484), 1481-1495
** Code available from: https://are.berkeley.edu/~mlanderson/ARE_Website/Research.html
////////////////////////////////////////////////////////////////////////////////
preserve
clear all
eststo clear
matrix vec=[.]
foreach y in $bin_out $after_out $cont_out {
estimates use "output/main/`y'_main.ster"
ereturn display
matrix R=r(table)
matrix vec=[vec \ R[4,1]]
matrix vec=[vec \ R[4,2]]
}
svmat vec
rename vec1 pval
drop if missing(pval)
do fdr_sharpened_qvalues // Call Anderson (2008) code, modified to run here
export delimited using "output/main/mht.csv", replace
restore
  
////////////////////////////////////////////////////////////////////////////////
// Heterogeneity 1: Use deciles of unexpected volumes instead of main exposure
// Goal: is there a threshold beyond which detrimental effects are observed?
////////////////////////////////////////////////////////////////////////////////

** Linear HDFE model
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
xtile dailyvol_avoid_q = unexp_avoid_na if e(sample), nq(10)
xtile dailyvol_nonavoid_q = unexp_nonavoid_na if e(sample), nq(10)
reghdfe `y' ib5.dailyvol_avoid_q ib5.dailyvol_nonavoid_q ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop dailyvol_avoid_q dailyvol_nonavoid_q
eststo `y'_decil
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/main/`y'_dec", replace
}

// Test differences between avoidable and non-avoidable on volume scale (not SD scale)
eststo clear
foreach y in $bin_out $after_out $cont_out {
display `"`y'"', `"`:var label `y''"', `"`:val label `y''"'
estimates use "output/main/`y'_dec.ster"
forvalues i=1/10 {
eststo `y'_`i': test `i'.dailyvol_avoid_q=`i'.dailyvol_nonavoid_q
estadd scalar wald_`i'=`r(p)'
}
}
esttab totdur_10 initdur_10 durinvtret_10 tot4hr_10 aecost_10 ninv1_10 ntrt1_10 admit_10 lwt_10 disch_nof_10 disch_gp_10 refclin_10 refoth_10 reatt7_all_10 dead30_ons_10 using "output/main/decwald.rtf", keep(*cons*) scalar(wald_1 wald_2 wald_3 wald_4 wald_6 wald_7 wald_8 wald_9 wald_10) p(%9.4f) label replace


////////////////////////////////////////////////////////////////////////////////
// ROBUSTNESS CHECKS
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Exogeneity of exposure 1: Model without covariates
////////////////////////////////////////////////////////////////////////////////

** Linear HDFE model
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_nocontrols
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/main/`y'_nocontrols", replace
}

////////////////////////////////////////////////////////////////////////////////
// Exogeneity of exposure 2: Patient characteristics on LHS
////////////////////////////////////////////////////////////////////////////////

gen nonwhite=0 if eth==0
replace nonwhite=1 if eth>0&eth<5
gen over65=.
replace over65=0 if arrivalage>0&arrivalage<66
replace over65=1 if arrivalage>64&arrivalage!=.
gen anyatt=.
replace anyatt=0 if nattsyr==0
replace anyatt=1 if nattsyr>0&nattsyr!=.
gen  anyemad=.
replace anyemad=0 if nemadyr==0
replace anyemad=1 if nemadyr>0&nemadyr!=.
global exolinear nattsyr nemadyr over65 

** Linear HDFE model
eststo clear
foreach y in $exolinear {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_exoav
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_exoav", replace
reghdfe `y' est_dailyvol_nonavoid if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_exonav
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_exonav", replace
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
eststo `y'_exoboth
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_exoboth", replace
drop est_dailyvol_avoid est_dailyvol_nonavoid
}

////////////////////////////////////////////////////////////////////////////////
// Supply side 1: Alternative "avoidable attendance" definitions
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Tammes definition instead of NHS 
// Transforming attendance volume to the SD scale
////////////////////////////////////////////////////////////////////////////////

bysort prov_id attdate: egen dailyvol_tamav=total(tammes_avoidable)
bysort prov_id attdate: egen dailyvol_tamnon=total(tammes_nonavoidable)
areg dailyvol_tamav if tammes_avoidable==0, absorb(prov_dow_fe)
predict exp_tamav, xbd
predict unexp_tamav, residual
areg dailyvol_tamnon if tammes_avoidable==0, absorb(prov_dow_fe)
predict exp_tamnon, xbd
predict unexp_tamnon, residual
corr unexp_tamav unexp_tamnon

// Linear HDFE model
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_tamav dailyvol_tamnon ${x1} if tammes_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_tamav if e(sample)
gen est_dailyvol_avoid=dailyvol_tamav/`r(sd)' if e(sample)
qui: sum unexp_tamnon if e(sample)
gen est_dailyvol_nonavoid=dailyvol_tamnon/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if tammes_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_tammes
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_tammes", replace
}
    
//////////////////////////////////////////////////////////////////////////////// 
// Ambulance arrival as definition of non-avoidable instead of NHS
// Transforming attendance volume to the SD scale
////////////////////////////////////////////////////////////////////////////////

bysort prov_id attdate: egen dailyvol_nonavoidamb=total(ambarrival)
gen dailyvol_avoidamb=dailyvol_total-dailyvol_nonavoidamb
areg dailyvol_avoidamb if ambarrival==1, absorb(prov_dow_fe)
predict exp_avoidamb, xbd
predict unexp_avoidamb, residual
areg dailyvol_nonavoidamb if ambarrival==1, absorb(prov_dow_fe)
predict exp_nonavoidamb, xbd
predict unexp_nonavoidamb, residual
corr unexp_nonavoidamb unexp_avoidamb

** Linear HDFE model
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoidamb dailyvol_nonavoidamb ${x1} if ambarrival==1, absorb(prov_dow_fe)
qui: sum unexp_avoidamb if e(sample)
gen est_dailyvol_avoidamb=unexp_avoidamb/`r(sd)' if e(sample)
qui: sum unexp_nonavoidamb if e(sample)
gen est_dailyvol_nonavoidamb=unexp_nonavoidamb/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoidamb est_dailyvol_nonavoidamb ${x1} if ambarrival==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoidamb est_dailyvol_nonavoidamb
eststo `y'_ambdef
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_ambdef", replace
}

////////////////////////////////////////////////////////////////////////////////
// Supply side 2: Heterogeneity across hospital sizes
// Explore heterogeneity across hospitals of different size by dividing 
// the sample into 3 groups based on total attendances
////////////////////////////////////////////////////////////////////////////////

xtile size=dailyvol_total, nq(3)
forvalues i=1/3 {

preserve
keep if size==`i'
** Linear HDFE model
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_size
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_size`i'", replace
}
restore	
}
       
////////////////////////////////////////////////////////////////////////////////
// Supply side 3: Heterogeneity between times of day (in-hours/out-of-hours)
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
** Linear HDFE model
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_ooh
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_ooh", replace
}
restore	


*** In-hours ***
preserve
keep if ooh==0
** Linear HDFE model
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_linear
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_ih", replace
}

////////////////////////////////////////////////////////////////////////////////
// Supply side 4: Heterogeneity across quintiles of hospital bed occupation
////////////////////////////////////////////////////////////////////////////////

forvalues i=1/5 {
preserve
keep if bedocc_yday_q==`i'
** Linear HDFE model
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_bed`i'
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_bed`i'", replace
}
restore	
}

////////////////////////////////////////////////////////////////////////////////
// Supply side 5: Role of measurement error during high pressure periods
// (i) Check if misclassified primary ED diagnosis, or having any treatment/invest.
// Is correlated with unexpected demand from avoidable and non-avoidable patients
// (ii) Check if proba of meeting criteria for NHS definition vary across unexp. demand
////////////////////////////////////////////////////////////////////////////////

// (i) Missclassifed and unexpected demand

gen misclass=.
replace misclass=0 if prim_diag!=.
replace misclass=1 if prim_diag==40

gen any_treat=.
replace any_treat=1 if ntrt1>0&ntrt1!=.
replace any_treat=0 if ntrt1==0

gen any_inv=.
replace any_inv=1 if ninv1>0&ninv1!=.
replace any_inv=0 if ninv1==0

global measurement misclass any_treat any_inv

foreach y in $measurement {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_error
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_error", replace
}

// (ii) Criteria for avoidable and unexpected demand
eststo clear
bysort unexp_avoid_na_q: eststo: quietly estpost summarize nhs_avoidable nhstreat nhsinvest nhsdisch nhsarrmode
esttab using "output/robustness/qvolume_avoid.rtf", cells("mean(fmt(3)) sd(fmt(3))") title("Quintile unexp avoidable") replace
eststo clear
bysort unexp_dailyvol_tot: eststo: quietly estpost summarize nhs_avoidable nhstreat nhsinvest nhsdisch nhsarrmode
esttab using "output/robustness/qvolume_total.rtf", cells("mean sd") title("Quintile unexp total") replace

////////////////////////////////////////////////////////////////////////////////
// Supply side 6: Exclude central London where ambulances may decide where to go
////////////////////////////////////////////////////////////////////////////////

// Flag London CCGs
gen london=0
replace london=1 if inlist(ccg_treatment,"07L","07M","07N","07P","07Q","07R","09A","07T")
replace london=1 if inlist(ccg_treatment, "07V","07W","07X","08A","08C","08D","08E","08F")
replace london=1 if inlist(ccg_treatment, "08G","07Y","08H","08J","08K","08L","08R","08M")
replace london=1 if inlist(ccg_treatment, "08N","08P","08Q","08T","08V","08W","08X","08Y")

// Run analysis excluding London
preserve
drop if london==1

** Linear HDFE model
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_london
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_nolondon", replace
}
restore   

////////////////////////////////////////////////////////////////////////////////
// Supply side 7: Analysis exploiting hourly variation
////////////////////////////////////////////////////////////////////////////////

// Generate FE for ED, hour-of-the-week and month
egen how=group(hour dow)
egen hourlyfe=group(prov_id how month) 
egen clusterhour=group(prov_id how)
egen hod=group(attdate hour)

** Linear HDFE model
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' hourvol_avoid hourvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(hourlyfe)
qui: sum hourvol_avoid if e(sample)
scalar m1=r(sd)
gen est_hourvol_avoid=hourvol_avoid/`r(sd)' if e(sample)
qui: sum hourvol_nonavoid if e(sample)
scalar m2=r(sd)
gen est_hourvol_nonavoid=hourvol_nonavoid/`r(sd)' if e(sample)
areg `y' est_hourvol_avoid est_hourvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(hourlyfe) vce(cluster prov_id)
drop est_*
eststo `y'_hourly
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_hourly", replace
display m1
display m2
}

// Test differences between avoidable and non-avoidable on volume scale (not SD scale)
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' hourvol_avoid hourvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(hourlyfe)
qui: sum hourvol_avoid if e(sample)
scalar sd1=`r(sd)' 
qui: sum hourvol_nonavoid if e(sample)
scalar sd2=`r(sd)'
display `"`y'"', `"`:var label `y''"', `"`:val label `y''"'
estimates use "output/robustness/`y'_hourly.ster"
eststo `y': test sd1*est_hourvol_avoid=sd2*est_hourvol_nonavoid
estadd scalar wald=`r(p)'
}
esttab $bin_out $after_out $cont_out using "output/robustness/wald_hour.rtf", keep(*cons*) scalar(wald) p(%9.4f) label replace


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
** Linear HDFE model
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_age`i'
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_agegr`i'", replace
}
restore		
}

////////////////////////////////////////////////////////////////////////////////
// Demand side 2: Results across discharge destinations
// Ruling out that effects on outcomes are influenced by discharge destination
////////////////////////////////////////////////////////////////////////////////            

global bin_outstab tot4hr
    
*** After-attendance outcomes by discharge group ***

*(a) leaving without treatment
preserve
keep if lwt==1
** Linear HDFE model
eststo clear
foreach y in $bin_outstab $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_lwt
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_lwt", replace
}
restore	

*(b) discharged without follow-up
preserve
keep if disch_nof==1
** Linear HDFE model
eststo clear
foreach y in $bin_outstab $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_nof
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_nof", replace
}
restore	
       
*(c) discharged with GP follow-up
preserve
keep if disch_gp==1
** Linear HDFE model
eststo clear
foreach y in $bin_outstab $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_gpfu
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_gpfu", replace
}
restore	

*(d) referred to other provider
preserve
keep if refclin==1|refoth==1
** Linear HDFE model
eststo clear
foreach y in $bin_outstab $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_refprov
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_refprov", replace
}
restore	
       
////////////////////////////////////////////////////////////////////////////////////////
// Demand side 3: Checking coefficient stability to adding more detailed severity controls 
// Goal: Check coefficient stability to inclusion of better risk-adjustment vars
// on the sub-sample of patients admitted to inpatient care for which we have
// more accurate ICD-10 codes
////////////////////////////////////////////////////////////////////////////////////////
       
////////////////////////////////////////////////////////////////////////////////
// Merge in primary diagnosis for admitted patients
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
*/

preserve
clear all
use admitted
compress

* Coefficient stability analysis on patients admitted to inpatient care *
global bin_outstab tot4hr

** WITHOUT ELIXHAUSERS COMORBIDITIES **

** Linear HDFE model
eststo clear
foreach y in $bin_outalt $after_out  {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x2} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x2} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_stab1
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_stab1", replace
}

** WITH ELIXHAUSERS COMORBIDITIES **

** Linear HDFE model
eststo clear
foreach y in $cont_out $bin_outstab $after_out  {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_stab2
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_stab2", replace
}
     
** WITH ELIXHAUSERS COMORBIDITIES AND INPATIENT ICD10 CODE**

** Linear HDFE model
eststo clear
foreach y in $cont_out $bin_outstab $after_out  {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} i.recode_icd10 if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe recode_icd10) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_stab3
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_stab3", replace
}
restore

////////////////////////////////////////////////////////////////////////////////
// Specification check: Run main models including ineraction term
////////////////////////////////////////////////////////////////////////////////

** Linear HDFE model
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
gen vol_inter=est_dailyvol_nonavoid*est_dailyvol_avoid
reghdfe `y' est_dailyvol_avoid est_dailyvol_nonavoid vol_inter ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid vol_inter
eststo `y'_interacted
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_interacted", replace
}

////////////////////////////////////////////////////////////////////////////////////////
// CHECK STABILITY IN YEAR 2017/2018
////////////////////////////////////////////////////////////////////////////////////////
clear all
*do prepare_altyear // Uncomment to run prepare_altyear.do saved in the same folder, which prepares the dataset from raw NHS data files
use "finaldata_altyear", clear
keep if nhs_nonavoidable==1
compress 

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

// Outcomes 
global cont_out ninv1 ntrt1 totdur initdur durinvtret
global bin_out tot4hr admit lwt disch_nof disch_gp refclin refoth 
global after_out reatt7_all dead30_ons 

// Estimating Linear/Poisson models with HDFEs and transforming both attendances
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
gen est_dailyvol_avoid=dailyvol_avoid/`r(sd)' if e(sample)
qui: sum unexp_nonavoid_na if e(sample)
gen est_dailyvol_nonavoid=dailyvol_nonavoid/`r(sd)' if e(sample)
areg `y' est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
drop est_dailyvol_avoid est_dailyvol_nonavoid
eststo `y'_altyear
sum `y' if e(sample)
estadd scalar mean = `r(mean)'
estimates save "output/robustness/`y'_altyear", replace
}

// Test differences between avoidable and non-avoidable on volume scale (not SD scale)
eststo clear
foreach y in $bin_out $after_out $cont_out {
qui: areg `y' dailyvol_avoid dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe)
qui: sum unexp_avoid_na if e(sample)
scalar sd1=`r(sd)' 
qui: sum unexp_nonavoid_na if e(sample)
scalar sd2=`r(sd)'
display `"`y'"', `"`:var label `y''"', `"`:val label `y''"'
estimates use "output/robustness/`y'_altyear.ster"
eststo `y': test sd1*est_dailyvol_avoid=sd2*est_dailyvol_nonavoid
estadd scalar wald=`r(p)'
}
esttab $bin_out $after_out $cont_out using "output/robustness/wald_altyear.rtf", keep(*cons*) scalar(wald) p(%9.4f) label replace
