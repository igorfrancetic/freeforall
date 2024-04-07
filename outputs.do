////////////////////////////////////////////////////////////////////////////////
// Do-file producing outputs for the paper
// Version: 28/02/2024
///////////////////////////////////////////////////////////////////////////////

/* Set paths */
if c(username)=="Add your username" {
	global path "Add your working directory"
}
cd "$path"
set scheme s1mono

// Open main dataset
use finaldata, clear

// Set global for outcomes
global elixyear nchf2yr narrhy2yr nvd2yr npcd2yr npvd2yr nhptn_uc2yr nhptn_c2yr npara2yr nothnd2yr ncopd2yr ndiab_uc2yr ndiab_c2yr nhypothy2yr nrf2yr nld2yr npud_nb2yr nhiv2yr nlymp2yr nmets2yr ntumor2yr nrheum_a2yr ncoag2yr nobesity2yr nwl2yr nfluid2yr nbla2yr nda2yr nalcohol2yr ndrug2yr npsycho2yr ndep2yr
global x1 i.agesex ib0.eth i.imd2015_cat ib41.prim_diag ib0.nattsyr ib0.nemadyr ${elixyear} ib0.nelixyr ib2.arrmode i.attcat ib80.patgrp ib1.refsource ib10.incloc ib1.bedocc_yday_q 

global cont_out aecost ninv1 ntrt1 totdur initdur durinvtret
global bin_out tot4hr admit lwt disch_nof disch_gp refclin refoth 
global after_out reatt7_all dead30_ons 
global tret_out aecost ninv1 ntrt1 
global wt_out totdur durinvtret initdur
global binno4_out refoth 
global restore est_dailyvol_avoid est_dailyvol_nonavoid agesex eth imd2015_cat prim_diag nattsyr nemadyr ${elixyear} nelixyr arrmode attcat patgrp refsource incloc bedocc_yday_q prov_dow_fe

global missingvars agesex eth imd2015_cat prim_diag nattsyr nemadyr ${elixyear} nelixyr arrmode attcat patgrp refsource incloc bedocc_yday_q

// Generate variables needed to generate outputs
gen est_dailyvol_avoid=dailyvol_avoid 
gen est_dailyvol_nonavoid=dailyvol_nonavoid 
gen vol_inter=.

// Relabel a few variables
label var est_dailyvol_avoid "Unexp. volume of avoid." 
label var est_dailyvol_nonavoid "Unexp. volume of non-avoid."
label var vol_inter "Interaction AV/NAV"

////////////////////////////////////////////////////////////////////////////////
// Main models, produce RTF tables (TeX also possible)
////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $bin_out {
estimates use "output/main/`y'_main.ster"
eststo `y'_linear
}
esttab using "output/main/bin.rtf", /// 
keep(est_dailyvol_avoid est_dailyvol_nonavoid) replace ci(%9.4f) b(%9.4f) label scalar(N r2_a mean) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons)
       
eststo clear
foreach y in $after_out  {
estimates use "output/main/`y'_main.ster"
eststo `y'_linear
}
esttab using "output/main/after.rtf", ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) replace ci(%9.4f) b(%9.4f) label scalar(N r2_a mean) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons)     

eststo clear
foreach y in $cont_out {
estimates use "output/main/`y'_main.ster"
eststo `y'_main
}
esttab using "output/main/cont.rtf", ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) replace ci(%9.4f) b(%9.4f) label scalar(N r2_a mean) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons)      

eststo clear
foreach y in $cont_out {
estimates use "output/robustness/`y'_poisson_main.ster"
eststo `y'_poisson
}
esttab using "output/robustness/cont_poisson.rtf", ///
replace ci(%9.4f) b(%9.4f) label scalar(N r2_a mean) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons) 

////////////////////////////////////////////////////////////////////////////////
// Outputs for Heterogeneity checks
// Models estimated on deciles of unexpected attendances, produce coefficient plots
////////////////////////////////////////////////////////////////////////////////

*Coefplots avoidable and nonavoidable
eststo clear
foreach y in $bin_out $after_out $cont_out {
local vla : variable label `y'
estimates use "output/main/`y'_dec.ster"
eststo `y'
coefplot (`y', keep(*.dailyvol_avoid_q) drop(5.dailyvol_avoid_q) label("Unexp. Avoidable ") pstyle(p3)) ///
	 (`y', keep(*.dailyvol_nonavoid_q) drop(5.dailyvol_nonavoid_q) label("Unexp. Non-avoidable") pstyle(p4)), ///
	 yline(0, lcolor(black) lwidth(thin) lpattern(dash)) vertical ///
         mlabel format(%6.5f) mlabsize(vsmall) mlabposition(12) mlabgap(*8) ylab(, angle(45) labsize(small)) ///
	 xtitle("Att. volume shock decile (vs. Median)") ///
	 ytitle("Marginal effect on outcome") title("`vla'") ///
	 ciopts(lwidth(*2)) levels(99) ///
	 coeflabels(1.dailyvol_avoid_q = "1" 2.dailyvol_avoid_q = "2" 3.dailyvol_avoid_q = "3" 4.dailyvol_avoid_q = "4" ///
	           6.dailyvol_avoid_q = "6" 7.dailyvol_avoid_q = "7" 8.dailyvol_avoid_q = "8" 9.dailyvol_avoid_q = "9" ///
	           10.dailyvol_avoid_q = "10" 1.dailyvol_nonavoid_q = "1" 2.dailyvol_nonavoid_q = "2" 3.dailyvol_nonavoid_q = "3" 4.dailyvol_nonavoid_q = "4" ///
	           6.dailyvol_nonavoid_q = "6" 7.dailyvol_nonavoid_q = "7" 8.dailyvol_nonavoid_q = "8" 9.dailyvol_nonavoid_q = "9" ///
	           10.dailyvol_nonavoid_q = "10")
graph save "figures/`y'_deccomb.gph", replace
graph export "figures/`y'_deccomb.png", replace
}

////////////////////////////////////////////////////////////////////////////////
// ROBUSTNESS CHECKS
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Models with no controls
////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $cont_out $bin_out $after_out {
eststo: estimates use "output/main/`y'_nocontrols.ster"
}
esttab using "output/robustness/nocontrols.rtf",
replace ci(%9.4f) b(%9.4f) label scalar(N r2_a) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons) keep(est_*) 

eststo clear
foreach y in $cont_out $bin_out $after_out {
local vla : variable label `y'	
estimates use "output/main/`y'_main"
eststo `y'_main
estimates use "output/main/`y'_nocontrols.ster"
eststo `y'_nocontrol
coefplot (`y'_main, label("Main") msymbol(D) mlcolor(magenta)) (`y'_nocontrol, label("No controls")), ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) ///
saving(figures/`y'_nocontrols, replace) level(99) legend(cols(3) span stack)
graph export "figures/`y'_nocontrols.png", replace
}

eststo clear
foreach y in $cont_out $bin_out $after_out {
local vla : variable label `y'	
estimates use "output/main/`y'_main"
eststo `y'_main
estimates use "output/main/`y'_nocontrols.ster"
eststo `y'_nocontrol
coefplot (`y'_main, label("Main") msymbol(D) mlcolor(magenta)) (`y'_nocontrol, label("No controls")), ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2)  onecell///
title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) ///
saving(figures/`y'_nocontrols, replace) level(99) legend(cols(3) span stack)
graph export "figures/`y'_nocontrols.svg", replace
}

////////////////////////////////////////////////////////////////////////////////
// Exogeneity of exposures
////////////////////////////////////////////////////////////////////////////////
eststo clear
eststo: est use output/robustness/nemadyr_exoboth
eststo: est use output/robustness/nattsyr_exoboth
eststo: est use output/robustness/over65_exoboth
esttab est1 est2 est3 using "output/robustness/exogeneity.rtf", /// 
star(* 0.01  ** 0.001 *** 0.0001 ) stats(N mean) ci replace ///
mtitle("Nr of Emergency Admissions in prev. year" "Nr of ED attendances in prev. year" "Aged 65 or above")

////////////////////////////////////////////////////////////////////////////////
// Measurement error
////////////////////////////////////////////////////////////////////////////////
eststo clear
eststo: est use output/robustness/misclass_error
esttab using "output/robustness/error.rtf", star(* 0.01  ** 0.001 *** 0.0001 ) stats(N r2_a) replace 

////////////////////////////////////////////////////////////////////////////////
// Supply side 1: Reverse causality in defining avoidable/non-avoidable patients
////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $after_out $bin_out $cont_out {
local vla : variable label `y'	
estimates use "output/main/`y'_main.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_tammes.ster"
eststo `y'_tammes
estimates use "output/robustness/`y'_ambdef.ster"
mat b=e(b)
mat colname b= "est_dailyvol_avoid" "est_dailyvol_nonavoid" ""
erepost b=b, rename
eststo `y'_ambdef
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_tammes, msymbol(T) label("Tammes"))  /// 
(`y'_ambdef, msymbol(O) label("Ambulance")), title("`vla'", size(medium)) ///
coeflabels(est_dailyvol_avoid = "Unexp. volume of avoid." est_dailyvol_nonavoid = "Unexp. volume of non-avoid." ) ///
level(99) legend(cols(3)) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est*) saving(figures/`y'_altdef, replace)
graph export "figures/`y'_altdef.png", replace
}

////////////////////////////////////////////////////////////////////////////////
// Supply side 2: Heterogeneity across hospital sizes
////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $after_out $bin_out $cont_out {
local vla : variable label `y'	
estimates use "output/main/`y'_main.ster"
eststo `y'_linear
foreach i in 1 2 3 {
estimates use "output/robustness/`y'_size`i'.ster"
eststo `y'size`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'size1, msymbol(T) label("Tertile 1 total vol."))  /// 
(`y'size2, msymbol(O) label("Tertile 2 total vol.")) (`y'size3, msymbol(+) label("Tertile 3 total vol.")), ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) ///
saving(figures/`y'_size, replace) level(99) legend(cols(3) span stack)
graph export "figures/`y'_size.png", replace
}

////////////////////////////////////////////////////////////////////////////////
// Supply side 3: Heterogeneity between times of day (in-hours/out-of-hours)
////////////////////////////////////////////////////////////////////////////////  

// Coefplot comparing main, oh, ih
eststo clear
foreach y in $after_out $bin_out $cont_out {
local vla : variable label `y'	
estimates use "output/main/`y'_main.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_ooh.ster"
eststo `y'_oh
estimates use "output/robustness/`y'_ih.ster"
eststo `y'_ih
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_oh, label("Out-of-hours"))  /// 
(`y'_ih, msymbol(O) label("In-hours")), title("`vla'", size(medium)) level(99) legend(cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) label ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) saving(figures/`y'_time, replace)
graph export "figures/`y'_time.png", replace
}

////////////////////////////////////////////////////////////////////////////////
// Supply side 4: Heterogeneity across quintiles of hospital bed occupation
////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $bin_out $cont_out $after_out {
local vla : variable label `y'	
estimates use "output/main/`y'_main.ster"
eststo `y'_linear
forvalues i=1/5 {
estimates use "output/robustness/`y'_bed`i'.ster"
eststo `y'_bed`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_bed1, msymbol(T) label("Quintile 1 bed. occ."))  /// 
(`y'_bed2, msymbol(O) label("Quintile 2 bed. occ.")) (`y'_bed3, msymbol(S) label("Quintile 3 bed. occ.")) ///
(`y'_bed4, msymbol(+) label("Quintile 4 bed. occ.")) (`y'_bed5, msymbol(X) label("Quintile 5 bed. occ.")), ///
 title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*1) ///
saving(figures/`y'_bed, replace) level(99)
graph export "figures/`y'_bed.png", replace
}

////////////////////////////////////////////////////////////////////////////////
// Supply side 6: Exclude central London where ambulances may decide where to go
////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $after_out $bin_out $cont_out {
local vla : variable label `y'	
estimates use "output/main/`y'_main.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_nolondon.ster"
eststo `y'_nolondon
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_nolondon, label("Excl. London")), ///
 title("`vla'", size(medium)) level(99) mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) saving(figures/`y'_nolondon, replace)
graph export "figures/`y'_nolondon.png", replace
}

////////////////////////////////////////////////////////////////////////////////
// Supply side 7: Analysis exploiting hourly variation
////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $cont_out $bin_out $after_out {
eststo: estimates use "output/robustness/`y'_hourly.ster"
}
esttab using "output/robustness/hourly.rtf", ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) replace ci(%9.4f) b(%9.4f) label scalar(N r2_a mean) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons) keep(est_*) 

////////////////////////////////////////////////////////////////////////////////
// Demand side 1: Heterogeneity across age groups
////////////////////////////////////////////////////////////////////////////////

gen agegr=1 if arrivalage>0&arrivalage<10
replace agegr=2 if arrivalage>10&arrivalage<20
replace agegr=3 if arrivalage>19&arrivalage<30
replace agegr=4 if arrivalage>29&arrivalage<60
replace agegr=5 if arrivalage>59&arrivalage<75
replace agegr=6 if arrivalage>74&arrivalage!=.

eststo clear
foreach y in $cont_out $bin_out $after_out  {
local vla : variable label `y'	
estimates use "output/main/`y'_main.ster"
eststo `y'_linear
forvalues i=1/6 {
estimates use "output/robustness/`y'_agegr`i'.ster"
eststo `y'_agegr`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_agegr1, msymbol(T) label("Aged 0-9"))  /// 
(`y'_agegr2, msymbol(O) label("Aged 10-19")) (`y'_agegr3, msymbol(+) label("Aged 20-29")) (`y'_agegr4, msymbol(S) label("Aged 30-59")) ///
(`y'_agegr5, msymbol(X) label("Aged 60-74")) (`y'_agegr6, msymbol(V) label("Aged 74+")), title("`vla'", size(medium)) legend(cols(4) stack) ///
mlabel format(%6.5f) mlabsize(vsmall) mlabposition(12) mlabgap(*1) level(99) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) saving(figures/`y'_agegr, replace)
graph export "figures/`y'_agegr.png", replace
}

////////////////////////////////////////////////////////////////////////////////
// Demand side 2: Checking unobserved severity across discharge destinations
////////////////////////////////////////////////////////////////////////////////

// Coefplot comparing main with disch destinations
eststo clear
foreach y in $cont_out tot4hr $after_out {
local vla : variable label `y'	
estimates use "output/main/`y'_main.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_lwt.ster"
eststo `y'_lwt
estimates use "output/robustness/`y'_nof.ster"
eststo `y'_nof
estimates use "output/robustness/`y'_gpfu.ster"
eststo `y'_gpfu
estimates use "output/robustness/`y'_refprov.ster"
eststo `y'_refprov
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_lwt, msymbol(T) label("Left without treatment")) ///
(`y'_nof, msymbol(O) label("Disch. with no follow-up")) (`y'_gpfu, msymbol(+) label("Disch. GP follow-up")) ///
(`y'_refprov, msymbol(X) label("Referred to clinic/other prov.")), ///
title("`vla'", size(medium)) level(99) legend(rows(2) cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) saving(figures/`y'_disch, replace)
graph export "figures/`y'_disch.png", replace
}

////////////////////////////////////////////////////////////////////////////////////////
// Demand side 3: Checking coefficient stability to adding more detailed severity controls 
////////////////////////////////////////////////////////////////////////////////////////
       
eststo clear
foreach y in $after_out tot4hr $cont_out {
local vla : variable label `y'	
estimates use "output/main/`y'_main.ster"
eststo `y'_linear
forvalues i=1/3 {
estimates use "output/robustness/`y'_stab`i'.ster"
eststo `y'_stab`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_stab1, msymbol(S) label("Without Elix. Comorb."))  /// 
(`y'_stab2, msymbol(O) label("With Elix. Comorb.")) (`y'_stab3, msymbol(X) label("With Elix. Comorb. + ICD10 code")) ///
, title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) level(99) legend(cols(2) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
msymbol(S) saving(figures/`y'_stab, replace)
graph export "figures/`y'_stab.png", replace
}

////////////////////////////////////////////////////////////////////////////////////////
// Specification check: models with interaction
////////////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $bin_out $after_out $cont_out {
estimates use "output/robustness/`y'_interacted.ster"
eststo `y'_interacted
}
esttab using "output/robustness/interacted.rtf", /// 
keep(est_dailyvol_avoid est_dailyvol_nonavoid *inter*) replace ci(%9.4f) b(%9.4f) label scalar(N r2_a) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons)

////////////////////////////////////////////////////////////////////////////////////////
// Results for year 2017/2018
////////////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $bin_out $after_out ninv1 ntrt1 totdur initdur durinvtret {
estimates use "output/robustness/`y'_altyear.ster"
eststo `y'_linear
}
esttab using "output/robustness/altyear.rtf", /// 
keep(est_dailyvol_avoid est_dailyvol_nonavoid) replace ci(%9.4f) b(%9.4f) label scalar(N r2_a mean) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons) stats(N r2_a mean)

////////////////////////////////////////////////////////////////////////////////
*------------------------------SUMMARY STATISTICS------------------------------*
////////////////////////////////////////////////////////////////////////////////

qui: areg admit dailyvol_avoid dailyvol_nonavoid ${x2} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
gen sample=e(sample)

// Gen descriptives for outcomes for Table 1
eststo clear
eststo: estpost sum $cont_out $bin_out $after_out if sample==1
esttab using "output/averages.rtf", replace cells("mean(fmt(3)) sd(fmt(3))") compress

// Create characteristics for descriptives
gen age015=0
replace age015=1 if arrivalage>=0&arrivalage<16
gen age1630=0
replace age1630=1 if arrivalage>=16&arrivalage<31
gen age3145=0
replace age3145=1 if arrivalage>=31&arrivalage<46
gen age4670=0
replace age4670=1 if arrivalage>=46&arrivalage<71
gen age71p=0
replace age71p=1 if arrivalage>=71&arrivalage<119
qui: tab eth, gen(deth)
qui: tab imd2015_cat, gen(dimd)
qui: tab nelixyr, gen(delixn)
gen attprevy0=0
replace attprevy0=1 if nattsyr==0
gen attprevy1=0
replace attprevy1=1 if nattsyr==1
gen attprevy2=0
replace attprevy2=1 if nattsyr>1&nattsyr!=.
gen refs0=0
replace refs0=1 if refsource==0
gen refs1=0
replace refs1=1 if refsource==1
gen refs2=0
replace refs2=1 if refsource==3
gen refs3=0
replace refs3=1 if refsource==4
gen refs4=0
replace refs4=1 if refsource==6
gen refs5=0
replace refs5=1 if refsource==7
gen refs6=0
replace refs6=1 if refsource==8
qui: tab bedocc_yday_q, gen(d_bedocc)
qui: tab attcat, gen(dattcatt)
qui: tab patgrp, gen(dpatgrp)

global summaryvars age015-age71p deth* dimd*  nattsyr nemadyr nelixyr refs1-refs6 ambarrival dattcatt1 dattcatt2 dpatgrp1-dpatgrp7 bedocc_yday dailyvol_avoid dailyvol_nonavoid

label var age015 "Age 0 to 15"
label var age1630 "Age 16 to 30"
label var age3145 "Age 31 to 45"
label var age4670 "Age 46 to 70"
label var age71p "Age over 70"
label var deth1 "White"
label var deth2 "Any Asian bckgr."
label var deth3 "Any Black bckgr."
label var deth4 "Chinese + Other"
label var deth5 "Missing/N.A."
label var dimd1 "Inc. depr. quint. 1"
label var dimd2 "Inc. depr. quint. 2"
label var dimd3 "Inc. depr. quint. 3"
label var dimd4 "Inc. depr. quint. 4"
label var dimd5 "Inc. depr. quint. 5"
label var ambarrival "Ambulance arrival"
label var delixn1 "0 or 1 morbidities"
label var delixn2 "2 morbidities"
label var delixn3 "3 morbidities"
label var delixn4 "4 morbidities"
label var delixn5 "5 morbidities"
label var delixn6 "6 morbidities"
label var delixn7 "7 morbidities"
label var delixn8 "8+ morbidities"
label var refs0 "GP"
label var refs1 "Self-referral"
label var refs2 "Emergency services"
label var refs3 "Work"
label var refs4 "Police"
label var refs5 "Health care provider"
label var refs6 "Other"
label var attprevy0 "No ED attendances in previous year"
label var attprevy1 "1 ED attendance in previous year"
label var attprevy2 "2 or more ED attendances in previous year"
label var dattcatt1 "First ED attendance"
label var dattcatt2 "Unplanned ED follow-up"
label var dpatgrp1 "Road traffic accident"
label var dpatgrp2 "Assault"
label var dpatgrp3 "Self-harm"
label var dpatgrp4 "Sports"
label var dpatgrp5 "Other accidents"
label var dpatgrp6 "Other"
label var dpatgrp7 "Not known"
eststo clear
eststo: estpost sum $summaryvars if sample==1
eststo: estpost sum $summaryvars if unexp_avoid_na_q==1&sample==1
eststo: estpost sum $summaryvars if unexp_avoid_na_q==5&sample==1
esttab using "output/summary_short.rtf", replace cells("mean(fmt(3)) sd(fmt(3))") label mtitles("Full NAV" "Q1 Unexp. AV" "Q5 Unexp. AV") compress

eststo clear
forvalues i=1/5 {
eststo: estpost sum $summaryvars if sample==1&unexp_avoid_na_q==`i'
}
esttab using "output/summary_full.rtf", replace cells("mean(fmt(3)) sd(fmt(3))") label mtitles("Q1 Unexp. AV" "Q2 Unexp. AV" "Q3 Unexp. AV" "Q4 Unexp. AV" "Q5 Unexp. AV") compress

fre prim_diag if sample==1

// Defined as avoidable across quintiles of volume
areg dailyvol_total, absorb(prov_dow_fe)
predict unexptot, residual
xtile qtls=unexptot, nq(5)
eststo clear
eststo: mean nhstreat, over(qtls) vce(cluster prov_id)
eststo: mean nhsinvest, over(qtls) vce(cluster prov_id)
eststo: mean nhsdisch, over(qtls) vce(cluster prov_id)
eststo: mean nhsarrmode, over(qtls) vce(cluster prov_id)
eststo: mean nhs_avoidable, over(qtls) vce(cluster prov_id)
coefplot (est1, asequation(Treatment) \ est2, asequation(Investigations) \ est3, asequation(Discharge) \ ///
est4, asequation(Arrival mode) \ est5, asequation(Avoidable) \ ), ///
coeflabels(*1.qtls = "Q1" *2.qtls = "Q2" *3.qtls = "Q3" *4.qtls = "Q4" *5.qtls = "Q5") ///
yline(0) ytitle("Probability of meeting criteria") xtitle("Quintile of unexpected attendance volume") ///
xlabel(,labsize(small)) xlabel(,grid) ci vertical saving(figures/probavoidable.gph, replace)
graph export figures/probavoidable.png, replace

// Distribution stats
set scheme s1mono
qui: sum dailyvol_tot if sample
hist dailyvol_tot if sample, normal saving(figures/totvol.gph, replace) percent xtitle("ED-level daily attendances") ///
note("Mean=`:display %05.3f `r(mean)'', SD=`:display %05.3f `r(sd)''") title("")
graph export figures/hvol_tot.png, replace
qui: sum unexp_avoid_na  if sample
hist unexp_avoid_na if sample, normal saving(figures/hunexp_avoid.gph, replace) percent xtitle("ED-level daily unexpected attendances") ///
note("Mean=`:display %05.3f `r(mean)'', SD=`:display %05.3f `r(sd)''") title("")
graph export figures/hunexp_avoid.png, replace
qui: sum unexp_nonavoid_na if sample
hist unexp_nonavoid_na if sample, normal saving(figures/hunexp_nonavoid.gph, replace) percent xtitle("ED-level daily unexpected attendances") ///
note("Mean=`:display %05.3f `r(mean)'', SD=`:display %05.3f `r(sd)''") title("")
graph export figures/hunexp_nonavoid.png, replace
qui: sum dailyvol_total if sample
hist dailyvol_avoid if sample, normal saving(figures/hvol_avoid.gph, replace) percent xtitle("ED-level daily attendances") ///
note("Mean=`:display %05.3f `r(mean)'', SD=`:display %05.3f `r(sd)''") title("")
graph export figures/hvol_avoid.png, replace
label var unexp_avoid_na "Unexp. volume of avoid." 
label var unexp_nonavoid_na "Unexp. volume of non-avoid."
sum dailyvol_total  unexp_dailyvol_tot unexp_avoid_na unexp_nonavoid_na if sample
corr unexp_avoid_na unexp_nonavoid_na if sample
scatter unexp_avoid_na unexp_nonavoid_na if sample, saving(figures/scatter, replace) ///
note("Correlation=`:display %05.3f `r(rho)''") title("") 
graph export figures/scatter.png, replace

// Seasonality plot
preserve
gen week=week(attdate)
bysort prov_id month: egen avgmonth=mean(dailyvol_total)
bysort prov_id week: egen avgweek=mean(dailyvol_total)
keep if prov_id==19
label var dailyvol_total "Daily volume"
label var avgweek "Weekly seasonal component"
label var avgmonth "Monthly seasonal component"
duplicates drop attdate, force
tsset attdate
twoway(tsline dailyvol_total) (tsline avgweek, lcolor(red)) (tsline avgmonth, lcolor(blue)), ///
saving(figures/dailyweekly, replace) yline(0)  tlabel(,angle(45))
graph export "figures/tsseaonality.png", replace
restore

// Unexpected pattern for hospital 19
bysort prov_id: egen avgdaily=mean(dailyvol_total)
sum avgdaily, det
label var unexp_avoid_na "Unexp. daily avoid. att."
label var attdate "Date"
preserve
keep if prov_id==19
duplicates drop attdate, force
tsset attdate
tsline unexp_avoid_na, title("") saving(figures/dailyshocks_p99, replace) yline(0)  tlabel(,angle(45))
graph export "figures/dailyshocks_p99.png", replace
restore
