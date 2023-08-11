////////////////////////////////////////////////////////////////////////////////
// Do-file producing outputs for the paper
// Version: 08/08/2023
///////////////////////////////////////////////////////////////////////////////

/* Set paths */
if c(username)=="Add your username" {
	global path "Add your working directory"
}
cd "$path"
set scheme s1mono

// Open main dataset
use finaldata, clear

// Set global for covariates
global elixyear nchf2yr narrhy2yr nvd2yr npcd2yr npvd2yr nhptn_uc2yr nhptn_c2yr npara2yr nothnd2yr ncopd2yr ndiab_uc2yr ndiab_c2yr nhypothy2yr nrf2yr nld2yr npud_nb2yr nhiv2yr nlymp2yr nmets2yr ntumor2yr nrheum_a2yr ncoag2yr nobesity2yr nwl2yr nfluid2yr nbla2yr nda2yr nalcohol2yr ndrug2yr npsycho2yr ndep2yr
global x1 i.agesex ib0.eth i.imd2015_cat ib41.prim_diag ib0.nattsyr ib0.nemadyr ${elixyear} ib0.nelixyr ib2.arrmode i.attcat ib80.patgrp ib1.refsource ib10.incloc ib1.bedocc_yday_q  

// Set global for outcomes
global cont_out aecost ninv1 ntrt1 totdur initdur durinvtret
global bin_out tot4hr admit lwt disch_nof disch_gp refclin refoth 
global after_out reatt7_all reatt30_all dead7_ons dead30_ons 
global tret_out aecost ninv1 ntrt1 
global wt_out totdur durinvtret initdur
global binno4_out admit lwt disch_nof disch_gp refclin refoth 

global restore est_dailyvol_avoid est_dailyvol_nonavoid agesex eth imd2015_cat prim_diag nattsyr nemadyr ${elixyear} nelixyr arrmode attcat patgrp refsource incloc bedocc_yday_q prov_dow_fe


// Generate variables needed to generate outputs
gen est_dailyvol_avoid=dailyvol_avoid 
gen est_dailyvol_nonavoid=dailyvol_nonavoid 
gen vol_inter=.

// Relabel a few variables
label var est_dailyvol_avoid "Unexp. volume of avoid." 
label var est_dailyvol_nonavoid "Unexp. volume of non-avoid."
label var vol_inter "Interaction AV/NAV"
label var aecost "ED reimbursement (£)"
label var durinvtret "Time between assessment and treatment"
label var reatt7_all "Re-attendance within 7 days"
label var reatt30_all "Re-attendance within 30 days"
label var dead7_ons "Death within 7 days"
label var dead30_ons "Death within 30 days"
label var refoth "Referred to other professional or provider"

////////////////////////////////////////////////////////////////////////////////
// Outputs for main results and heterogeneity checks
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Main models, produce RTF tables (TeX also possible)
////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $bin_out {
estimates use "output/main/`y'_allhdfe.ster"
eststo `y'_linear
}
esttab using "tables/bin.rtf", /// 
keep(est_dailyvol_avoid est_dailyvol_nonavoid) replace ci(%9.4f) b(%9.4f) label scalar(N r2_a) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons)
       
eststo clear
foreach y in $after_out  {
estimates use "output/main/`y'_allhdfe.ster"
eststo `y'_linear
}
esttab using "tables/after.rtf", ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) replace ci(%9.4f) b(%9.4f) label scalar(N r2_a) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons)     

eststo clear
foreach y in $cont_out {
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_poisson
}
esttab using "tables/cont.rtf", replace ci(%9.4f) b(%9.4f) label scalar(N r2_a) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons)      

////////////////////////////////////////////////////////////////////////////////
// Models estimated on deciles of unexpected attendances, produce coefficient plots
////////////////////////////////////////////////////////////////////////////////

*Coefplots avoidable and nonavoidable
eststo clear
foreach y in $bin_out $after_out {
local vla : variable label `y'
estimates use "output/main/`y'_allhdfe_dec.ster"
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

*Coefplots avoidable and nonavoidable
eststo clear
foreach y in $cont_out {
local vla : variable label `y'
estimates use "output/main/`y'_poisson_dec.ster"
eststo `y'_poisson
coefplot (`y'_poisson, keep(*.dailyvol_avoid_q) label("Unexp. Avoidable") pstyle(p3)) ///
	 (`y'_poisson, keep(*.dailyvol_nonavoid_q) label("Unexp. Non-avoidable") pstyle(p4)), /// 
	 yline(0, lcolor(black) lwidth(thin) lpattern(dash)) vertical ///
	 mlabel format(%6.5f) mlabsize(vsmall) mlabposition(12) mlabgap(*8) ylab(, angle(45) labsize(small)) ///
	 ciopts(lwidth(*2)) levels(99) ///
	 xtitle("Att. volume shock decile (vs. Median)") ///
	 ytitle("Marginal effect on outcome") title("`vla'") ///
	 coeflabels(r1vs5.dailyvol_avoid_q = "1" r2vs5.dailyvol_avoid_q = "2" r3vs5.dailyvol_avoid_q = "3" r4vs5.dailyvol_avoid_q = "4" ///
	           r6vs5.dailyvol_avoid_q = "6" r7vs5.dailyvol_avoid_q = "7" r8vs5.dailyvol_avoid_q = "8" r9vs5.dailyvol_avoid_q = "9" ///
	           r10vs5.dailyvol_avoid_q = "10" r1vs5.dailyvol_nonavoid_q = "1" r2vs5.dailyvol_nonavoid_q = "2" r3vs5.dailyvol_nonavoid_q = "3" ///
	           r4vs5.dailyvol_nonavoid_q = "4" r6vs5.dailyvol_nonavoid_q = "6" r7vs5.dailyvol_nonavoid_q = "7" r8vs5.dailyvol_nonavoid_q = "8" ///
	           r9vs5.dailyvol_nonavoid_q = "9" r10vs5.dailyvol_nonavoid_q = "10")
graph save "figures/`y'_deccomb.gph", replace
graph export "figures/`y'_deccomb.png", replace
}

////////////////////////////////////////////////////////////////////////////////
// Model estimated across tertiles of NAV sharet, produce coefficient plots
////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $after_out $binno4_out {
local vla : variable label `y'	
estimates use "output/main/`y'_allhdfe.ster"
eststo `y'_linear
foreach i in 1 2 3 {
estimates use "output/main/`y'_ratio`i'.ster"
eststo `y'ratio`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'ratio1, label("Tertile 1 NA share")) /// 
(`y'ratio2, label("Tertile 2 NA share")) (`y'ratio3, label("Tertile 3 NA share")), level(99) legend(cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) msymbol(S) saving(figures/`y'_ratio, replace)
graph export "figures/`y'_ratio.png", replace
}

eststo clear
foreach y in $tret_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
foreach i in 1 2 3 {
estimates use "output/main/`y'_ratio`i'.ster"
eststo `y'ratio`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'ratio1, label("Tertile 1 NA share")) /// 
(`y'ratio2, label("Tertile 2 NA share")) (`y'ratio3, label("Tertile 3 NA share")), level(99) legend(cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) msymbol(S) saving(figures/`y'_ratio, replace)
graph export "figures/`y'_ratio.png", replace
}

eststo clear
foreach y in $wt_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
foreach i in 1 2 3 {
estimates use "output/main/`y'_ratio`i'.ster"
eststo `y'ratio`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'ratio1, label("Tertile 1 NA share")) /// 
(`y'ratio2, label("Tertile 2 NA share")) (`y'ratio3, label("Tertile 3 NA share")), level(99) legend(cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/`y'_ratio, replace)
graph export "figures/`y'_ratio.png", replace
}

eststo clear
local vla : variable label tot4hr
estimates use "output/main/tot4hr_allhdfe.ster"
eststo tot4hr_linear
foreach i in 1 2 3 {
estimates use "output/main/tot4hr_ratio`i'.ster"
eststo tot4hrratio`i'	
}
coefplot (tot4hr_linear, label("Main") msymbol(D) mlcolor(magenta)) (tot4hrratio1, label("Tertile 1 NA share AV"))  /// 
(tot4hrratio2, label("Tertile 2 NA share")) (tot4hrratio3, label("Tertile 3 NA share")), level(99) legend(cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/tot4hr_ratio, replace)
graph export "figures/tot4hr_ratio.png", replace

////////////////////////////////////////////////////////////////////////////////
// Outputs for robustness checks
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Supply side 1: Hospital size, produce coefficient plots
////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $after_out $binno4_out {
local vla : variable label `y'	
estimates use "output/main/`y'_allhdfe.ster"
eststo `y'_linear
foreach i in 1 2 3 {
estimates use "output/robustness/`y'_size`i'.ster"
eststo `y'size`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'size1, label("Tertile 1 total vol."))  /// 
(`y'size2, label("Tertile 2 total vol.")) (`y'size3, label("Tertile 3 total vol.")), ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) msymbol(S) saving(figures/`y'_size, replace) level(99) legend(cols(3) span stack)
graph export "figures/`y'_size.png", replace
}

eststo clear
foreach y in $tret_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
foreach i in 1 2 3 {
estimates use "output/robustness/`y'_size`i'.ster"
eststo `y'size`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'size1, label("Tertile 1 total vol."))  /// 
(`y'size2, label("Tertile 2 total vol.")) (`y'size3, label("Tertile 3 total vol.")), ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) msymbol(S) saving(figures/`y'_size, replace) level(99) legend(cols(3) span stack)
graph export "figures/`y'_size.png", replace
}

eststo clear
foreach y in $wt_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
foreach i in 1 2 3 {
estimates use "output/robustness/`y'_size`i'.ster"
eststo `y'size`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'size1, label("Tertile 1 total vol."))  /// 
(`y'size2, label("Tertile 2 total vol.")) (`y'size3, label("Tertile 3 total vol.")), ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/`y'_size, replace) level(99) legend(cols(3) span stack)
graph export "figures/`y'_size.png", replace
}

eststo clear
local vla : variable label tot4hr
estimates use "output/main/tot4hr_allhdfe.ster"
eststo tot4hr_linear
foreach i in 1 2 3 {
estimates use "output/robustness/tot4hr_size`i'.ster"
eststo tot4hrsize`i'	
}
coefplot (tot4hr_linear, label("Main") msymbol(D) mlcolor(magenta)) (tot4hrsize1, label("Tertile 1 total vol."))  /// 
(tot4hrsize2, label("Tertile 2 total vol.")) (tot4hrsize3, label("Tertile 3 total vol.")), ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/tot4hr_size, replace) level(99) legend(cols(3) span stack)
graph export "figures/tot4hr_size.png", replace

eststo clear
local vla : variable label lwt
estimates use "output/main/lwt_allhdfe.ster"
eststo lwt_linear
foreach i in 1 2 3 {
estimates use "output/robustness/lwt_size`i'.ster"
eststo lwt`i'	
}
coefplot (lwt_linear, label("Main") msymbol(D) mlcolor(magenta)) (lwt1, label("Tertile 1 total vol."))  /// 
(lwt2, label("Tertile 2 total vol.")) (lwt3, label("Tertile 3 total vol.")), ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/lwt_size, replace) level(99) legend(cols(3) span stack)
graph export "figures/lwt_size.png", replace

////////////////////////////////////////////////////////////////////////////////
// Supply side 2: Time of day (in-hours/out-of-hours), produce coefficient plots
////////////////////////////////////////////////////////////////////////////////  

// Coefplot comparing main, oh, ih
eststo clear
foreach y in $after_out $binno4_out {
local vla : variable label `y'	
estimates use "output/main/`y'_allhdfe.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_ooh.ster"
eststo `y'_oh
estimates use "output/robustness/`y'_ih.ster"
eststo `y'_ih
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_oh, label("Out-of-hours"))  /// 
(`y'_ih, label("In-hours")), title("`vla'", size(medium)) level(99) legend(cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) msymbol(S) saving(figures/`y'_time, replace)
graph export "figures/`y'_time.png", replace
}

eststo clear
foreach y in $tret_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_ooh.ster"
eststo `y'_oh
estimates use "output/robustness/`y'_ih.ster"
eststo `y'_ih
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_oh, label("Out-of-hours"))  /// 
(`y'_ih, label("In-hours")), title("`vla'", size(medium)) level(99) legend(cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) msymbol(S) saving(figures/`y'_time, replace)
graph export "figures/`y'_time.png", replace
}

eststo clear
foreach y in $wt_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_ooh.ster"
eststo `y'_oh
estimates use "output/robustness/`y'_ih.ster"
eststo `y'_ih
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_oh, label("Out-of-hours"))  /// 
(`y'_ih, label("In-hours")), title("`vla'", size(medium)) level(99) legend(cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/`y'_time, replace)
graph export "figures/`y'_time.png", replace
}

eststo clear
local vla : variable label tot4h	
estimates use "output/main/tot4hr_allhdfe.ster"
eststo tot4hr_linear
estimates use "output/robustness/tot4hr_ooh.ster"
eststo tot4hr_oh
estimates use "output/robustness/tot4hr_ih.ster"
eststo tot4hr_ih
coefplot (tot4hr_linear, label("Main") msymbol(D) mlcolor(magenta)) (tot4hr_oh, label("Out-of-hours"))  /// 
(tot4hr_ih, label("In-hours")), title("`vla'", size(medium)) level(99) legend(cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/tot4hr_time, replace)
graph export "figures/tot4hr_time.png", replace

////////////////////////////////////////////////////////////////////////////////
// Supply side 3: Defining avoidable attendances, produce coefficient plots
////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $after_out $binno4_out {
local vla : variable label `y'	
estimates use "output/main/`y'_allhdfe.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_tammes.ster"
eststo `y'_tammes
estimates use "output/robustness/`y'_ambdef.ster"
eststo `y'_ambdef
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_tammes, label("Tammes"))  /// 
(`y'_ambdef, label("Ambulance")), title("`vla'", size(medium)) ///
coeflabels(est_dailyvol_avoid = "Unexp. volume of avoid." est_dailyvol_nonavoid = "Unexp. volume of non-avoid." ///
est_dailyvol_avoidamb = "Unexp. volume of avoid." est_dailyvol_nonavoidamb = "Unexp. volume of non-avoid.") level(99) legend(cols(3)) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est*) xline(0) msymbol(S) saving(figures/`y'_altdef, replace)
graph export "figures/`y'_altdef.png", replace
}

eststo clear
foreach y in $tret_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_poisson_tammes.ster"
eststo `y'_tammes
estimates use "output/robustness/`y'_poisson_ambdef.ster"
eststo `y'_ambdef
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_tammes, label("Tammes"))  /// 
(`y'_ambdef, label("Ambulance")), title("`vla'", size(medium)) ///
coeflabels(est_dailyvol_avoid = "Unexp. volume of avoid." est_dailyvol_nonavoid = "Unexp. volume of non-avoid." ///
est_dailyvol_avoidamb = "Unexp. volume of avoid." est_dailyvol_nonavoidamb = "Unexp. volume of non-avoid.") level(99) legend(cols(3)) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est*) xline(0) msymbol(S) saving(figures/`y'_altdef, replace)
graph export "figures/`y'_altdef.png", replace
}

eststo clear
foreach y in $wt_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_poisson_tammes.ster"
eststo `y'_tammes
estimates use "output/robustness/`y'_poisson_ambdef.ster"
eststo `y'_ambdef
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_tammes, label("Tammes"))  /// 
(`y'_ambdef, label("Ambulance")), title("`vla'", size(medium)) ///
coeflabels(est_dailyvol_avoid = "Unexp. volume of avoid." est_dailyvol_nonavoid = "Unexp. volume of non-avoid." ///
est_dailyvol_avoidamb = "Unexp. volume of avoid." est_dailyvol_nonavoidamb = "Unexp. volume of non-avoid.") level(99) legend(cols(3)) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid est_dailyvol_avoidamb est_dailyvol_nonavoidamb) msymbol(S) saving(figures/`y'_altdef, replace)
graph export "figures/`y'_altdef.png", replace
}

eststo clear
local vla : variable label tot4hr	
estimates use "output/main/tot4hr_allhdfe.ster"
eststo tot4hr_linear
estimates use "output/robustness/tot4hr_tammes.ster"
eststo tot4hr_tammes
estimates use "output/robustness/tot4hr_ambdef.ster"
eststo tot4hr_ambdef
coefplot (tot4hr_linear, label("Main") msymbol(D) mlcolor(magenta)) (tot4hr_tammes, label("Tammes"))  /// 
(tot4hr_ambdef, label("Ambulance")), title("`vla'", size(medium)) ///
coeflabels(est_dailyvol_avoid = "Unexp. volume of avoid." est_dailyvol_nonavoid = "Unexp. volume of non-avoid." ///
est_dailyvol_avoidamb = "Unexp. volume of avoid." est_dailyvol_nonavoidamb = "Unexp. volume of non-avoid.") level(99) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid est_dailyvol_avoidamb est_dailyvol_nonavoidamb) msymbol(S) saving(figures/tot4hr_altdef, replace) legend(cols(3))
graph export "figures/tot4hr_altdef.png", replace

eststo clear
local vla : variable label lwt
estimates use "output/main/lwt_allhdfe.ster"
eststo lwt_linear
estimates use "output/robustness/lwt_tammes.ster"
eststo lwt_tammes
estimates use "output/robustness/lwt_ambdef.ster"
eststo lwt_ambdef
coefplot (lwt_linear, label("Main") msymbol(D) mlcolor(magenta)) (lwt_tammes, label("Tammes"))  /// 
(lwt_ambdef, label("Ambulance")), title("`vla'", size(medium)) ///
coeflabels(est_dailyvol_avoid = "Unexp. volume of avoid." est_dailyvol_nonavoid = "Unexp. volume of non-avoid." ///
est_dailyvol_avoidamb = "Unexp. volume of avoid." est_dailyvol_nonavoidamb = "Unexp. volume of non-avoid.") ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) level(99) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid est_dailyvol_avoidamb est_dailyvol_nonavoidamb) msymbol(S) saving(figures/lwt_altdef, replace) legend(cols(3))
graph export "figures/lwt_altdef.png", replace
       
////////////////////////////////////////////////////////////////////////////////
// Supply side 4: Exclude central London, produce coefficient plots
////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $after_out $bin_out {
local vla : variable label `y'	
estimates use "output/main/`y'_allhdfe.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_nolondon.ster"
eststo `y'_nolondon
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_nolondon, label("Excl. London")), title("`vla'", size(medium)) level(99)  /// 
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) msymbol(S) saving(figures/`y'_nolondon, replace)
graph export "figures/`y'_nolondon.png", replace
}

eststo clear
foreach y in $tret_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_nolondon.ster"
eststo `y'_nolondon
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_nolondon, label("Excl. London")), title("`vla'", size(medium)) level(99) /// 
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) msymbol(S) saving(figures/`y'_nolondon, replace)
graph export "figures/`y'_nolondon.png", replace
}

eststo clear
foreach y in $wt_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_nolondon.ster"
eststo `y'_nolondon
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_nolondon, label("Excl. London")), title("`vla'", size(medium)) level(99) /// 
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/`y'_nolondon, replace)
graph export "figures/`y'_nolondon.png", replace
}

eststo clear
local vla : variable label tot4hr	
estimates use "output/main/tot4hr_allhdfe.ster"
eststo tot4hr_linear
estimates use "output//robustness/tot4hr_nolondon.ster"
eststo tot4hr_nolondon
coefplot (tot4hr_linear, label("Main") msymbol(D) mlcolor(magenta)) (tot4hr_nolondon, label("Excl. London")), title("`vla'", size(medium)) level(99)  /// 
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/tot4hr_nolondon, replace)
graph export "figures/tot4hr_nolondon.png", replace

eststo clear
local vla : variable label lwt
estimates use "output/main/lwt_allhdfe.ster"
eststo lwt_linear
estimates use "output//robustness/lwt_nolondon.ster"
eststo lwt_nolondon
coefplot (lwt_linear, label("Main") msymbol(D) mlcolor(magenta)) (lwt_nolondon, label("Excl. London")), title("`vla'", size(medium)) level(99) /// 
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/lwt_nolondon, replace)
graph export "figures/lwt_nolondon.png", replace

////////////////////////////////////////////////////////////////////////////////
// Supply side 5: Quintiles of hospital bed occupation, produce coefficient plots
////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $binno4_out  $after_out {
local vla : variable label `y'	
estimates use "output/main/`y'_allhdfe.ster"
eststo `y'_linear
forvalues i=1/5 {
estimates use "output/robustness/`y'_bed`i'.ster"
eststo `y'_bed`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_bed1, label("Quintile 1 bed. occ."))  /// 
(`y'_bed2, label("Quintile 2 bed. occ.")) (`y'_bed3, label("Quintile 3 bed. occ.")) (`y'_bed4, label("Quintile 4 bed. occ."))  ///
(`y'_bed5, label("Quintile 5 bed. occ.")), title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*1) ///
msymbol(S) saving(figures/`y'_bed, replace) level(99)
graph export "figures/`y'_bed.png", replace
}

eststo clear
foreach y in $tret_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
forvalues i=1/5 {
estimates use "output/robustness/`y'_bed`i'.ster"
eststo `y'_bed`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_bed1, label("Quintile 1 bed. occ."))  /// 
(`y'_bed2, label("Quintile 2 bed. occ.")) (`y'_bed3, label("Quintile 3 bed. occ.")) (`y'_bed4, label("Quintile 4 bed. occ."))  ///
(`y'_bed5, label("Quintile 5 bed. occ.")), title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) level(99) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*1) ///
xline(0) msymbol(S) saving(figures/`y'_bed, replace)
graph export "figures/`y'_bed.png", replace
}

eststo clear
foreach y in $wt_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
forvalues i=1/5 {
estimates use "output/robustness/`y'_bed`i'.ster"
eststo `y'_bed`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_bed1, label("Quintile 1 bed. occ."))  /// 
(`y'_bed2, label("Quintile 2 bed. occ.")) (`y'_bed3, label("Quintile 3 bed. occ.")) (`y'_bed4, label("Quintile 4 bed. occ.")) ///
(`y'_bed5, label("Quintile 5 bed. occ.")), title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid)  ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*1) ///
msymbol(S) saving(figures/`y'_bed, replace) level(99)
graph export "figures/`y'_bed.png", replace
}

eststo clear
local vla : variable label tot4hr	
estimates use "output/main/tot4hr_allhdfe.ster"
eststo tot4hr_linear
forvalues i=1/5 {
estimates use "output/robustness/tot4hr_bed`i'.ster"
eststo tot4hr_bed`i'	
}
coefplot (tot4hr_linear, label("Main") msymbol(D) mlcolor(magenta)) (tot4hr_bed1, label("Quintile 1 bed. occ."))  /// 
(tot4hr_bed2, label("Quintile 2 bed. occ.")) (tot4hr_bed3, label("Quintile 3 bed. occ.")) (tot4hr_bed4, label("Quintile 4 bed. occ.")) ///
(tot4hr_bed5, label("Quintile 5 bed. occ.")), title("`vla'", size(medium))  /// 
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*1) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/tot4hr_bed, replace) level(99)
graph export "figures/tot4hr_bed.png", replace

////////////////////////////////////////////////////////////////////////////////
// Demand side 1: Heterogeneity across age groups, produce coefficient plots
////////////////////////////////////////////////////////////////////////////////

gen agegr=1 if arrivalage>0&arrivalage<10
replace agegr=2 if arrivalage>10&arrivalage<20
replace agegr=3 if arrivalage>19&arrivalage<30
replace agegr=4 if arrivalage>29&arrivalage<60
replace agegr=5 if arrivalage>59&arrivalage<75
replace agegr=6 if arrivalage>74&arrivalage!=.

eststo clear
foreach y in $binno4_out $after_out {
local vla : variable label `y'	
estimates use "output/main/`y'_allhdfe.ster"
eststo `y'_linear
forvalues i=1/6 {
estimates use "output/robustness/`y'_agegr`i'.ster"
eststo `y'_agegr`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_agegr1, label("Aged 0-9"))  /// 
(`y'_agegr2, label("Aged 10-19")) (`y'_agegr3, label("Aged 20-29")) (`y'_agegr4, label("Aged 30-59")) ///
(`y'_agegr5, label("Aged 60-74")) (`y'_agegr6, label("Aged 74+")), title("`vla'", size(medium)) legend(cols(4) stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*1) level(99) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) msymbol(S) saving(figures/`y'_agegr, replace)
graph export "figures/`y'_agegr.png", replace
}

eststo clear
foreach y in $tret_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
forvalues i=1/6 {
estimates use "output/robustness/`y'_agegr`i'.ster"
eststo `y'_agegr`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_agegr1, label("Aged 0-9"))  /// 
(`y'_agegr2, label("Aged 10-19")) (`y'_agegr3, label("Aged 20-29")) (`y'_agegr4, label("Aged 30-59")) ///
(`y'_agegr5, label("Aged 60-74")) (`y'_agegr6, label("Aged 74+")), title("`vla'", size(medium)) legend(cols(4) stack) /// 
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*1) level(99) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) msymbol(S) saving(figures/`y'_agegr, replace)
graph export "figures/`y'_agegr.png", replace
}

eststo clear
foreach y in $wt_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
forvalues i=1/6 {
estimates use "output/robustness/`y'_agegr`i'.ster"
eststo `y'_agegr`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_agegr1, label("Aged 0-9"))  /// 
(`y'_agegr2, label("Aged 10-19")) (`y'_agegr3, label("Aged 20-29")) (`y'_agegr4, label("Aged 30-59")) ///
(`y'_agegr5, label("Aged 60-74")) (`y'_agegr6, label("Aged 74+")), title("`vla'", size(medium)) legend(cols(4) stack) /// 
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*1) level(99) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/`y'_agegr, replace)
graph export "figures/`y'_agegr.png", replace
}

eststo clear
local vla : variable label tot4hr	
estimates use "output/main/tot4hr_allhdfe.ster"
eststo tot4hr_linear
forvalues i=1/6 {
estimates use "output/robustness/tot4hr_agegr`i'.ster"
eststo tot4hr_agegr`i'	
}
coefplot (tot4hr_linear, label("Main") msymbol(D) mlcolor(magenta)) (tot4hr_agegr1, label("Aged 0-9"))  /// 
(tot4hr_agegr2, label("Aged 10-19")) (tot4hr_agegr3, label("Aged 20-29")) (tot4hr_agegr4, label("Aged 30-59")) ///
(tot4hr_agegr5, label("Aged 60-74")) (tot4hr_agegr6, label("Aged 74+")), title("`vla'", size(medium)) legend(cols(4) stack) /// 
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*1) level(99) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/tot4hr_agegr, replace)
graph export "figures/tot4hr_agegr.png", replace

////////////////////////////////////////////////////////////////////////////////
// Demand side 2: Discharge destinations, produce coefficient plots
////////////////////////////////////////////////////////////////////////////////

// Coefplot comparing main with disch destinations
eststo clear
foreach y in $after_out {
local vla : variable label `y'	
estimates use "output/main/`y'_allhdfe.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_lwt.ster"
eststo `y'_lwt
estimates use "output/robustness/`y'_nof.ster"
eststo `y'_nof
estimates use "output/robustness/`y'_gpfu.ster"
eststo `y'_gpfu
estimates use "output/robustness/`y'_refprov.ster"
eststo `y'_refprov
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_lwt, label("Left without treatment"))  (`y'_nof, label("Disch. with no follow-up")) /// 
(`y'_gpfu, label("Disch. GP follow-up")) (`y'_refprov, msymbol(X) label("Referred to clinic/other prov."))  ///
, title("`vla'", size(medium)) level(99) legend(rows(2) cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) msymbol(S) saving(figures/`y'_disch, replace)
graph export "figures/`y'_disch.png", replace
}

eststo clear
foreach y in $tret_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_lwt.ster"
eststo `y'_lwt
estimates use "output/robustness/`y'_nof.ster"
eststo `y'_nof
estimates use "output/robustness/`y'_gpfu.ster"
eststo `y'_gpfu
estimates use "output/robustness/`y'_refprov.ster"
eststo `y'_refprov
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_lwt, label("Left without treatment"))  (`y'_nof, label("Disch. with no follow-up")) /// 
(`y'_gpfu, label("Disch. GP follow-up")) (`y'_refprov, msymbol(X) label("Referred to clinic/other prov."))  ///
, title("`vla'", size(medium)) level(99) legend(rows(2) cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) msymbol(S) saving(figures/`y'_disch, replace)
graph export "figures/`y'_disch.png", replace
}

eststo clear
foreach y in $wt_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
estimates use "output/robustness/`y'_lwt.ster"
eststo `y'_lwt
estimates use "output/robustness/`y'_nof.ster"
eststo `y'_nof
estimates use "output/robustness/`y'_gpfu.ster"
eststo `y'_gpfu
estimates use "output/robustness/`y'_refprov.ster"
eststo `y'_refprov
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_lwt, label("Left without treatment"))  (`y'_nof, label("Disch. with no follow-up")) /// 
(`y'_gpfu, label("Disch. GP follow-up")) (`y'_refprov, msymbol(X) label("Referred to clinic/other prov."))  ///
, title("`vla'", size(medium)) level(99) legend(rows(2) cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/`y'_disch, replace)
graph export "figures/`y'_disch.png", replace
}

local vla : variable label tot4hr	
estimates use "output/main/tot4hr_allhdfe.ster"
eststo tot4hr_linear
estimates use "output/robustness/tot4hr_lwt.ster"
eststo tot4hr_lwt
estimates use "output/robustness/tot4hr_nof.ster"
eststo tot4hr_nof
estimates use "output/robustness/tot4hr_gpfu.ster"
eststo tot4hr_gpfu
estimates use "output/robustness/tot4hr_refprov.ster"
eststo tot4hr_refprov
coefplot (tot4hr_linear, label("Main") msymbol(D) mlcolor(magenta)) (tot4hr_lwt, label("Left without treatment"))  (tot4hr_nof, label("Disch. with no follow-up")) (tot4hr_gpfu, label("Disch. GP follow-up")) /// 
(tot4hr_refprov, msymbol(X) label("Referred to clinic/other prov.")), ///
title("`vla'", size(medium)) level(99) legend(rows(2) cols(3) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) msymbol(S) saving(figures/tot4hr_disch, replace)
graph export "figures/tot4hr_disch.png", replace

////////////////////////////////////////////////////////////////////////////////////////
// Demand side 3: Coefficient stability, produce coefficient plots
////////////////////////////////////////////////////////////////////////////////////////
       
eststo clear
foreach y in $after_out {
local vla : variable label `y'	
estimates use "output/main/`y'_allhdfe.ster"
eststo `y'_linear
forvalues i=1/3 {
estimates use "output/robustness/`y'_stab`i'.ster"
eststo `y'_stab`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_stab1, label("Without Elix. Comorb."))  /// 
(`y'_stab2, label("With Elix. Comorb.")) (`y'_stab3, label("With Elix. Comorb. + ICD10 code")) ///
, title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) level(99) legend(cols(2) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
msymbol(S) saving(figures/`y'_stab, replace)
graph export "figures/`y'_stab.png", replace
}

eststo clear
foreach y in $tret_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
forvalues i=1/3 {
estimates use "output/robustness/`y'_stab`i'.ster"
eststo `y'_stab`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_stab1, label("Without Elix. Comorb."))  /// 
(`y'_stab2, label("With Elix. Comorb.")) (`y'_stab3, label("With Elix. Comorb. + ICD10 code")) ///
, title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) xline(0) level(99) legend(cols(2) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
msymbol(S) saving(figures/`y'_stab, replace)
graph export "figures/`y'_stab.png", replace
}

eststo clear
foreach y in $wt_out {
local vla : variable label `y'	
estimates use "output/main/`y'_poisson_allhdfe.ster"
eststo `y'_linear
forvalues i=1/3 {
estimates use "output/robustness/`y'_stab`i'.ster"
eststo `y'_stab`i'	
}
coefplot (`y'_linear, label("Main") msymbol(D) mlcolor(magenta)) (`y'_stab1, label("Without Elix. Comorb."))  /// 
(`y'_stab2, label("With Elix. Comorb.")) (`y'_stab3, label("With Elix. Comorb. + ICD10 code")) ///
, title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) level(99) legend(cols(2) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
msymbol(S) saving(figures/`y'_stab, replace)
graph export "figures/`y'_stab.png", replace
}

eststo clear
local vla : variable label tot4hr	
estimates use "output/main/tot4hr_allhdfe.ster"
eststo tot4hr_linear
forvalues i=1/3 {
estimates use "output/robustness/tot4hr_stab`i'.ster"
eststo tot4hr_stab`i'	
}
coefplot (tot4hr_linear, label("Main") msymbol(D) mlcolor(magenta)) (tot4hr_stab1, label("Without Elix. Comorb."))  /// 
(tot4hr_stab2, label("With Elix. Comorb.")) (tot4hr_stab3, label("With Elix. Comorb. + ICD10 code")) ///
, title("`vla'", size(medium)) keep(est_dailyvol_avoid est_dailyvol_nonavoid) level(99) legend(cols(2) span stack) ///
mlabel format(%6.5f) mlabsize(small) mlabposition(12) mlabgap(*2) ///
msymbol(S) saving(figures/tot4hr_stab, replace)
graph export "figures/tot4hr_stab.png", replace

////////////////////////////////////////////////////////////////////////////////////////
// Models including interactions, produce RTF tables (TeX also possible)
////////////////////////////////////////////////////////////////////////////////////////

eststo clear
foreach y in $bin_out {
estimates use "output/interacted/`y'_allhdfe.ster"
eststo `y'_linear
}
esttab using "tables/bin_interact.rtf", /// 
keep(est_dailyvol_avoid est_dailyvol_nonavoid) replace ci(%9.4f) b(%9.4f) label scalar(N r2_a) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons)
       
eststo clear
foreach y in $after_out  {
estimates use "output/interacted/`y'_allhdfe.ster"
eststo `y'_linear
}
esttab using "tables/after_interact.rtf", ///
keep(est_dailyvol_avoid est_dailyvol_nonavoid) replace ci(%9.4f) b(%9.4f) label scalar(N r2_a) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons)     

eststo clear
foreach y in $cont_out {
estimates use "output/interacted/`y'_poisson_allhdfe.ster"
eststo `y'_poisson
}
esttab using "tables/cont_interact.rtf", replace ci(%9.4f) b(%9.4f) label scalar(N r2_a) ///
star(* 0.01 ** 0.001 *** 0.0001) level(99) order(_cons) 


////////////////////////////////////////////////////////////////////////////////
// Various tables and figures of descriptive statistics throughout the paper
////////////////////////////////////////////////////////////////////////////////

// Define a sample of reference (= sample used in the analysis on probability of admission)
qui: areg admit est_dailyvol_avoid est_dailyvol_nonavoid ${x1} if nhs_nonavoidable==1, absorb(prov_dow_fe) vce(cluster prov_id)
gen sample=e(sample)

// Create variables for descriptives
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

// Label new variables
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

// Define global for all variables to summarise
global summaryvars age015-age71p deth* dimd*  nattsyr nemadyr nelixyr refs1-refs6 ambarrival dattcatt1 dattcatt2 dpatgrp1-dpatgrp7 bedocc_yday dailyvol_avoid dailyvol_nonavoid

// Generate Table 2
eststo clear
eststo: estpost sum $summaryvars if sample==1
eststo: estpost sum $summaryvars if unexp_avoid_na_q==1&sample==1
eststo: estpost sum $summaryvars if unexp_avoid_na_q==5&sample==1
esttab using "output/summary_short.rtf", replace cells("mean(fmt(3)) sd(fmt(3))") label mtitles("Full NAV" "Q1 Unexp. AV" "Q5 Unexp. AV") compress

// Generate Table A1
eststo clear
forvalues i=1/5 {
eststo: estpost sum $summaryvars if sample==1&unexp_avoid_na_q==`i'
}
esttab using "output/summary_full.rtf", replace cells("mean(fmt(3)) sd(fmt(3))") label mtitles("Q1 Unexp. AV" "Q2 Unexp. AV" "Q3 Unexp. AV" "Q4 Unexp. AV" "Q5 Unexp. AV") compress

// Generate Table A2
fre prim_diag if sample==1

// Generate Figure 1
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
twoway(tsline dailyvol_total) (tsline avgweek, lcolor(red)) (tsline avgmonth, lcolor(blue)), saving(figures/dailyweekly, replace) yline(0)  tlabel(,angle(45))
graph export "figures/tsseaonality.png", replace
restore

// Generate boxes for Figure 2
set scheme s1mono
qui: sum dailyvol_tot
hist dailyvol_tot, normal saving(figures/hunexp_tot.gph, replace) percent xtitle("ED-level daily attendances") note("Mean=`r(mean)', SD=`r(sd)'") title("")
graph export figures/hvol_tot.png, replace
qui: sum unexp_dailyvol_tot
hist unexp_avoid_na, normal saving(figures/hunexp_avoid.gph, replace) percent xtitle("ED-level daily unexpected attendances") note("Mean=`r(mean)', SD=`r(sd)'") title("")
graph export figures/hunexp_avoid.png, replace
qui: sum unexp_nonavoid_na
hist unexp_nonavoid_na, normal saving(figures/hunexp_nonavoid.gph, replace) percent xtitle("ED-level daily unexpected attendances") note("Mean=`r(mean)', SD=`r(sd)'") title("")
graph export figures/hunexp_nonavoid.png, replace
qui: sum dailyvol_total
hist dailyvol_avoid if nhs_nonavoidable==1, normal saving(figures/hvol_avoid.gph, replace) percent xtitle("ED-level daily attendances") note("Mean=`r(mean)', SD=`r(sd)'") title("")
graph export figures/hvol_avoid.png, replace
label var unexp_avoid_na "Unexp. volume of avoid." 
label var unexp_nonavoid_na "Unexp. volume of non-avoid."
sum dailyvol_total  unexp_dailyvol_tot unexp_avoid_na unexp_nonavoid_na
corr unexp_avoid_na unexp_nonavoid_na
scatter unexp_avoid_na unexp_nonavoid_na, saving(figures/scatter, replace) note("Correlation=`r(rho)'") title("") 
graph export figures/scatter.png, replace
** Unexpected pattern for hospital 19
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

// Generate Figure 5
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