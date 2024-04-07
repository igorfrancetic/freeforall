////////////////////////////////////////////////////////////////////////////////
// Do-file defining avoidable patients and counting their volumes
// Follows the do-file named preliminary.do
// Version: 28/02/2024
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Define avoidable attendances according to NHS definition 
////////////////////////////////////////////////////////////////////////////////

/*
Conditions to apply following NHS Digital
See https://github.com/nhsengland/ESA_Avoidable_ED_Attendances
**************
INVEST_N: 06 = Urinalysis=, 21 = Pregnancy test, 22 = Dental investigation, or 24 = None;
TREAT_N: 221 = Guidance/advice only - written, 222 = Guidance/advice only - verbal, 30 = Recording vital signs, /// 
56 = Dental treatment, 57 = Prescription/medicines prepared to take away, 99 = None, or 07 = Prescriptions;
AEATTENDDISP: 02 = Discharged - follow-up treatment to be provided by general practitioner, /// 
03 = Discharged - did not require any follow-up treatment, or 12 = Left department before being treated;
AEARRIVALMODE: all codes except for 1 = non-ambulance arrivals. 
*/

// Here I implement these conditions in our dataset

*Flag treatments
foreach x in 1 2 3 4 5 6 7 8 9 10 11 12 {
                gen taga_`x'=inlist(treat_`x', "99", "", ".", "22", "221", "222", "57", "56", "30")
                gen tagb_`x'=inlist(treat_`x', "07", "999", ".30")
}
egen row=rsum(taga*)
egen crow=rsum(tagb*)
drop tag*
gen nhstreat=0
replace nhstreat=1 if row+crow==12
drop *row*

*Flag investigations
foreach x in 1 2 3 4 5 6 7 8 9 10 11 12 {
                gen tag_`x'=inlist(invest_`x', "06", "", ".", "21", "22", "24")
}
egen row=rsum(tag*)
drop tag*
gen nhsinvest=0
replace nhsinvest=1 if row>0&row==12
drop *row*

*Flag discharge destinations
gen nhsdisch=inlist(disp1a, 2, 3, 12)

*Flag arrival model
gen nhsarrmode=1 if arrmode!=.
replace nhsarrmode=0 if arrmode==1

*Create NHS definition
gen nhs_avoidable=nhstreat*nhsinvest*nhsdisch*nhsarrmode
gen nhs_nonavoidable=1 if nhs_avoidable!=.
replace nhs_nonavoidable=0 if nhs_avoidable==1

// Compute daily volumes by hospital and day
bysort prov_id arrivaldate: egen dailyvol_avoid=total(nhs_avoidable)
bysort prov_id arrivaldate: egen dailyvol_nonavoid=total(nhs_nonavoidable)
gen one=1
bysort prov_id arrivaldate: egen dailyvol_total=total(one)

// Generate daily unexpected volumes of avoidable and non-avoidable patients for graphs etc
areg dailyvol_avoid if nhs_avoidable==0, absorb(prov_dow_fe)
predict exp_dailyvol_avoid_na, xbd
predict unexp_avoid_na, residual
areg dailyvol_nonavoid if nhs_avoidable==0, absorb(prov_dow_fe)
predict exp_nondailyvol_avoid_na, xbd
predict unexp_nonavoid_na, residual
corr unexp_avoid_na unexp_nonavoid_na
areg dailyvol_total if nhs_avoidable==0, absorb(prov_dow_fe)
predict unexp_dailyvol_tot, residual
su unexp_dailyvol_tot, detail

xtile unexp_avoid_na_q=unexp_avoid_na, nq(5)
xtile unexp_nonavoid_na_q=unexp_nonavoid_na, nq(5)
xtile unexp_dailyvol_q=unexp_dailyvol_tot, nq(5)

// Compute the same but for hour of day
egen day=group(attdate)
egen how=group(hour dow)
egen clusterhour=group(prov_id how)
egen hod=group(attdate hour)

// Compute daily volumes for av and nav patients
bysort prov_id arrivaldate hour: egen hourvol_avoid=total(nhs_avoidable)
bysort prov_id arrivaldate hour: egen hourvol_nonavoid=total(nhs_nonavoidable)

////////////////////////////////////////////////////////////////////////////////
// Define avoidable attendances according to Tammes et al. (2016) definition
// Follows https://doi.org/10.1370/afm.2136
////////////////////////////////////////////////////////////////////////////////

/*AEREFSOURCE 01 = Self-referral;
AEATTENDDISP 03 = Discharged  did not require any follow-up treatment, 12 =
Left department before being treated, 13 = Left department having refused
treatment, or 02 = Discharged  follow-up treatment to be provided by general
practitioner.*/

*Flag discharge destinations
gen tammes_disch=inlist(disp1a, 2, 3, 12, 13)

*Flag arrival model
gen tammes_refsource=0 if refsource!=.
replace tammes_refsource=1 if refsource==1

*Create Tammes definition
gen tammes_avoidable=tammes_disch*tammes_refsource
gen tammes_nonavoidable=1 if tammes_avoidable!=.
replace tammes_nonavoidable=0 if tammes_avoidable==1

// Compute daily volumes by hospital and day
bysort prov_id attdate: egen dailyvol_tamav=total(tammes_avoidable)
bysort prov_id attdate: egen dailyvol_tamnon=total(tammes_nonavoidable)
areg dailyvol_tamav if tammes_avoidable==0, absorb(prov_dow_fe)
predict exp_tamav, xbd
predict unexp_tamav, residual
areg dailyvol_tamnon if tammes_avoidable==0, absorb(prov_dow_fe)
predict exp_tamnon, xbd
predict unexp_tamnon, residual
corr unexp_tamav unexp_tamnon

////////////////////////////////////////////////////////////////////////////////
// Define avoidable attendances based on ambulance arrival ***
////////////////////////////////////////////////////////////////////////////////

gen ambarrival=0
replace ambarrival=1 if arrmode==1

// Compute daily volumes by hospital and day
bysort prov_id attdate: egen dailyvol_nonavoidamb=total(ambarrival)
gen dailyvol_avoidamb=dailyvol_total-dailyvol_nonavoidamb
areg dailyvol_avoidamb if ambarrival==1, absorb(prov_dow_fe)
predict exp_avoidamb, xbd
predict unexp_avoidamb, residual
areg dailyvol_nonavoidamb if ambarrival==1, absorb(prov_dow_fe)
predict exp_nonavoidamb, xbd
predict unexp_nonavoidamb, residual
corr unexp_nonavoidamb unexp_avoidamb

////////////////////////////////////////////////////////////////////////////////
// Save final version ready for analysis.do
////////////////////////////////////////////////////////////////////////////////
save finaldata.dta, replace
