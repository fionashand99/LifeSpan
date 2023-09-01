********************************************************************************
* Count models - output for publication supplementary Table                    *
********************************************************************************

global path "«path»"
cd "${path}"
* Open log file
capture log close
set linesize 132
log using "«log file name»", replace smcl
use "«data file name», clear 

/* Variables relevant to analysis 

site - LifeSpan intervention site; NSW coded as 5
time - Time in months from baseline
month_sin - sin(month-of-year-9)
P - phase of intervention 

sex - Male/Female of record
n - number of events 
erp10e5 - Estimated Resident Population in units of 100,000

n_total - number of events for males+females  
erp_total10e5 - ERP for males+females in units of 100,000

*/

******************************************************************************************
**# SETUP
******************************************************************************************
* Extract site labels as globals for labelling
forvalues L=1/5 {
global site_ID_`L' : label (site) `L'
}

* Extract phase labels as globals for labelling
forvalues p=0/6 {
global P_`L' : label P_label `p'
}

* Create binary indicator - in/not in (pre or post) intervention phase
recode P (2/4=2), gen(P_collapse)
label define P_collapse_label 0 "Pre" 1 "Trans" 2 "Intervention" 5 "Post 1" 6 "Post 2"
label values P_collapse  P_collapse_label 
label variable P_collapse "Phase (collapsed)"

* End of standard setup

******************************************************************************************
**# MODELLING INDIVIDUAL SITES USING 8-MONTH PERIODS AND TRANSITION
* - GENERATES SUPPLEMENTARY TABLE 1
******************************************************************************************
* "Final" model - just for checking
nbreg n ibn.site i.site#c.time c.month_sin   i.P  if site<5, exposure(erp10e5) nolog noconstant irr
estimates store model_all

* Final model in each site
local sex_labels `""Male" "Female""'
forvalues site = 1/5 {
forvalues sx =1/2 {
local sex_label : word `sx' of `sex_labels'
	di _n "************************* ${site_ID_`site'}: `sex_label' ***********************"
nbreg n time c.month_sin i.P if site==`site'&sex==`sx', exposure(erp10e5)  irr difficult
estimates store model_site`site'_`sex_label'
}
}

* Suppl Table 1
esttab model_site1_Male model_site1_Female model_site2_Male model_site2_Female model_site3_Male ///
       model_site3_Female model_site4_Male model_site4_Female model_site5_Male model_site5_Female  ///
	   using "«Supplementary table»" ///
       ,b(3) not  ci(3) eform compress tab type replace constant	   
	   
******************************************************************************************
**# TEST PRIMARY HYPOTHESIS - BOTTOM ROW OF SUPPL TABLE
******************************************************************************************
* Generates entries for "Intervention overall" cells of Suppl Table 1
local sex_labels `""Male" "Female""'
forvalues site = 1/5 {
forvalues sx =1/2 {
local sex_label : word `sx' of `sex_labels'
	di _n "************************* ${site_ID_`site'}: `sex_label' ***********************"
estimates restore model_site`site'_`sex_label'
quietly: nbreg, irr
* Test of 1y hypo - average of three intervention phases == pre-intervention
quietly: lincom _b[2.P]/3+_b[3.P]/3+_b[4.P]/3-_b[0.P]
di "par=" %6.3f r(estimate) ", 95% CI: " %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) ", p=" %6.3f r(p)

quietly: lincom _b[2.P]/3+_b[3.P]/3+_b[4.P]/3-_b[0.P],irr
di "IRR=" %6.3f r(estimate) ", 95% CI: " %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) ", p=" %6.3f r(p)
di as result  %6.3f r(estimate) "(" %5.3f  r(lb) "-" %5.3f r(ub) ")"
scalar ch = 100*(1-r(estimate))
di as result "% change=" %5.2f ch "%"
}
}
scalar drop ch

log close

