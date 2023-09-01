********************************************************************************
* Count models - output for publication Table 2                                *
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

* Extract site labels as globals for labelling
forvalues L=1/5 {
global site_ID_`L' : label (site) `L'
}
label values P P_label 
* Extract phase labels as globals for labelling
forvalues p=0/6 {
global P_`L' : label P_label `p'
}

* Create binary indicator - in/not in (pre or post) intervention phase
recode P (2/4=2), gen(P_collapse)
label define P_collapse_label 0 "Pre" 1 "Trans" 2 "Intervention" 5 "Post 1" 6 "Post 2"
label values P_collapse  P_collapse_label 
label variable P_collapse "Phase (collapsed)"

******************************************************************************************
* MODELLING INDIVIDUAL SITES USING 8-MONTH PERIODS AND TRANSITION
******************************************************************************************

* "Final" model
nbreg n ibn.site i.site#c.time c.month_sin i.P  if site<5, exposure(erp10e5) nolog noconstant irr
estimates store model_all

* Final model in each site
forvalues site = 1/5 {
di _n "************************* ${site_ID_`site'} ***********************"
nbreg n ibn.site i.site#c.time c.month_sin i.P if site==`site', exposure(erp10e5) nolog noconstant irr 
estimates store model_site`site'
}

* Final model for males and females
di _n "************************* Males ***********************"
nbreg n ibn.site i.site#c.time c.month_sin i.P if site<5&sex==1, exposure(erp10e5) noconstant nolog irr 
estimates store model_Male
di _n "************************* Females ***********************"
nbreg n ibn.site i.site#c.time c.month_sin i.P if site<5&sex==2, exposure(erp10e5) noconstant nolog irr 
estimates store model_Female

******************************************************************************************
**# TEST PRIMARY HYPOTHESIS FOR ALL SITES
******************************************************************************************
di _n "************************* ALL SITES ***********************"
estimates restore model_all
`switch_quiet' nbreg, irr
* Test of 1y hypo - average of three intervention phases == pre-intervention
`switch_quiet' lincom _b[2.P]/3+_b[3.P]/3+_b[4.P]/3-_b[0.P]
di "par=" %6.3f r(estimate) ", 95% CI: " %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) ", p=" %6.4f r(p)

`switch_quiet' lincom _b[2.P]/3+_b[3.P]/3+_b[4.P]/3-_b[0.P],irr
di "IRR=" %6.3f r(estimate) ", 95% CI: " %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) ", p=" %6.4f r(p)
di  %6.3f r(estimate) "[" %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) "]"
scalar ch = 100*(1-r(estimate))
di "% change=" %5.2f ch "%"

di "Alternative - collapse three intervention phases - not reported"
`switch_quiet' nbreg n ibn.site i.site#c.time c.month_sin i.P_collapse if site<5, exposure(erp10e5) nolog noconstant irr 

* Comparison of pre-intervention phase with aggregated intervention phase
`switch_quiet' lincom _b[2.P_collapse] - _b[0.P_collapse]
di "par=" %6.3f r(estimate) ", 95% CI: " %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) ", p=" %6.4f r(p)
`switch_quiet' lincom _b[2.P_collapse] - _b[0.P_collapse], irr
di "IRR=" %6.3f r(estimate) ", 95% CI: " %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) ", p=" %6.4f r(p)
scalar ch = 100*(1-r(estimate))
di "% change=" %5.2f ch "%"

******************************************************************************************
**# TEST PRIMARY HYPOTHESIS FOR EACH SITE AND NSW
******************************************************************************************
forvalues site = 1/5 {
	di _n "************************* ${site_ID_`site'} ***********************"
estimates restore model_site`site'
`switch_quiet' nbreg, irr
* Test of 1y hypo - average of three intervention phases == pre-intervention
`switch_quiet' lincom _b[2.P]/3+_b[3.P]/3+_b[4.P]/3-_b[0.P]
di "par=" %6.3f r(estimate) ", 95% CI: " %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) ", p=" %6.4f r(p)

`switch_quiet' lincom _b[2.P]/3+_b[3.P]/3+_b[4.P]/3-_b[0.P],irr
di "IRR=" %6.3f r(estimate) ", 95% CI: " %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) ", p=" %6.4f r(p)
di  %6.3f r(estimate) "[" %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) "]"
scalar ch = 100*(1-r(estimate))
di "% change=" %5.2f ch "%"

di "Alternative - collapse three intervention phases - not reported"
`switch_quiet' nbreg n c.time c.month_sin i.P_collapse if site==`site', exposure(erp10e5) nolog difficult
* Comparison of pre-intervention phase with aggregated intervention phase
`switch_quiet' lincom _b[2.P_collapse] - _b[0.P_collapse]
di "par=" %6.3f r(estimate) ", 95% CI: " %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) ", p=" %6.4f r(p)
`switch_quiet' lincom _b[2.P_collapse] - _b[0.P_collapse], irr
di "IRR=" %6.3f r(estimate) ", 95% CI: " %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) ", p=" %6.4f r(p)
scalar ch = 100*(1-r(estimate))
di "% change=" %5.2f ch "%"
}

******************************************************************************************
* TEST PRIMARY HYPOTHESIS FOR MALES AND FEMALES 
******************************************************************************************
local sex_labels `""Male" "Female""'
forvalues sx =1/2 {
local sex_label : word `sx' of `sex_labels'
	di _n "*************************  `sex_label' ***********************"
estimates restore model_`sex_label'
quietly: nbreg, irr
* Test of 1y hypo - average of three intervention phases == pre-intervention
quietly: lincom _b[2.P]/3+_b[3.P]/3+_b[4.P]/3-_b[0.P]
di "par=" %6.3f r(estimate) ", 95% CI: " %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) ", p=" %6.4f r(p)

quietly: lincom _b[2.P]/3+_b[3.P]/3+_b[4.P]/3-_b[0.P],irr
di "IRR=" %6.3f r(estimate) ", 95% CI: " %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) ", p=" %6.4f r(p)
di  %6.3f r(estimate) "[" %6.3f  r(lb) ustrunescape("\u2013") %6.3f r(ub) "]"
scalar ch = 100*(1-r(estimate))
di "% change=" %5.2f ch "%"

}
scalar drop ch

******************************************************************************************
**# GENERATE TABLES FOR PAPER 
******************************************************************************************
esttab model_all model_Male model_Female model_site1 model_site2 model_site3  model_site4 model_site5 ///
 using "«Table 2»" ///
,b(3) not  ci(3) compress eform tab type replace constant

log close
