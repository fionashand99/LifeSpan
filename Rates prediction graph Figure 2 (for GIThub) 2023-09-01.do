********************************************************************************
* Graph of actual and predicted rates - Figure 2                               *
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


******************************************************************************************
* FINAL OVERALL MODELLING USING 8-MONTH PERIODS AND TRANSITION PERIOD
* Site specific intercepts and means, common seasonality and intervention effect
******************************************************************************************

* Detrend and set 
summ month_sin
scalar mean_month_sin = r(mean) // get average value (it's nearly 0)

* Rates per 100,000
nbreg n i.site##c.time c.month_sin i.P if site<5, exposure(erp10e5) nolog
estimates store model_site_all

recode P (0=0)  (1/4=1) (5/6=0), gen(lifespan_inter)

* Fit model to each site and NSW
forvalues site =1/5 {
di _n "************************* ${site_ID_`site'} ***********************"
nbreg n time month_sin i.P if site==`site', exposure(erp10e5) nolog
estimates store model_site_`site'
}

******************************************************************************************
**# MODEL for NSW
******************************************************************************************

* Constrain regional site parameters to equality - like site-only model
* Seasonality
constraint define 1 1.site#c.month_sin = 2.site#c.month_sin // Seasonality
constraint define 2 1.site#c.month_sin = 3.site#c.month_sin // Seasonality
constraint define 3 1.site#c.month_sin = 4.site#c.month_sin // Seasonality
* Intervention phase
constraint define 4 1.site#1.P_average = 2.site#1.P_average // Intervention - transition
constraint define 5 1.site#1.P_average = 3.site#1.P_average // Intervention - transition
constraint define 6 1.site#1.P_average = 4.site#1.P_average // Intervention - transition
*
constraint define 7 1.site#2.P_average = 2.site#2.P_average // Intervention - early
constraint define 8 1.site#2.P_average = 3.site#2.P_average // Intervention - early
constraint define 9 1.site#2.P_average = 4.site#2.P_average // Intervention - early
*
constraint define 10 1.site#3.P_average = 2.site#3.P_average // Intervention - mid
constraint define 11 1.site#3.P_average = 3.site#3.P_average // Intervention - mid
constraint define 12 1.site#3.P_average = 4.site#3.P_average // Intervention - mid
*
constraint define 13 1.site#4.P_average = 2.site#4.P_average // Intervention - late
constraint define 14 1.site#4.P_average = 3.site#4.P_average // Intervention - late
constraint define 15 1.site#4.P_average = 4.site#4.P_average // Intervention - late
*
constraint define 16 1.site#5.P_average = 2.site#5.P_average // Intervention - post 1
constraint define 17 1.site#5.P_average = 3.site#5.P_average // Intervention - post 1
constraint define 18 1.site#5.P_average = 4.site#5.P_average // Intervention - post 1
*
constraint define 19 1.site#6.P_average = 2.site#6.P_average // Intervention - post 2
constraint define 20 1.site#6.P_average = 3.site#6.P_average // Intervention - post 2
constraint define 21 1.site#6.P_average = 4.site#6.P_average // Intervention - post 

* Model includes: NSW & individual site intercept/slopes 
*                 NSW & common site seasonality (sine)
*                 NSW & common site intervention phases
*                 NSW is reference condition
nbreg n ib5.site##c.time c.month_sin ib5.site#c.month_sin ib5.site#i.P_average , ///
	constraint(1/21) exposure(erp10e5) nolog
estimates store model_site_nsw
predict rate_pred_xb_nsw , xb
gen rate_pred_deseason_nsw = exp(rate_pred_xb_nsw-_b[month_sin]*month_sin)/erp10e5

* This is to create date labels for time of the form mm/yy without 
*ssc install sencode
gen datestr=strofreal(month,"%2.0f")+"/"+strofreal(year-2000,"%2.0f")
gen ID=_n
sort year month
sencode datestr, gsort(year month) gen(date) label(date_label) // Date descriptors are now labels!
sort ID
label define date_label 0 "6/12" ,add
label value time date_label

******************************************************************************************
* PLOTS OF TRAJECTORIES AS RATES
******************************************************************************************

* Generate predicted n events
estimates restore model_site_all
capture drop rate_pred 
predict rate_pred , ir
predict rate_pred_xb , xb
predict n_pred
gen rate_pred_xb_exp = exp(rate_pred_xb)
gen rate_pred_xb_nu = exp(rate_pred_xb)/erp10e5

gen rate_pred_deseason = exp(rate_pred_xb-_b[month_sin]*month_sin)/erp10e5
gen rate_raw = n/erp10e5

gen rate_pred_deseason_cf=.
forvalues s = 1/4 {
di _n "**** `s'"	
replace rate_pred_deseason_cf=exp(_b[_cons]+_b[`s'.site]+_b[time]*time+_b[`s'.site#time]*time) if site==`s'
corr rate_pred_deseason_cf rate_pred_deseason if site==`s'&P==0
summ rate_pred_deseason_cf rate_pred_deseason if site==`s'&P==0
}

*Calc total (M+F) rates
gen rate_raw_agg = n_total/erp_total10e5  if sex==2

colorpalette cblind, cblind(.5) n(6) select(6 4 3 5 2 1) local(Site*,suffix(_col))
colorpalette cblind, cblind(.5) n(6) select(6 4 3 5 2 1) local(Site*,suffix(_int)) intensity(.4)

preserve
local active lifespan_trans_inter
drop if sex==1 // Plot only male results to avoid duplication of lines
twoway /// NB For some bizarre reason the intervention and fitted plots don't exactly coincide for all sites - fixed in Illustrator
 (line rate_raw_agg time if site==1, lw(vthin)  lcolor(`Site1_int')) /// observed
 (line rate_raw_agg time if site==2, lw(vthin)  lcolor(`Site2_int')) ///
 (line rate_raw_agg time if site==3, lw(vthin)  lcolor(`Site3_int')) ///
 (line rate_raw_agg time if site==4, lw(vthin)  lcolor(`Site4_int')) ///
 (line rate_raw_agg time if site==5, lw(vthin)  lcolor(`Site5_int')) ///
 ///
 (line rate_pred_deseason_cf time if site==1&P>0, lw(medium) lpattern(shortdash) lcolor(`Site1_col')) /// Counterfactuals
 (line rate_pred_deseason_cf time if site==2&P>0, lw(medium) lpattern(shortdash) lcolor(`Site2_col')) ///
 (line rate_pred_deseason_cf time if site==3&P>0, lw(medium) lpattern(shortdash) lcolor(`Site3_col')) ///
 (line rate_pred_deseason_cf time if site==4&P>0, lw(medium) lpattern(shortdash) lcolor(`Site4_col')) ///
 ///
 (line rate_pred_deseason time if site==1, lw(medium) lpattern(dash) lcolor(`Site1_col')) /// Fitted plots
 (line rate_pred_deseason time if site==2, lw(medium) lpattern(dash) lcolor(`Site2_col')) ///
 (line rate_pred_deseason time if site==3, lw(medium) lpattern(dash) lcolor(`Site3_col')) ///
 (line rate_pred_deseason time if site==4, lw(medium) lpattern(dash) lcolor(`Site4_col')) ///
 (line rate_pred_deseason_nsw time if site==5, lw(medium) lpattern(dash) lcolor(`Site5_col')) ///
 ///
 (line rate_pred_deseason time if site==1&lifespan_inter==1, lw(medthick) lcolor(`Site1_col')) /// 15 Thick==intervention
 (line rate_pred_deseason time if site==2&lifespan_inter==1, lw(medthick) lcolor(`Site2_col')) /// 16
 (line rate_pred_deseason time if site==4&lifespan_inter==1, lw(medthick) lcolor(`Site4_col')) /// 17
 (line rate_pred_deseason time if site==3&lifespan_inter==1, lw(medthick) lcolor(`Site3_col')) /// 18
 (scatteri . .  if site==5&lifespan_inter==1, recast(line) lw(medthick) lcolor(`Site5_col'))   /// 19 Null plot required for legend
 /// 
 (scatteri . .  if site==5&lifespan_inter==1, recast(line) lw(vthin) lcolor(`Site6_col'))  /// Generic - observations
 (scatteri . .  if site==5&lifespan_inter==1, recast(line) lw(medium) lpattern(dash) lcolor(`Site6_col'))  /// Generic - model deasoned
 (scatteri . .  if site==5&lifespan_inter==1, recast(line) lw(medthick) lcolor(`Site6_col'))  /// Generic - model intervention
 (scatteri . .  if site==5&lifespan_inter==1, recast(line) lw(medium) lpattern(shortdash) lcolor(`Site6_col'))  /// Generic - counterfactual
,xtitle("Time (month/year)" ,size(medsmall) margin(bottom)) /// X axis
 xscale(range()) ///
 xmtick(##6) ///
 xlabel(0(6)102 , labels labsize(vsmall) angle(forty_five) format(%6s) valuelabel nogrid) ///
 ///
 ytitle("Self-harm rate per 100,000 population", size(small) margin(2 2 2 2)) /// Y axis
 yscale(range(30 30)) ylabel(#6,labsize(small)) /// 
 ylabel(,nogrid) ///
 legend(order( 15 "$site_ID_1" 16 "$site_ID_2" 18 "$site_ID_3" 17 "$site_ID_4" 19 "$site_ID_5" /// legend - use last five graphs
               20 "Observations" 21 "Models (deseasoned)" 22 "Intervention phases" 23 "Counterfactuals") ///  
 pos(6) col(5) symxsize(*.75) size(vsmall) region(margin(1 1 2 0) color(white) lwidth(thick)) bmargin(zero))  /// 
 /// 
 title("") xsize(18) ysize(12) name(rate_events_ds,replace)
restore
 
graph export "«Figure 2».eps", name(rate_events_ds) replace

log close
