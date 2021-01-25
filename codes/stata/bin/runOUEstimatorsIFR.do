/* This file produces the estimators for the OU parameters as used in the paper.
It assumes that firms are sampled i.i.d. (all firms are lumped together) */


/* Uncomment if running as stand-alone
global masterpath "/Users/andreamanera/Dropbox (MIT)/Papers Fall 18/Automation/Data/Calibration"

global out  "$masterpath/out"
global data "$masterpath/dataRaw"
*/



set more off

local dataType "Hamilton"

foreach dT of local dataType {

 
* use detrended or non-detrended residuals
use "$out/`dT'CompustatResidIFR.dta", clear

* restrict to firms where we can reliably reject ADF test
keep if N_firm_df>20 & adf_p_val<.1

preserve
keep ifrCode industry_ifr19 
duplicates drop ifrCode industry_ifr19 , force
save "$out/`dT'XWalkIFRCodeIndustry.dta", replace
restore

if ("`dT'"=="HPF"|"`dT'"=="Hamilton"){
	gen p = tfp_cycle
}
else{
	gen p = tfp
}

/*plot the data*/
sort  ifrCode gvkey fyear 

* Winsorize at desired percentile

winsor2 p, replace cuts(1 99) by(ifrCode)

outsheet using "$out/`dT'EstimationSample.csv", replace
*sepscatter p fyear if gvkey<2000, separate(gvkey) recast(connect) legend(off) name(trends)
* Follows Tang and Chen (2009)
xtset gvkey fyear

gen lag_p = l.p
gen lag_sq = lag_p^2
gen p_times_lag = p*lag_p


* drop firms with missing value of the lag
keep if  lag_p!=.

bys gvkey: gen n_gvkey = _N

preserve


collapse (sum) sum_p = p ///
  sum_lag_p = lag_p ///
  sum_p_times_lag = p_times_lag ///
  sum_lag_sq = lag_sq ///
  (max) n_gvkey ifrCode , by(gvkey)


/*extract number of observations*/
egen n_sample_obs = max(_n)



bys ifrCode: egen N_ifrCode = total(n_gvkey)

collapse (sum) sum_p sum_lag_p sum_p_times_lag sum_lag_sq ///
  (max) N_ifrCode , by(ifrCode)
  
  
gen beta_1 = ( N_ifrCode^(-1)*sum_p_times_lag -  N_ifrCode^(-2)*sum_p*sum_lag_p )/( N_ifrCode^(-1)*sum_lag_sq - N_ifrCode^(-2)*sum_lag_p^2 )
  

gen beta_2 = N_ifrCode^(-1)*( sum_p - beta_1*sum_lag_p)/(1 - beta_1)

  
  
save "$out/`dT'OUBetaIFR.dta", replace

restore

sort ifrCode
merge m:1 ifrCode using "$out/`dT'OUBetaIFR.dta", keep(match)

gen resid_3 = (p - beta_1*lag_p - beta_2 * (1- beta_1))^2

collapse (sum) resid_3  (max)  beta_1 (max) beta_2 (max) N_ifrCode,  by(ifrCode)

gen beta_3 = N_ifrCode^(-1) * resid_3


* Vasicek Parameters
gen theta_log = -log(beta_1)
rename beta_2 mu_log
gen sigma_log = sqrt(2 * theta_log * beta_3 /(1 - beta_1^2)) 

* Obtain parameters of the exponential Ornstein-Uhlenbeck
gen mu_p = mu_log - sigma_log^2/( 2 * theta_log)
rename sigma_log sigma
rename theta_log theta_p
gen sigma_p = sqrt(sigma^2/(2*theta_p))


merge 1:1 ifrCode using "$out/`dT'XWalkIFRCodeIndustry.dta", keep(match)

keep sigma_p theta_p industry_ifr19


if "`dT'" == "detrended" {
	rename sigma_p sigma_p_detrend
	rename theta_p theta_p_detrend
}	

if "`dT'" == "HPF" {
	rename sigma_p sigma_p_HPF
	rename theta_p theta_p_HPF
}	


save "$out/`dT'OUEstimatorsIFR.dta", replace

}
