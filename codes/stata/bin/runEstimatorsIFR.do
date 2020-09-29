
/* Uncomment if running as stand-alone
global masterpath "/Users/andreamanera/Dropbox (MIT)/Papers Fall 18/Automation/Data/Calibration"

global out  "$masterpath/out"
global data "$masterpath/dataRaw"
*/


set more off

local dataType "raw detrended"

foreach dT of local dataType {

 
* use detrended or non-detrended residuals
use "$out/`dT'CompustatResidIFR.dta", clear

gen p = exp(tfp)
/*plot the data*/
sort  ifrCode gvkey fyear 

* Winsorize at desired percentile

winsor2 p, replace cuts(1 99) by(ifrCode)
*sepscatter p fyear if gvkey<2000, separate(gvkey) recast(connect) legend(off) name(trends)
xtset gvkey fyear
gen lag_p = l.p
gen inv_l_p = 1/lag_p
gen growth_p = p/l.p
gen diff_p_firm = p-lag_p

* drop firms with missing value of the shifter parameter
keep if  growth_p!=.

bys gvkey: gen n_gvkey = _N

preserve

collapse (sum) sum_p_firm = p  sum_g_p_firm = growth_p   diff_p_firm   sum_inv_l_p_firm = inv_l_p  sum_l_p_firm = lag_p (max) n_gvkey ifrCode , by(gvkey)


/*extract number of observations*/
egen n_sample_obs = max(_n)



bys ifrCode: egen N_ifrCode = total(n_gvkey)

collapse (sum) ssum_p = sum_p_firm ssum_g_p = sum_g_p_firm ssum_diff_p = diff_p_firm  ssum_inv_l_p = sum_inv_l_p_firm ssum_l_p = sum_l_p_firm (max) N_ifrCode , by(ifrCode)

* estimate the diffusion coefficients by max likelihood
gen alph = ( N_ifrCode *ssum_p - ssum_g_p*ssum_l_p)/( N_ifrCode ^2-ssum_l_p*ssum_inv_l_p) 
gen beta = ( N_ifrCode ^2 -  N_ifrCode *ssum_g_p + ssum_diff_p*ssum_inv_l_p)/( N_ifrCode ^2-ssum_l_p*ssum_inv_l_p)


save "$out/`dT'AlphaBetaIFR.dta", replace

restore

sort ifrCode
merge m:1 ifrCode using "$out/`dT'AlphaBetaIFR.dta", keep(match)

* generate the diffusion residuals
sort gvkey fyear	
bys gvkey : gen z = (p - l.p - (alph - beta * l.p))/sqrt(l.p)
bys gvkey : gen z_2 = z^2


collapse (mean) beta alph z_2, by(industry_ifr19)

* obtain diffusion mean
gen p_d = alph/beta


* create standardized version of the process so that p_d ==1
gen theta_p = beta
gen sigma = sqrt(z_2)/sqrt(p_d)
gen sigma_p = sqrt(sigma^2/(2*theta_p))
/*
sepscatter sigma_p ifrCode, separate(ifrCode) name(`dT'sigma) title("sigma")
sepscatter theta_p ifrCode, separate(ifrCode) name(`dT'theta) title("theta")
*/
keep sigma_p theta_p industry_ifr19

if "`dT'" == "detrended" {
	rename sigma_p sigma_p_detrend
	rename theta_p theta_p_detrend
}	


save "$out/`dT'CIREstimatorsIFR.dta", replace

}
