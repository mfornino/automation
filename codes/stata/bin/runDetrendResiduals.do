/* This file removes a firm-specific trend from the productivity residual 
variables  contained in rawCompustatResidIFR*/


/* Uncomment if running as stand-alone
global masterpath "/Users/andreamanera/Dropbox (MIT)/Papers Fall 18/Automation/Data/Calibration"

global out  "$masterpath/out"
global data "$masterpath/dataRaw"
*/

set more off


use "$out/rawCompustatResidIFR.dta", clear


xtset gvkey fyear
sort ifrCode gvkey fyear
*generate a counter for the years each firm appears
bys gvkey: gen T_gvkey = _N

*eliminate all firms that appear less than three years (or trend kills everything)
drop if T_gvkey<3

rename tfp tfp_raw
gen p_raw = exp(tfp_raw)

* drop useless variables
keep gvkey tfp_raw p_raw fyear ifrCode sic2 industry_ifr19

* run a firm-specific time trend and generate residuals
levelsof gvkey, l(firmnames)
gen tfp = .
scalar j = 0
local count_values_num: word count of `firmnames'

foreach firm in `r(levels)'{
	scalar j = j + 1
	scalar percent = j/`count_values_num' * 100
	di percent
	capture reg tfp_raw fyear if gvkey==`firm'
	if _rc == 0{
		predict resid, residuals
		replace tfp = resid if gvkey==`firm'
		drop resid
	}
	else 
}

save "$out/detrendedCompustatResidIFR.dta", replace 
