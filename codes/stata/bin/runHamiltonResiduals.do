/* This file removes a firm-specific trend from the productivity residual 
variables  contained in rawCompustatResidIFR*/


* Uncomment if running as stand-alone
global masterpath "/Users/andreamanera/Dropbox (MIT)/Papers Fall 18/Automation/Data/Calibration"

global out  "$masterpath/out"
global data "$masterpath/dataRaw"
*

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

gen tfp = .
scalar j = 0
gen adf_p_val = .
gen adf_p_val_trend = .
gen N_firm_df = .
gen trend_stationary = .
gen stationary = .
xtset gvkey fyear
bys gvkey: gen N_firm = _N
*keep if N_firm>50
bysort gvkey: ipolate tfp_raw fyear, generate(tfp_ip)
drop if tfp_ip==.

* compute the spells for various gvkey observations
tsspell, f(L.fyear==.)
* select observations to avoid using different spells with gaps
* Find length of spells
bys gvkey _spell: egen max_spell= max(_seq)
* Find maximum length of spells by gvkey and keep only longest spells
bys gvkey: egen max_max_spell= max(max_spell)
* find highest progressive spell id in remaining series and keep most recent
bys gvkey: egen max_id_spell= max(_spell)
keep if _spell==max_id_spell

xtset gvkey fyear
* Residualize using hamilton on panel data
*gen tfp_cycle = tfp_ip - L2.tfp_ip
gen tfp_cycle = .

levelsof gvkey, l(firmnames)
local count_values_num: word count of `firmnames'

levelsof gvkey, l(firmnames)
foreach firm in `r(levels)'{
*foreach firm in `firm'{
	scalar j = j + 1
	scalar percent = j/`count_values_num' * 100
	di percent
	cap reg tfp_ip L2.tfp_ip L3.tfp_ip L4.tfp_ip L5.tfp_ip if gvkey==`firm'
	cap drop tfp_resid
	cap predict tfp_resid, resid
	cap replace tfp_cycle = tfp_resid if gvkey==`firm'
	capture dfuller tfp_cycle if gvkey==`firm', lags(1) trend 
	if _rc == 0{
		qui replace adf_p_val_trend = r(p) if gvkey==`firm'
		qui replace N_firm_df = r(N) if gvkey==`firm'
		qui replace trend_stationary = (adf_p_val_trend <.05) if gvkey==`firm' & adf_p_val_trend!=.
	}
	capture dfuller tfp_cycle if gvkey==`firm'
	if _rc == 0{
		qui replace adf_p_val= r(p) if gvkey==`firm'
		qui replace stationary = (adf_p_val <.05) if gvkey==`firm' & adf_p_val!=.
	}
else 
}

drop if N_firm_df==.

sum adf_p_val if N_firm_df>20, detail
*sum adf_p_val if N_firm>50, detail
*collapse  adf_p_val adf_p_val_trend N_firm_df, by(gvkey)

*keep if N_firm>25

save "$out/HamiltonCompustatResidIFR.dta", replace 

