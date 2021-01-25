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

drop ifrCode
egen ifrCode = group(industry_ifr19)

sort ifrCode gvkey fyear
*generate a counter for the years each firm appears
bys gvkey: gen T_gvkey = _N

*eliminate all firms that appear less than three years (or trend kills everything)

rename tfp tfp_raw

* drop useless variables
keep gvkey tfp_raw  fyear ifrCode sic2 industry_ifr19


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

drop tfp_raw max_id_spell _spell _end _seq sic2 N_firm max_spell
rename max_max_spell N_firm

outsheet * using "$out/RawResiduals.csv", comma replace


