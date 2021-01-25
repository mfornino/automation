
/* Uncomment if running as stand-alone
global masterpath "/Users/andreamanera/Dropbox (MIT)/Papers Fall 18/Automation/Data/Calibration"

global out  "$masterpath/out"
global data "$masterpath/dataRaw"
*/


cd "$masterpath"


* Import Compustat query
use "$masterpath/dataRaw/Compustatquery0726.dta", clear
destring gvkey, replace
xtset gvkey fyear
set matsize 11000
set more off

preserve
clear 
* Import fixed investment deflator and save dta
import delimited "$masterpath/dataRaw/FIXEDINVDEF.csv"
rename date fyear
save "$masterpath/dataRaw/FIXEDINVDEF.dta",replace
clear 
restore
merge m:1 fyear using "$masterpath/dataRaw/FIXEDINVDEF.dta"
destring sic, replace
drop _merge
drop if sic<2000 | sic>3999
drop if at==.|emp==.|sale==.|ppegt==.|ppent==.


/* generate the first value of gross capital stock */
bys gvkey : gen countnonmissing = sum(!missing(ppegt)) if !missing(ppegt)
bysort gvkey (countnonmissing) : gen ppegt_first = ppegt[1]/fixedinvdef[1]
bysort gvkey (countnonmissing) : gen ppegt_year_first = fyear[1]
bysort gvkey: drop if fyear<ppegt_year_first

/* interpolate missing net capital stock*/
bysort gvkey: ipolate ppent fyear, generate(ippent)
label var ippent "Interpolated net capital stock"

/* generate net investment */
sort gvkey fyear
bys gvkey: gen realnetinv =(ippent - l.ippent)/fixedinvdef
bys gvkey: gen kstock = ppegt_first + sum(realnetinv)
replace kstock = ppegt_first if kstock ==.

/* further cleaning */
drop if sale <= 0 | emp <=0
/*drop if aqc/fixedinvdef > .05*kstock*/

replace sale = sale/gdpdef



gen l_sale = log(sale)
gen l_emp = log(emp)
gen l_lkstock = log(l.kstock)
* Include lags for inventories
gen lag_l_emp = l.l_emp
gen lag_l_lkstock = l.l_lkstock

gen sic2 = floor(sic/100)
rename sic sic87

save "$masterpath/out/dataCleanedCompustat.dta", replace

* merge crosswalk
use "$masterpath/dataRaw/CrossWalkSICIFR.dta", clear

drop notes
rename sic87dd sic87
merge 1:m sic87 using "$masterpath/out/dataCleanedCompustat.dta", keep(match)
drop _merge

save "$masterpath/out/dataCleanedCompustat.dta", replace

use "$masterpath/dataRaw/apr_measures_ifr19.dta", clear
keep industry* apr_us_lv*
merge 1:m industry_ifr19 using "$masterpath/out/dataCleanedCompustat.dta"
drop _merge

encode industry_ifr19, gen(ifrCode)

* Obtain TFP residuals by IFR sector
xtreg l_sale i.fyear#i.ifrCode i.ifrCode#c.l_emp i.ifrCode#c.l_lkstock, fe  
*i.ifrCode#c.lag_l_emp  i.ifrCode#c.lag_l_lkstock 

predict tfp, e



save "$masterpath/out/rawCompustatResidIFR.dta", replace
