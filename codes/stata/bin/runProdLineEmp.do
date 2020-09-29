/* This code obtains the number of production line employees by IFR code using 
data from the BLS for 1989 */

/* Uncomment if running as stand-alone
global masterpath "/Users/andreamanera/Dropbox (MIT)/Papers Fall 18/Automation/Data/Calibration"

global out  "$masterpath/out"
global data "$masterpath/dataRaw"
*/
cd "$masterpath"
clear
* Import fixed investment deflator and save dta
import delimited "$masterpath/dataRaw/BLSProdLineEmp/mf89d3.csv"
keep sic occ_code emp occ_tit

keep if occ_code == 80000
drop if sic<2000 | sic>3999
replace sic = sic/10

save "$masterpath/dataRaw/BLSProdLineEmp.dta",replace

use "$masterpath/dataRaw/CrossWalkSICIFR.dta", clear

encode industry_ifr19, gen(ifrCode)
gen sic = floor(sic87dd/10)
bys sic ifrCode: gen ifrCount = _N
keep sic ifrCode ifrCount industry_ifr19
bys sic: egen maxCount = max(ifrCount)
duplicates drop 
* Attribute sic sector 3-digit to the IFR code containing most 4-digit sectors
bys sic: keep if ifrCount == maxCount
drop *Count

merge 1:1 sic using "$masterpath/dataRaw/BLSProdLineEmp.dta", keep(match)
drop _merge
drop occ_tit occ_code

collapse (sum) emp, by(industry_ifr19) 
replace emp = emp/1000

save "$masterpath/out/ProdLineEmpIFR.dta", replace
