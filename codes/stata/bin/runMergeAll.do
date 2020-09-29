


/* Uncomment if running as stand-alone
global masterpath "/Users/andreamanera/Dropbox (MIT)/Papers Fall 18/Automation/Data/Calibration"

global out  "$masterpath/out"
global data "$masterpath/dataRaw"
*/


* Import raw value-added data 
import delimited "$data/VA_DATA14.csv", varn(1) clear

* merge with the data from AR '19 replication kit
merge 1:1 industry_ifr19 using "$data/apr_measures_ifr19.dta"
drop _merge

* encode the industry_ifr variable
encode industry_ifr19, gen(ifrCode)

* Merge with the 2004 levels from table A.1 in AR '19 to transform changes into levels 
merge 1:1 ifrCode using "$data/AR19JPE_APR_level04", keep(match)
drop _merge

rename apr_lv04 apr_us_lv04
keep apr_us_lv* emp_us_ ifrCode industry_ifr19 shares

* generated the levels of unadjusted apr
gen apr_us_lv07 = apr_us_lv04 + apr_us_lv04_07_
gen apr_us_lv10 = apr_us_lv04 + apr_us_lv04_10_
gen apr_us_lv14 = apr_us_lv04 + apr_us_lv04_14_

rename apr_us_lv04 apr_lv04 
rename apr_us_lv07 apr_lv07
rename apr_us_lv10 apr_lv10
rename apr_us_lv14 apr_lv14

keep apr_lv* emp_us_ ifrCode industry_ifr19 shares

* Merge with the DRS parameters Theta
merge 1:1 industry_ifr19 using "$out/ThetaProdIFR.dta", keep(match)
drop _merge

* Merge with estimated stochastic process parameters
merge 1:1 industry_ifr19 using "$out/rawCIREstimatorsIFR.dta", keep(match)
drop _merge

* Merge with estimated stochastic process parameters from detrended series (actually used in Matlab)
merge 1:1 industry_ifr19 using "$out/detrendedCIREstimatorsIFR.dta", keep(match)
drop _merge

* Merge with labor share of production line employees
merge 1:1 industry_ifr19 using "$out/ProdLineEmpIFR.dta", keep(match)
drop _merge

drop emp_us_
rename emp emp_us_

order industry_ifr19 apr_lv04  apr_lv07 apr_lv10 apr_lv14 thetaProd theta_p sigma_p theta_p_detrend  sigma_p_detrend emp_us_ vadd_1989 shares

gen ifrString = "Automotive"
replace ifrString = "Electronics" if industry_ifr19 == "electronics"
replace ifrString = "Food and Beverages" if industry_ifr19 == "food"
replace ifrString = "Wood and Furniture" if industry_ifr19 == "furniture"
replace ifrString = "Miscellaneous Manufacturing" if industry_ifr19 == "manufacturing_other"
replace ifrString = "Basic Metals" if industry_ifr19 == "metal_basic"
replace ifrString = "Industrial Machinery" if industry_ifr19 == "metal_machinery"
replace ifrString = "Metal Products" if industry_ifr19 == "metal_products"
replace ifrString = "Clay Glass and Minerals" if industry_ifr19 == "mineral"
replace ifrString = "Paper and Publishing" if industry_ifr19 == "paper"
replace ifrString = "Plastics Chemicals and Pharmaceuticals" if industry_ifr19 == "petrochemicals"
replace ifrString = "Apparel and Textiles" if industry_ifr19 == "textiles"
replace ifrString = "Shipbuilding and Aerospace" if industry_ifr19 == "vehicles_other"

order industry_ifr19 ifrString apr_lv04  apr_lv07 apr_lv10 apr_lv14 thetaProd theta_p sigma_p theta_p_detrend  sigma_p_detrend emp_us_ vadd_1989 shares

label var apr_lv04 "Robots per thousand employees (2004)"
label var apr_lv07 "Robots per thousand employees (2007)"
label var apr_lv10 "Robots per thousand employees (2010)"
label var apr_lv14 "Robots per thousand employees (2014)"
label var thetaProd "DRS parameter theta"
label var theta_p "CIR parameter theta (raw firm residuals)"
label var sigma_p "Revenue-shock standard deviation (raw firm residuals)"
label var theta_p_detrend "CIR parameter theta (detrended firm residuals)"
label var sigma_p_detrend "Revenue-shock standard deviation (detrended firm residuals)"
label var emp_us_ "Sector Production-Line Employees in 1989 (thousands)"
label var vadd_1989 "Sector Value Added in 1989 (millions USD)"

outsheet * using "$out/Statistics.csv", comma replace


gen ballparkW = thetaProd/(emp* 1000/vadd_1989 / 1000000)
