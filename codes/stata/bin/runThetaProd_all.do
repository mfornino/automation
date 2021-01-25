

/* Uncomment and change if running as stand-alone
global masterpath "/Users/andreamanera/Dropbox (MIT)/Papers Fall 18/Automation/Data/oldcalibration"
*/


use "$masterpath/dataRaw/sic5811.dta", clear



/*
preserve
clear 
import delimited "$masterpath/dataRaw/CPIAUCSL.csv"
rename date year
save "$masterpath/dataRaw/CPIAUCSL.dta",replace
clear 
restore
merge m:1 year using "$masterpath/dataRaw/CPIAUCSL.dta"
drop _merge
rename cpiaucsl_nbd19870101 cpiu
replace cpiu = cpiu/100
*/

/* Creates industry residuals by year */
set matsize 800
sort sic year
xtset sic year
destring sic, replace

gen sic2 = floor(sic/100)
drop if sic2<20 | sic2>39


gen share_vship = prodw/vship
gen share_vadd = prodw/vadd
* gen share_labor_vship = pay/cpiu/vship


label define siclabel 20 "Food & Kindred Products" ///
	21 "Tobacco Products" ///
	22 "Textile Mill Products" ///
	23 "Apparel & Other Textile Products" ///
	24 "Lumber & Wood Products" ///
	25 "Furniture & Fixtures" ///
	26 "Paper & Allied Products" ///
	27 "Printing & Publishing" ///
	28 "Chemical & Allied Products" ///
	29 "Petroleum & Coal Products" ///
	30 "Rubber & Miscellaneous Plastics Products" ///
	31 "Leather & Leather Products" ///
	32 "Stone, Clay, & Glass Products" ///
	33 "Primary Metal Industries" ///
	34 "Fabricated Metal Products" ///
	35 "Industrial Machinery & Equipment" ///
	36 "Electronic & Other Electric Equipment" ///
	37 "Transportation Equipment" ///
	38 "Instruments & Related Products" ///
	39 "Miscellaneous Manufacturing Industries"
	
label define siclabelshort 20 "Food" ///
	21 "Tobacco" ///
	22 "Textile Mill" ///
	23 "Apparel & Other Textile" ///
	24 "Lumber & Wood" ///
	25 "Furniture & Fixtures" ///
	26 "Paper" ///
	27 "Printing & Publishing" ///
	28 "Chemical" ///
	29 "Petroleum & Coal" ///
	30 "Rubber & Plastics" ///
	31 "Leather" ///
	32 "Stone, Clay, & Glass" ///
	33 "Primary Metal" ///
	34 "Fabricated Metal" ///
	35 "Industrial Machinery" ///
	36 "Electronic & Other Electric" ///
	37 "Transportation Equipment" ///
	38 "Instruments" ///
	39 "Miscellaneous"

label values sic2 siclabelshort






rename sic sic87dd

* Merge with Crosswalk IFR

merge m:1 sic87dd using "$masterpath/dataRaw/CrossWalkSICIFR.dta", keep(match)
drop _merge
drop notes
encode sic87dd_desc, gen(sicDesc)


* collapse by IFR sector
collapse  piship (sum) vship (sum) pay (sum) vadd (sum) prodw, by(year) 

gen share_sums_vadd_prodw = prodw/vadd


tsset year

* HP-filter shares to remove fluctuations
tsfilter hp share_prodw_c = share_sums_vadd_prodw, smooth(6.25) trend(share_prodw_t)

* Tale select shares for specific years for collapse
// bysort ifrCode: egen vadd_1989 = total(cond(year == 1989, vadd, .))
// bysort ifrCode: egen share_prodw_t_04 = total(cond(year == 2004, share_prodw_t, .))
// bysort ifrCode: egen share_prodw_t_10 = total(cond(year == 2010, share_prodw_t, .))

// gen d_share_prodw_t = share_prodw_t_10 - share_prodw_t_04


* Collapse by IFR by year
// collapse vadd_1989  d_share* share_prodw_t*, by(industry_ifr19 year)

// encode industry_ifr19, gen(ifrCode)

/* Uncomment for graph
sepscatter share_prodw_t  year , separate(ifrCode) recast(connect) name(seriesThetaProd) title("Production Worker Shares")
*/


* Compute average shares of production-line employees in each sector before 1980
collapse  thetaProd = share_prodw_t if year<1980 

save "$masterpath/out/ThetaProdIFR_all.dta", replace


