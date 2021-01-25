/* This File generates the input data for the Matlab calibration routine 
contained in "Statistics.csv" */

/* Choose a masterpath containing the subfolders "dataRaw", "bin", "out */


global masterpath "/Users/andreamanera/Dropbox (MIT)/Papers Fall 18/Automation/Data/Calibration"

global out  "$masterpath/out"
global data "$masterpath/dataRaw"
global bin "$masterpath/bin"

* Obtain the DRS parameters as prod.-line employees wage bill over VA
do "$bin/runThetaProd.do"

* Obtain the DRS parameters for the representative sector
do "$bin/runThetaProd_all.do"

* Obtain prod.-line employees in each IFR sector
do "$bin/runProdLineEmp.do"

* Generate sales residuals from Compustat
do "$bin/runCompustatResidIFR.do"

* Outsheet Raw Residuals for Non-stationary Calibration
do "$bin/outsheetResiduals.do"

* Obtain Stationary Residuals Using the Hamilton Filter
do "$bin/runHamiltonResiduals.do"

* Obtain estimates of the standardized OU parameters
do "$bin/runOUEstimatorsIFR.do"

* Obtain estimates of the standardized OU parameters for the representative sector
do "$bin/runOUEstimatorsIFR_allsectors.do"

* Merge all datasets and produce csv for Matlab
do "$bin/runMergeAll.do"
