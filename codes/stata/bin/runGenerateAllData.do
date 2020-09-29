/* This File generates the input data for the Matlab calibration routine 
contained in "Statistics.csv" */

/* Choose a masterpath containing the subfolders "dataRaw", "bin", "out */

global masterpath "/Users/andreamanera/Dropbox (MIT)/Papers Fall 18/Automation/⁨submission⁩/codes/stata"

global out  "$masterpath/out"
global data "$masterpath/dataRaw"
global bin "$masterpath/bin"

* Obtain the DRS parameters as prod.-line employees wage bill over VA
do "$bin/runThetaProd.do"

* Obtain prod.-line employees in each IFR sector
do "$bin/runProdLineEmp.do"

* Generate sales residulas from Compustat
do "$bin/runCompustatResidIFR.do"

* Detrend Compustat Residuals
do "$bin/runDetrendResiduals.do"

* Obtain estimates of the standardized CIR process parameters
do "$bin/runEstimatorsIFR.do"

* Obtain estimates of the standardized CIR process parameters
do "$bin/runMergeAll.do"
