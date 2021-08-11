# Automation and the Future of Work
Repository for my paper
"Automation and the Future of Work: Assessing the Role of Labor Flexibility" |
Review of Economic Dynamics | joint with Andrea Manera

## Abstract
We study the economic incentives for automation when labor and machines are perfect substitutes. Labor may still be employed in production, even when it is a costlier input than robots on a productivity-adjusted basis. This occurs if firms face uninsurable idiosyncratic risk, adjusting the stock of machines is costly, and workers can be hired and fired quickly enough. Even though labor survives, jobs become less stable, as workers are hired in short-lived bursts to cope with shocks. We calibrate a general equilibrium, multi-industry version of our model to match data on robot adoption in US manufacturing sectors, and use it to compute the employment and labor share consequences of progress in automation technology. A fall in the relative price of robots leads to relatively few jobs losses, while reductions in adjustment costs, or improvements in relative robot productivity, can be far more disruptive. The model-implied semi-elasticity of aggregate employment to robot penetration (number of robots per thousand employees) ranges between 0.01% and 0.12%, depending on the underlying source of increased robot adoption, consistent with findings in the empirical literature. In an extension, we show that reduced-form hiring and firing costs unambiguously depress long-run employment.

[Paper](https://github.com/mfornino/automation/blob/master/paper.pdf) | [Slides](https://github.com/mfornino/automation/blob/master/slides.pdf)

## Replication
The starting point for the replication is to obtain the empirical moments that are used to calibrate the model. The user should start from running the Stata code in runGenerateAllData.do, which is a wrapper for all the routines that both assemble intermediate datasets, and run empirical tests. For the convenience of the reader, we also include intermediate datasets that significantly speed up the execution (but no Compustat data). Please refer to the PDF [instructions](https://github.com/mfornino/automation/blob/master/codes/stataReadme.pdf) file for additional details and for an explanation on how to obtain a suitable Compustat extract.

The simulation and calibration of the model is performed in MATLAB. The user should start from the file Run.m, with more details provided in the PDF [instructions](https://github.com/mfornino/automation/blob/master/codes/matlabReadme.pdf) file. Random numbers are only used to generate Figure 3, and the seed is reported in the Run.m file at the beginning of the specific section. We ran these numerical exercises on a workstation with Linux kernel 5.9. The CPU is a Ryzen Threadripper 2950X with 16 cores and the RAM is 64 GB of DDR4-3200. The expected completion time with these resources is about 12 hours for the Matlab simulation and 2 hours for the Stata exercise.

The codes folder is organized as follows:

```bash
├── matlab
│   ├── data
│   │   ├── emp_BEA_IFR.csv
│   │   ├── GBMStatistics.csv
│   │   ├── HamiltonEstimationSample.csv
│   │   ├── OUStatistics.csv
│   │   └── RawResiduals.csv
│   ├── graphs
│   ├── guesses
│   │   ├── Gamma.csv
│   │   ├── wpath.csv
│   │   └── x0.csv
│   ├── tables
│   ├── computeLoss.m
│   ├── DiscretizeDiffusion.m
│   ├── DiscretizeGBM.m
│   ├── EOU_cs_workspace.mat
│   ├── GBM_cs_workspace.mat
│   ├── generalEquilibrium.m
│   ├── generalEquilibrium_onesector.m
│   ├── GetGrids.m
│   ├── GetGrids_rigidlabor.m
│   ├── ind2sub_vec.m
│   ├── LaborDemand_linear.m
│   ├── LaborDemand_rigidlabor.m
│   ├── LaborDemand_trapz.m
│   ├── num2hms.m
│   ├── partialEquilibrium.m
│   ├── RunEOUCalibration.m
│   ├── RunGBMCalibration.m
│   ├── Run.m
│   ├── SetParametersGE.m
│   ├── SetParameters.m
│   ├── solvePathpR.m
│   ├── SolvePolicy_rigidlabor.m
│   ├── SolvePolicy_trapz.m
│   ├── SolveTransition.m
│   ├── SolveTransition_rigidlabor.m
│   ├── SolveTransitionWrapper.m
│   ├── sub2ind_vec.m
│   └── table2latex.m
├── stata
│   ├── bin
│   │   ├── outsheetResiduals.do
│   │   ├── runCompustatResidIFR.do
│   │   ├── runGenerateAllData.do
│   │   ├── runHamiltonResiduals.do
│   │   ├── runMergeAll.do
│   │   ├── runOUEstimatorsIFR_allsectors.do
│   │   ├── runOUEstimatorsIFR.do
│   │   ├── runProdLineEmp.do
│   │   ├── runThetaProd_all.do
│   │   └── runThetaProd.do
│   ├── dataRaw
│   │   ├── BLSmeanwages
│   │   │   ├── oesm04nat
│   │   │   │   ├── field_descriptions.xls
│   │   │   │   └── national_may2004_dl.xls
│   │   │   ├── oesm05nat
│   │   │   │   ├── field_descriptions.xls
│   │   │   │   └── national_may2005_dl.xls
│   │   │   ├── oesm07nat
│   │   │   │   ├── field_descriptions.xls
│   │   │   │   └── national_May2007_dl.xls
│   │   │   ├── oesm10nat
│   │   │   │   ├── field_descriptions.xls
│   │   │   │   └── national_M2010_dl.xls
│   │   │   ├── oesm14nat
│   │   │   │   ├── field_descriptions.xlsx
│   │   │   │   └── national_M2014_dl.xlsx
│   │   │   ├── oesm17nat
│   │   │   │   ├── field_descriptions.xlsx
│   │   │   │   └── national_M2017_dl.xlsx
│   │   │   ├── mf95d2.xls
│   │   │   └── national_1997_dl.xls
│   │   ├── BLSProdLineEmp
│   │   │   └── mf89d3.csv
│   │   ├── apr_measures_ifr19.dta
│   │   ├── AR19JPE_APR_level04.csv
│   │   ├── AR19JPE_APR_level04.dta
│   │   ├── BLSProdLineEmp.dta
│   │   ├── CPIAUCSL.csv
│   │   ├── CPIAUCSL.dta
│   │   ├── CrossWalkSICIFR.dta
│   │   ├── FIXEDINVDEF.csv
│   │   ├── FIXEDINVDEF.dta
│   │   ├── sic5811.dta
│   │   ├── VA_DATA14.csv
│   │   └── VALUEADDEDBEA.xls
│   └── out
│       ├── dataCleanedCompustat.dta
│       ├── HamiltonCompustatResidIFR.dta
│       ├── HamiltonEstimationSample.csv
│       ├── HamiltonOUBetaIFR.dta
│       ├── HamiltonOUBetaIFR_total.dta
│       ├── HamiltonOUEstimatorsIFR.dta
│       ├── HamiltonOUEstimatorsIFR_total.dta
│       ├── HamiltonXWalkIFRCodeIndustry.dta
│       ├── OUStatistics.csv
│       ├── ProdLineEmpIFR.dta
│       ├── rawCompustatResidIFR.dta
│       ├── RawResiduals.csv
│       ├── ThetaProdIFR_all.dta
│       ├── ThetaProdIFR.dta
│       └── ThetaProdTsIfr.dta
├── matlabReadme.pdf
└── stataReadme.pdf
```

## Support

Please don't hesitate to contact me should you need any information on the replication.
