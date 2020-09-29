# Automation and the Future of Work
Repository for my Job Market Paper "Automation and the Future of Work: Assessing the Role of Labor Flexibility" (joint with Andrea Manera).

## Abstract
We study the economic incentives for automation when labor and machines are perfect substitutes. Labor may still be employed in production, even when it is a costlier input than robots on a productivity-adjusted basis. This occurs if firms face uninsurable idiosyncratic risk, adjusting the stock of machines is costly, and workers can be hired and fired quickly enough. Even though labor survives, jobs become less stable, as workers are hired in short-lived bursts to cope with shocks. We calibrate a general equilibrium, multi-industry version of our model to match data on robot adoption in US manufacturing sectors, and use it to compute the employment and labor share consequences of progress in automation technology. A fall in the relative price of robots leads to relatively few jobs losses, while reductions in adjustment costs, or improvements in relative robot productivity, can be far more disruptive. The model-implied semi-elasticity of aggregate employment to robot penetration ranges between 0.01% and 0.1%, depending on the underlying source of increased robot adoption. Adding reduced-form hiring and firing costs to our benchmark model reveals that the scare of automation is justified when regulations impose substantial rigidity on employment relations.

[Paper](https://github.com/mfornino/automation/blob/master/paper.pdf) | [Slides](https://github.com/mfornino/automation/blob/master/slides.pdf)

## Replication
The various versions of the model are simulated in MATLAB. Please refer to the [instructions](https://github.com/mfornino/automation/blob/master/codes/matlabReadme.pdf).

The empirical analysis was carried out in Stata. Please refer to the [instructions](https://github.com/mfornino/automation/blob/master/codes/stataReadme.pdf).

The codes folder is organized as follows:

```bash
├── matlab
│   ├── generalEquilibrium.m
│   ├── GetGrids.m
│   ├── num2hms.m
│   ├── partialEquilibrium.m
│   ├── SetParametersGE.m
│   ├── SolvePolicy_rigidlabor.m
│   ├── SolveTransition.m
│   ├── Statistics.csv
│   ├── table2latex.m
│   ├── x0.csv
│   └── x0_gamma.csv
├── matlabReadme.pdf
├── stata
│   ├── bin
│   │   ├── runCompustatResidIFR.do
│   │   ├── runDetrendResiduals.do
│   │   ├── runEstimatorsIFR.do
│   │   ├── runGenerateAllData.do
│   │   ├── runMergeAll.do
│   │   ├── runProdLineEmp.do
│   │   └── runThetaProd.do
│   ├── dataRaw
│   └── out
└── stataReadme.pdf
```

## Support

Please don't hesitate to contact me should you need any information on the replication.
