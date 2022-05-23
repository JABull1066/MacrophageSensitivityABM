# MacrophageSensitivityABM
Code associated with "Macrophage sensitivity to microenvironmental cues influences spatial heterogeneity of tumours" (Bull and Byrne, 2022)


Running this code requires installation of the open source Chaste (Cancer, Heart and Soft Tissue Environment) software, available at https://www.chaste.cs.ox.ac.uk/. This repository contains additional code required to run the main simulation executable, Exe_2022_RunSimulation.cpp. This executable can be called with the following command line arguments, which alter the parameters of the simulation:


- "ID" - ID string for the simulation
- "ACCD_T" - averageCellCycleDuration_Tumour: Average cell cycle duration for tumour cells ($\tau_{tum}$)
- "ACCD_S" - averageCellCycleDuration_Stroma: Average cell cycle duration for stromal cells ($\tau_{str}$)
- "CS_M2CSF" - chemotacticSensitivity_macrophageToCSF: chi value for macrophage chemotaxis sensitivity to CSF ($\chi_{c}^{m}$)
- "CS_M2CXCL12" - chemotacticSensitivity_macrophageToCXCL12: chi value for macrophage chemotaxis sensitivity to CXCL12 ($\chi_{\xi}^{m}$)
- "CS_TC2EGF" - chemotacticSensitivity_tumourCellToEGF: chi value for tumour cell chemotaxis sensitivity to EGF ($\chi_{\epsilon}^{T}$)
- "MPI" - macrophagePhenotypeIncrementPerHour: amount to increment macrophage phenotype by each hour when above TGF threshold ($\Delta P$)
- "DTBV" - distanceToBloodVessels: initial distance from tumour cells to blood vessels ($R_{B}$)
- "NOBV" - numberOfBloodVessels: target number of blood vessels to include in the simulation ($N_{B}$)
- "HMECSF" - halfMaximalExtravasationCsf1Conc: concentration of CSF1 at which the extravasation rate will be half maximal ($c_{1/2}$)
- "MPOE" - maximumProbOfExtravasationPerHour: maximum rate of extravasation from a blood vessel ($P^{\star}$)
- "TGF_TFPS" - tgfThresholdForPhenotypeSwitch: TGF threshold above which macrophage will begin to polarise ($g_{crit}$)
- "IT" - Iteration Number
