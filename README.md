# mCSPC-NMA-publication
Code and data for mCSPC NMA publication


# Instructions on the fractional polynomial models
The fractional polynomials (PF) folder contains 2 subfolders: data and code.
The data folder contains the treatment codes (mHSPC OS trt codes.xlsx and mHSPC PFS trt codes.xlsx), which are only used for network plotting and to confirm that the correct treatment is used in the analyses. 
The FP_data_OS.csv and FP_data_PFS.csv files contain reconstructed individual patient data (IPD) from Kaplan-Meier (KM) curves, ready to be analysed. Each file is a data frame with 6 columns (s, a, dt r, z, time_j). s correspond to the study indicator, from 1 to 13. a is the correspondent arm in each study. dt is the time interval (in months) where a Piecewise Constant Hazard as approximate Likelihood is implemented. r is the number of events recorded during the time interval and z is the patients at risk at the beginning of the interval. time_j records the time at which the time interval ends (also in months).
The code folder implements the FP models on both PFS and OS outcomes. The OpenBUGS models are in the txt files, including fixed- and random-effects, first- and second-order FP models.
The OS_FP_models.R and the PFS_FP_models.R implement all the FP models explored during the analyses. OS_model_averaging.R and PFS_model_averaging.R implement model averaging from the 5 models with the lowest deviation information criteria (DIC) from the previous analyses. Finally, the OS_survival_curves.R and the PFS_survival_curves.R scripts implement survival extrapolation per treatment.