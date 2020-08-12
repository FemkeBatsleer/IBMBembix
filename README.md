# IBMBembix
source code and data for manuscript

## Microhabitat model
- 'Microhabitat Final Model with iterations - scales.R': source code for the microhabitat model with INLA, including 10 iterations for cross-validation, and choice of scale.
- 'Microhabitat Final Model visualisations of iterations.R': source code for making graphs of the microhabitat model. Datafiles needed in this script are made by 'Microhabitat Final Model with iterations.R'.
- MicrohabitatModel_FullData_04062020.xlsx: raw datafile to start microhabitat model.
- RegPoints_04062020.txt: raw datafile of the regular points to make predictions of the model, the habitat suitability map.

Raw GIS data (rasters: raw data - CIR (band 1 NIR, band 3 R) DEM and derived - NDVI, insolation, slope) can be found on https://drive.google.com/drive/folders/12TNHrsqpTIvBvV4oNfkj-UeUmv68ufLx?usp=sharing

## IBM
- Bembix_model_simple.py: actual IBM, used by all Bembix_*_s.py
- Bembix_*_s.py: loops for all scenarios to make several runs for a scenario (input: begin number, end number e.g. 0 1000: makes 1000 runs)
- Hist_distances.txt, suitabilityINLAv2.jpg: input for Bembix_model_simple.py
- directory 'analyse runs': two scripts to extract the parameters and calculate summary statistics: spatial pattern statistics and network analyses from the runs of the IBM; H_* scripts are helper scripts used in the two main scripts. The map 'Outputs' gives the outputs for these of the IBM runs made on the HPC@UGent. These are used in the actual ABC analysis.

Note: some names of scenarios were changed while writing the manuscript and so are not consistent between the scripts and manuscript: 'bimodal' became 'individual flexibility IF', 'CF' became 'Producers Scroungers PS', 'PEBU' became 'Local site fidelity LSF'.

## ABC
- directory 'analyse field data': script (.Rmd) to calculate summary statistics: spatial pattern statistics and network analyses. Input are 'Distances field.txt' and 'Output field.txt'. 'Summary_stats_field.txt' is the output of the analysis.
- directory 'ABC analysis' contains the actual ABC analysis. It uses output files from 'IBM/analyse runs/Outputs' and 'ABC/analyse field data'. 'ABC_Bembix_INLA.Rmd' is the overall ABC analysis, 'Detailed_model_selection.Rmd' is the code to reproduce table 2 from the manuscript.



Note: a selection of 'real nests' are made within the scripts, these are nests that have multiple visits on several days, where prey were brought to, or where prolonged digging was seen. This species is known for doing test diggings, and we wanted to exclude these from the analyses.
