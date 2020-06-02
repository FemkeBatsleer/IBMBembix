# IBMBembix
source code and data for manuscript

## Microhabitat model
- 'Microhabitat Final Model with iterations.R': source code for the microhabitat model with INLA, including 10 iterations for cross-validation.
- 'Microhabitat Final Model visualisations of iterations.R': source code for making graphs of the microhabitat model. Datafiles needed in this script are made by 'Microhabitat Final Model with iterations.R' or given.
- MicrohabitatmodelComplete.xlsx: datafile needed
- RegPoints.txt

## IBM
- Bembix_model_simple.py --> actual IBM, used by all Bembix_*_s.py
- Bembix_*_s.py: loops for all scenarios to make several runs for a scenario (input: begin number, end number e.g. 0 1000: makes 1000 runs)
- Hist_distances.txt, grey_plot50.jpg --> input for Bembix_model_simple.py
- directory 'analyse runs': two scripts to extract the parameters and calculate summary statistics: spatial pattern statistics and network analyses from the runs of the IBM; H_* scripts are helper scripts used in the two main scripts. The map 'Outputs' gives the outputs for these of the IBM runs made on the HPC@UGent. These are used in the actual ABC analysis.

## ABC
- directory 'analyse field data': script (.Rmd) to calculate summary statistics: spatial pattern statistics and network analyses. Input are 'Distances field.txt' and 'Output field.txt'. 'Summary_stats_field.txt' is the output of the analysis.
- directory 'ABC analysis' contains the actual ABC analysis.



Note: a selection of 'real nests' are made within the scripts, these are nests that have multiple visits on several days, where prey were brought to, or where prolonged digging was seen. This species is known for doing test diggings, and we wanted to exclude these from the analyses.
