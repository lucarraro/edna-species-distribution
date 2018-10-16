# edna-species-distribution
Matlab code and data files supporting Carraro, L. et al. (2018) "Estimating species distribution and abundance in river networks using environmental DNA", Proceedings of the National Academy of Science.

List of files contained:

1) Data
- data_wigger.mat: morphological data for the Wigger catchment
- eDNA_data.mat: measured eDNA concentrations for both species
- data_explanation.xslx: definition of all data contained in data_wigger.mat, eDNA_data.mat

2) Main codes
- RUN_MODEL.m: it reads data for F. sultana, executes the model and performs the Metropolis-within-Gibbs sampling. Run a second time by changing 'Fs' to Tb' to run the model for T. bryosalmonae.
- ANALYSE_DATA.m: it reads 'results_Fs.mat' and 'results_Tb.mat' (result files produced by RUN_MODEL.m) and draws graphs

3) Auxiliary functions:
- cov_matrix.m: build matrix of covariates used to calculate eDNA productions
- eval_concups.m: solve Eq. (1) of the main text
- eval_posterior.m: evaluate posterior distribution by using the unreferenced equation of the Materials and Methods section
- TruncNormRnd.m: draw values from a truncated normal distribution
- v2struct.m: pack/pnpack variables to/from a scalar structure (credits to N. Adi)
