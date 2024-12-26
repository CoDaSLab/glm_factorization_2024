# Supporting code for ''Considerations for missing data, outliers and transformations in permutation testing for ANOVA  with multivariate responses''

*Oliver Polushkina Merchanskaya, Michael D. Sorochan Armstrong, Carolina Gómez Llorente, Patricia Ferrer, Sergi Fernandez-Gonzalez, Miriam Perez-Cruz, María Dolores Gómez-Roig, José Camacho*

This repository requires the MEDA toolbox v1.4 at https://github.com/josecamachop/MEDA-Toolbox/releases/tag/v1.4


Code related to the simulations for missing data in Section 5 of the paper. As a function of an increasing percentage of missing data, a comparison between unconditional mean replacement (UMR), conditional mean replacement (CMR) and permutational conditional mean replacement (pCMR) are considered using different missingness patterns and mechanisms.

## example_table

Code related to the simulatino of two factors, with one known _a-priori_ to be significant and another insignificant. Comparison of the drops in Sum of Squares (SS) as a function of the different missing value imputation methods for 5% missing data.

## missing_data_mcar_rand

Results in Fig. 1. Comparison of the relative performance of UMR, CMR, and pCMR as a function of the percent of missing data in synthetic data. Missingness is introduced following a haphazard pattern following an MCAR mechanism.

## missing_data_mar_rand

Results in Fig. 2. Comparison of the relative performance of UMR, CMR, and pCMR as a function of the percent of missing data in synthetic data. Missingness is introduced following a haphazard pattern following an MAR mechanism.

## missing_data_mar_mono

Results in Fig. 3. Comparison of the relative performance of UMR, CMR, and pCMR as a function of the percent of missing data in synthetic data. Missingness is introduced following a monotone pattern following an MCAR mechanism, closely resembling attrition.