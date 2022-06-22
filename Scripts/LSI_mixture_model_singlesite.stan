/*
This Stan program models LSI values from sites as a multilevel model
with variation due to measurement type (anatomical variation) and sites.
LSI values are modeled as a mixture model of 3 normal distributions,
representing adult males, adult females, and young animals (killed before reaching adult body size).

Modeling mu (mean) parameters:
Only female_LSI_mean is directly estimated, to prevent label-switching:
Adult male body size is defined as: female_LSI_mean + delta_LSI_male
Young body size is defined as: female_LSI_mean - delta_LSI_immature

Modeling sigma (sd) parameters:
Sigma parameters (Immature_LSI_sd, female_LSI_sd, male_LSI_sd) are all allowed to vary independently

Modeling theta (mixing proportion) parameters:
Theta parameters are all modeled together using the "stick-breaking" method,
meaning only two of the parameters are directly modeled.
For stick-breaking explanation, see Stan Manual for Unit Simplex Inverse Transform
(https://mc-stan.org/docs/2_24/reference-manual/simplex-transform-section.html)

These parameters vary based on the specific measurement type and the site of the measurements.
---

This model was written as part of the following paper:
[!!!PAPER]
(DOI: [!!!DOI])
Associated R code and other Stan scripts can be found at the following OSF Project: [!!!DOI LINK]

Learn more about model development with Stan at:
    http://mc-stan.org/users/interfaces/rstan.html
    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
*/

/* The input data:
This block defines the input data types for a mixture model of LSI values that uses
measurement error and multilevel modeling to relate parameters from different measurement types
and sites.

Analysis counts:
These variables define the various analytical units in the model fit:
N_Specimens: number of individual measured bones (each bone can have multiple measurements)
N_Measurements: number of specific measurements used in the model
N_Element_Portions: number of different element types (can separate proximal/distal), different measurements from the same element will be averaged per specimen
N_Measurement_Sets: number of different measurement types

Specimen observations:
Variables that provide information about the specimens that help structure the mixture model and multilevel structure
Element: what element is this specimen
Immature: could the specimen be considered immature? If not, will have a 0% probability of being immature, otherwise probability will be based on size and mixing component (theta[1])

Measurement observations:
These variables record the actual metrical observations used to calculate specimen-specific LSI values that are the base of the mixture models.
Includes code for measurement error on the observation of measured specimens (both analytical and reference specimens).
Measurement_obs, Reference_obs: the observed value of the measurement/reference value
Measurement_sd, Reference_sd: uncertainty statement on the value of the measurement/reference value
Measurement_Set: What measurement set the measurement refers to
Specimen: what specimen the measurement belongs to

Demographic observations:
These data are assemblage-level observations of relevant demographic parameters,
which are defined as different combinations of the mixing parameters:
Immature: theta[1] / (theta[2] + theta[3])
Female: theta[2] / theta[3]

Immature_obs, Female_obs: The number of young (or female) specimens, defined by fusion state or morphology
Immature_obs_n, Female_obs_n: The total number of ageable (or sexable) specimens that the observation comes from

NOTE: To have no observations at all, you must still include 1 observation of 0/0
Example:
Immature_obs = 0,
Immature_obs_n = 0
Female_obs = 0,
Female_obs_n = 0

This allows the model to run but does not include any observations of the proportion of immature specimens or adult sex ratio

User-defined prior distributions:
These data define the prior distributions of the biometric data for
the model hyper-parameters (not their transformed values)

All distributions are defined as normal distributions: c(mean, sigma).
For delta and sigma parameters, this means taking the (natural) log transformation of the distribution.
See associated R script for code on creating prior distributions and lists of prior distributions.

*/
data {
  //Analysis counts
  int<lower = 1> N_Specimens; //number of specimens
  int<lower = 1> N_Measurements; //number of total measurements
  int<lower = 1> N_Element_Portions; //number of elements
  int<lower = 1> N_Measurement_Sets; //number of measurement types
  real<lower = 0, upper = 1> Immature_Proportion[N_Element_Portions]; //proportion of potentially young specimens
  //Specimen observations
  int<lower = 1, upper = N_Element_Portions> Element_Portion[N_Specimens]; //what element is the specimen?
	int<lower = 0, upper = 1> Immature[N_Specimens]; //whether the specimen could be considered young
	//Measurement observations
	real<lower = 0> Measurement_obs[N_Measurements]; //observed measurements
	real<lower = 0> Measurement_sd[N_Measurements]; //uncertainty around that measurement
	real<lower = 0> Reference_obs[N_Measurement_Sets]; //reference specimen measurements
	real<lower = 0> Reference_sd[N_Measurement_Sets]; //uncertainty around the measurement
  int<lower = 1, upper = N_Measurement_Sets> Measurement_Set[N_Measurements]; //which measurement type is this?
	int<lower = 1, upper = N_Specimens> Specimen[N_Measurements]; // which specimen is this measurement coming from?
  //Demographic observations
	int<lower = 0> Immature_obs;
	int<lower = 0> Immature_obs_n;
	int<lower = 0> Female_obs;
	int<lower = 0> Female_obs_n;
  //User-defined prior distributions
  real prior_theta_raw_1[2];
  real prior_theta_raw_2[2];
  real prior_mu_female[2];
  real prior_logdelta_immature[2];
  real prior_logdelta_male[2];
  real prior_logsigma_immature[2];
  real prior_logsigma_female[2];
  real prior_logsigma_male[2];
}
/* The model parameters defined for the model:
These parameters define the overall relationship
between LSI parameters (mu, sigma) for different measurements.
The model uses a non-centered parameterization to provide
computational stability in the estimation of the posterior distribution
(Betancourt 2017).
There are two terms (measurement and site) and a third interaction term
that allows for parameters to vary between different sites and measurements
and for those differences to vary between sites
(e.g., mean LSI for Scapula_GLP will not always be larger than mean LSI for
Humerus_Bd at all sites)
*/
parameters {
  //observations (incorporating measurement error)
  vector<lower = 0>[N_Measurements] Measurement;
  vector<lower = 0>[N_Measurement_Sets] Reference;
  vector[N_Specimens] LSI;
  //global hyper-parameters of the mixture model
  vector[2] grand_theta_raw; //will be constrained to theta
  real grand_logdelta_immature; //log-scale
  real grand_mu_female;
  real grand_logdelta_male; //log-scale
  real grand_logsigma_immature; //log-scale
  real grand_logsigma_female; //log-scale
  real grand_logsigma_male; //log-scale
  //non-centered parameterization of the model parameters
  //anatomic variation
  vector<lower = 0>[8] element_sigma;
  matrix[8, N_Element_Portions] z_element;
  cholesky_factor_corr[8] L_Rho_element;
}
/*
Transformations of the parameters:
These transformations use the non-centered
parameterization to calculate the measurement-specific
parameter values.
*/
transformed parameters {
  //element-level variation
  matrix[N_Element_Portions, 8] v_element;
  matrix[8, 8] Rho_element;
  //output variables (un-transformed)
  vector[N_Element_Portions] theta_raw[2];
  vector[N_Element_Portions] mu_female;
  vector[N_Element_Portions] logdelta_immature;
  vector[N_Element_Portions] logdelta_male;
  vector[N_Element_Portions] logsigma_immature;
  vector[N_Element_Portions] logsigma_female;
  vector[N_Element_Portions] logsigma_male;
  //output variables (transformed)
  vector[N_Element_Portions] theta[3];
  vector[N_Element_Portions] mu_immature;
  vector[N_Element_Portions] mu_male;
  vector[N_Element_Portions] sigma_immature;
  vector[N_Element_Portions] sigma_female;
  vector[N_Element_Portions] sigma_male;
  vector[3] grand_theta;
  real p_immature;
  real theta_female;
  
  //calculating the offsets
  v_element = (diag_pre_multiply(element_sigma, L_Rho_element) * z_element)';
  Rho_element = L_Rho_element * L_Rho_element';

  //calculating the output variables (untransformed)
  theta_raw[1, ] = grand_theta_raw[1] + col(v_element, 1);
  theta_raw[2, ] = grand_theta_raw[2] + col(v_element, 2);
  mu_female = grand_mu_female + col(v_element, 3);
  logdelta_immature = grand_logdelta_immature + col(v_element, 4);
  logdelta_male = grand_logdelta_male + col(v_element, 5);
  logsigma_immature = grand_logsigma_immature + col(v_element, 6);
  logsigma_female = grand_logsigma_female + col(v_element, 7);
  logsigma_male = grand_logsigma_male + col(v_element, 8);
  
  //transformations
  mu_immature = mu_female - exp(logdelta_immature);
  mu_male = mu_female + exp(logdelta_male);
  sigma_immature = exp(logsigma_immature);
  sigma_female = exp(logsigma_female);
  sigma_male = exp(logsigma_male);
  //theta transformations
  grand_theta[1] = inv_logit(grand_theta_raw[1] + log(0.5));
  grand_theta[2] = (1 - grand_theta[1]) * inv_logit(grand_theta_raw[2] + log(1));
  grand_theta[3] = 1 - sum(grand_theta[1:2]);
  for(i in 1:N_Element_Portions) {
    theta[1, i] = inv_logit(theta_raw[1, i] + log(0.5));
    theta[2, i] = (1 - theta[1, i]) * inv_logit(theta_raw[2, i] + log(1));
    theta[3, i] = 1 - sum(theta[1:2, i]);
  }
  p_immature = grand_theta[1];
  theta_female = grand_theta[2] / sum(grand_theta[2:3]);
}
/* The estimated model:
First, this block models hyper-parameters based
on broad prior distributions.
Then this block models the observations to incorporate measurement error.
This includes the demographic observations that can inform site_theta parameters.
Finally, this block models specimen-specific LSI values as coming from a element- and site-specific
mixture model of normal distributions. If the specimen is considered potentially young, this is a
three-member mixture model, otherwise it is a two-member mixture model.
*/
model {
  real LSI_measurement; //dummy variable used to help calculate specimen-specific LSI values
  //distributions of hyper-parameters and multilevel model parameters
  grand_theta_raw[1] ~ normal(prior_theta_raw_1[1], prior_theta_raw_1[2]);
  grand_theta_raw[2] ~ normal(prior_theta_raw_2[1], prior_theta_raw_2[2]);
	grand_logdelta_immature ~ normal(prior_logdelta_immature[1], prior_logdelta_immature[2]);
	grand_mu_female ~ normal(prior_mu_female[1], prior_mu_female[2]);
	grand_logdelta_male ~ normal(prior_logdelta_male[1], prior_logdelta_male[2]);
	grand_logsigma_immature ~ normal(prior_logsigma_immature[1], prior_logsigma_immature[2]);
	grand_logsigma_female ~ normal(prior_logsigma_female[1], prior_logsigma_female[2]);
	grand_logsigma_male ~ normal(prior_logsigma_male[1], prior_logsigma_male[2]);
	//variability between analytical units
	element_sigma[1:2] ~ normal(0, 1); //less reason to think that mixture proportions are the same across elements
	element_sigma[3] ~ normal(0, 0.5);
	element_sigma[4:8] ~ normal(0, 0.5); //0.25 shouldn't be as much variation between elements for these parameters, though
	to_vector(z_element) ~ normal(0, 1);
	L_Rho_element ~ lkj_corr_cholesky(2);
	
  //measurement error
  Measurement_obs ~ normal(Measurement, Measurement_sd);
  Reference_obs ~ normal(Reference, Reference_sd);
  for(i in 1:N_Measurements) {
    LSI_measurement = log(Measurement[i] / Reference[Measurement_Set[i]]);
    LSI_measurement ~ normal(LSI[Specimen[i]], 0.02); //based on intra-element variation from Popkin, et al. (2012)
  }
  //observation of age and sex parameters
  Immature_obs ~ binomial(Immature_obs_n, p_immature);
  Female_obs ~ binomial(Female_obs_n, theta_female);

  //fitting the mixture model
  for(i in 1:N_Specimens) {
    if(Immature[i] == 1) {
      vector[3] lps;
      lps[1] = log(theta[1, Element_Portion[i]] / Immature_Proportion[Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_immature[Element_Portion[i]], sigma_immature[Element_Portion[i]]);
      lps[2] = log(theta[2, Element_Portion[i]] * Immature_Proportion[Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_female[Element_Portion[i]], sigma_female[Element_Portion[i]]);
      lps[3] = log(theta[3, Element_Portion[i]] * Immature_Proportion[Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_male[Element_Portion[i]], sigma_male[Element_Portion[i]]);
      target += log_sum_exp(lps);
    }
    if(Immature[i] == 0) {
      vector[2] lps;
      lps[1] = log(theta[2, Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_female[Element_Portion[i]], sigma_female[Element_Portion[i]]);
      lps[2] = log(theta[3, Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_male[Element_Portion[i]], sigma_male[Element_Portion[i]]);
      target += log_sum_exp(lps);
    }
  }
}
/* Generated quantities:
This block can produce any derived quantities that may be of interest.
Currently, it produces transformed versions of model parameters (not accounting for element or interactions)
and specimen-level probabilities of being young, female, or male (specimen_probability)
*/
generated quantities {
  //grand (population-level) parameters
  real grand_mu_immature;
  real grand_mu_male;
  real grand_sigma_immature;
  real grand_sigma_female;
  real grand_sigma_male;
  //specimen-level probabilities of being young, female, or male
  simplex[3] specimen_prob[N_Specimens];

  //population-level parameter estimates
  grand_mu_immature = grand_mu_female - exp(grand_logdelta_immature);
  grand_mu_male = grand_mu_female + exp(grand_logdelta_male);
  grand_sigma_immature = exp(grand_logsigma_immature);
  grand_sigma_female = exp(grand_logsigma_female);
  grand_sigma_male = exp(grand_logsigma_male);

  //individual specimen mixture probabilities
  for(i in 1:N_Specimens) {
    if(Immature[i] == 1) {
      vector[3] lps;
      lps[1] = log(theta[1, Element_Portion[i]] / Immature_Proportion[Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_immature[Element_Portion[i]], sigma_immature[Element_Portion[i]]);
      lps[2] = log(theta[2, Element_Portion[i]] * Immature_Proportion[Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_female[Element_Portion[i]], sigma_female[Element_Portion[i]]);
      lps[3] = log(theta[3, Element_Portion[i]] * Immature_Proportion[Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_male[Element_Portion[i]], sigma_male[Element_Portion[i]]);
      specimen_prob[i][1] = exp(lps[1] - log_sum_exp(lps));
      specimen_prob[i][2] = exp(lps[2] - log_sum_exp(lps));
      specimen_prob[i][3] = exp(lps[3] - log_sum_exp(lps));
    }
    if(Immature[i] == 0) {
      vector[2] lps;
      lps[1] = log(theta[2, Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_female[Element_Portion[i]], sigma_female[Element_Portion[i]]);
      lps[2] = log(theta[3, Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_male[Element_Portion[i]], sigma_male[Element_Portion[i]]);
      specimen_prob[i][1] = 0;
      specimen_prob[i][2] = exp(lps[1] - log_sum_exp(lps));
      specimen_prob[i][3] = exp(lps[2] - log_sum_exp(lps));
    }
  }
}
