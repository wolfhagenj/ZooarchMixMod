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
Sigma parameters (immature_LSI_sd, female_LSI_sd, male_LSI_sd) are all allowed to vary independently

Modeling theta (mixing proportion) parameters:
Theta parameters are all modeled together using the "stick-breaking" method,
meaning only two of the parameters are directly modeled.
For stick-breaking explanation, see Stan Manual for Unit Simplex Inverse Transform
(https://mc-stan.org/docs/2_24/reference-manual/simplex-transform-section.html)

These parameters vary based on the specific measurement type and the site of the measurements.
---

This model was written as part of the following paper:
Wolfhagen, J. L. (2023). "Estimating the Ontogenetic Age and Sex Composition of Faunal Assemblages with Bayesian Multilevel Mixture Models." Journal of Archaeological Method and Theory.
(DOI: 10.1007/s10816-023-09611-y)
Associated R code and other Stan scripts can be found at the following OSF Project: 10.17605/OSF.IO/4H9W6
Also found in the following GitHub page: https://github.com/wolfhagenj/ZooarchMixMod

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
N_Sites: number of sites
N_Element_Portions: number of different element types (can separate proximal/distal), different measurements from the same element will be averaged per specimen
N_Dimensions: number of different measuremed dimensions
Immature_Proportion: matrix of the proportion of potentially immature specimens per element portion

Specimen observations:
Variables that provide information about the specimens that help structure the mixture model and multilevel structure
Site: what site does the specimen come from
Element_Portion: what element is this specimen
Immature: could the specimen be considered young? If not, will have a 0% probability of being young, otherwise probability will be based on size and mixing component (theta[1])

Measurement observations:
These variables record the actual metrical observations used to calculate specimen-specific LSI values that are the base of the mixture models.
Includes code for measurement error on the observation of measured specimens (both analytical and reference specimens).
Measurement_obs, Reference_obs: the observed value of the measurement/reference value
Measurement_sd, Reference_sd: uncertainty statement on the value of the measurement/reference value
Dimension: What measurement type the measurement refers to
Specimen: what specimen the measurement belongs to

Demographic observations:
These data are assemblage-level observations of relevant demographic parameters,
which are defined as different combinations of the mixing parameters:
Immature: theta[1] / (theta[1] + theta[2] + theta[3])
Female: theta[2] / (theta[2] + theta[3])

N_Immature_obs, N_Female_obs: The number of relevant observations (there can be multiple observations for a site, or 0 observations for a site)
Immature_obs_sites, Female_obs_sites: What site the observation relates to
Immature_obs, Female_obs: The number of young (or female) specimens, defined by fusion state or morphology
Immature_obs_n, Female_obs_n: The total number of ageable (or sexable) specimens that the observation comes from

NOTE: To have no observations at all, you must still include 1 observation of 0/0
Example:
N_Immature_obs = 1,
Immature_obs_site = 1,
Immature_obs = 0,
Immature_obs_n = 0
N_Female_obs = 1,
Female_obs_site = 1,
Female_obs = 0,
Female_obs_n = 0

This allows the model to run but does not include any observations of the young proportion or adult sex ratio

User-defined prior distributions:
These data define the prior distributions of the biometric data for
component means and sigmas, as earlier defined (one mu and two deltas).
These parameter definitions are from known-age/known-sex specimens,
particularly the Popkin, et al. (2012) Shetland sheep population.

All distributions are defined as normal distributions: c(mean, sigma).
For delta and sigma parameters, this means taking the (natural) log transformation of the distribution.
See associated R script for code on creating prior distributions and lists of prior distributions.

*/
data {
  //Analysis counts
  int<lower = 1> N_Sites; //number of sites
  int<lower = 1> N_Specimens; //number of specimens
  int<lower = 1> N_Measurements; //number of total measurements
  int<lower = 1> N_Element_Portions; //number of elements
  int<lower = 1> N_Dimensions; //number of measurement types
  //Specimen observations
  array[N_Specimens] int<lower = 1, upper = N_Sites> Site; //what site does the specimen come from?
  array[N_Specimens] int<lower = 1, upper = N_Element_Portions> Element_Portion; //what element is the specimen?
	array[N_Specimens] int<lower = 0, upper = 1> Immature; //whether the specimen could be considered young
	matrix<lower = 0, upper = 1>[N_Sites, N_Element_Portions] Immature_Proportion; //the proportion of specimens that are potentially young for each element (impacts the prob(Young) for specimens)
	//Measurement observations
	array[N_Measurements] real<lower = 0> Measurement_obs; //observed measurements
	array[N_Measurements] real<lower = 0> Measurement_sd; //uncertainty around that measurement
	array[N_Dimensions] real<lower = 0> Reference_obs; //reference specimen measurements
	array[N_Dimensions] real<lower = 0> Reference_sd; //uncertainty around the measurement
  array[N_Measurements] int<lower = 1, upper = N_Dimensions> Dimension; //which measurement type is this?
	array[N_Measurements] int<lower = 1, upper = N_Specimens> Specimen; // which specimen is this measurement coming from?
  //Demographic observations
  int<lower = 1> N_Immature_obs;
	array[N_Immature_obs] int<lower = 1, upper = N_Sites> Immature_obs_site;
	array[N_Immature_obs] int<lower = 0> Immature_obs;
	array[N_Immature_obs] int<lower = 0> Immature_obs_n;
	int<lower = 1> N_Female_obs;
	array[N_Female_obs] int<lower = 1, upper = N_Sites> Female_obs_site;
	array[N_Female_obs] int<lower = 0> Female_obs;
	array[N_Female_obs] int<lower = 0> Female_obs_n;
  //User-defined prior distributions
  array[2] real prior_theta_raw_1;
  array[2] real prior_theta_raw_2;
  array[2] real prior_mu_female;
  array[2] real prior_logdelta_immature;
  array[2] real prior_logdelta_male;
  array[2] real prior_logsigma_immature;
  array[2] real prior_logsigma_female;
  array[2] real prior_logsigma_male;
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
  vector<lower = 0>[N_Dimensions] Reference;
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
  //site-level parameters
  vector<lower = 0>[8] site_sigma;
  matrix[8, N_Sites] z_site;
  cholesky_factor_corr[8] L_Rho_site;
  //interaction between the two
  vector<lower = 0>[8] interaction_sigma;
  matrix[8, N_Element_Portions * N_Sites] z_interaction;
  cholesky_factor_corr[8] L_Rho_interaction;
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
  //site-level variation
  matrix[N_Sites, 8] v_site;
  matrix[8, 8] Rho_site;
  //interaction of site and anatomy
  matrix[N_Element_Portions * N_Sites, 8] v_interaction;
  matrix[8, 8] Rho_interaction;
  //output variables (un-transformed)
  array[2] vector[N_Sites] site_theta_raw;
  array[2] matrix[N_Sites, N_Element_Portions] theta_raw;
  matrix[N_Sites, N_Element_Portions] mu_female;
  matrix[N_Sites, N_Element_Portions] logdelta_immature;
  matrix[N_Sites, N_Element_Portions] logdelta_male;
  matrix[N_Sites, N_Element_Portions] logsigma_immature;
  matrix[N_Sites, N_Element_Portions] logsigma_female;
  matrix[N_Sites, N_Element_Portions] logsigma_male;
  //output variables (transformed)
  array[3] vector[N_Sites] site_theta;
  vector[N_Sites] site_p_immature;
  vector[N_Sites] site_theta_female;
  array[3] matrix[N_Sites, N_Element_Portions] theta;
  matrix[N_Sites, N_Element_Portions] mu_immature;
  matrix[N_Sites, N_Element_Portions] mu_male;
  matrix[N_Sites, N_Element_Portions] sigma_immature;
  matrix[N_Sites, N_Element_Portions] sigma_female;
  matrix[N_Sites, N_Element_Portions] sigma_male;
  
  //calculating the offsets
  v_element = (diag_pre_multiply(element_sigma, L_Rho_element) * z_element)';
  Rho_element = L_Rho_element * L_Rho_element';
  v_site = (diag_pre_multiply(site_sigma, L_Rho_site) * z_site)';
  Rho_site = L_Rho_site * L_Rho_site';
  v_interaction = (diag_pre_multiply(interaction_sigma, L_Rho_interaction) * z_interaction)';
  Rho_interaction = L_Rho_interaction * L_Rho_interaction';

  //calculating the output variables (untransformed)
  for(i in 1:N_Sites) {
    site_theta_raw[1, i] = grand_theta_raw[1] + v_site[i, 1];
    site_theta_raw[2, i] = grand_theta_raw[2] + v_site[i, 2];
    for(j in 1:N_Element_Portions) {
      theta_raw[1, i, j] = grand_theta_raw[1] + v_site[i, 1] + v_element[j, 1] + v_interaction[N_Element_Portions * (i - 1) + j, 1];
      theta_raw[2, i, j] = grand_theta_raw[2] + v_site[i, 2] + v_element[j, 2] + v_interaction[N_Element_Portions * (i - 1) + j, 2];
      mu_female[i, j] = grand_mu_female + v_site[i, 3] + v_element[j, 3] + v_interaction[N_Element_Portions * (i - 1) + j, 3];
      logdelta_immature[i, j] = grand_logdelta_immature + v_site[i, 4] + v_element[j, 4] + v_interaction[N_Element_Portions * (i - 1) + j, 4];
      logdelta_male[i, j] = grand_logdelta_male + v_site[i, 5] + v_element[j, 5] + v_interaction[N_Element_Portions * (i - 1) + j, 5];
      logsigma_immature[i, j] = grand_logsigma_immature + v_site[i, 6] + v_element[j, 6] + v_interaction[N_Element_Portions * (i - 1) + j, 6];
      logsigma_female[i, j] = grand_logsigma_female + v_site[i, 7] + v_element[j, 7] + v_interaction[N_Element_Portions * (i - 1) + j, 7];
      logsigma_male[i, j] = grand_logsigma_male + v_site[i, 8] + v_element[j, 8] + v_interaction[N_Element_Portions * (i - 1) + j, 8];
    }
  }
  mu_immature = mu_female - exp(logdelta_immature);
  mu_male = mu_female + exp(logdelta_male);
  sigma_immature = exp(logsigma_immature);
  sigma_female = exp(logsigma_female);
  sigma_male = exp(logsigma_male);
  //theta transformations
  for(i in 1:N_Sites) {
    site_theta[1, i] = inv_logit(site_theta_raw[1, i] + log(0.5));
    site_theta[2, i] = (1 - site_theta[1, i]) * inv_logit(site_theta_raw[2, i] + log(1));
    site_theta[3, i] = 1 - sum(site_theta[1:2, i]);
    site_p_immature[i] = site_theta[1, i];
    site_theta_female[i] = site_theta[2, i] / sum(site_theta[2:3, i]);
    for(j in 1:N_Element_Portions) {
      theta[1, i, j] = inv_logit(theta_raw[1, i, j] + log(0.5));
      theta[2, i, j] = (1 - theta[1, i, j]) * inv_logit(theta_raw[2, i, j] + log(1));
      theta[3, i, j] = 1 - sum(theta[1:2, i, j]);
    }
  }
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
	element_sigma[1:2] ~ normal(0, 0.5); //less reason to think that mixture proportions are the same across elements
	element_sigma[3:8] ~ normal(0, 0.05);
	to_vector(z_element) ~ normal(0, 1);
	L_Rho_element ~ lkj_corr_cholesky(2);
	//
	site_sigma[1:2] ~ normal(0, 0.5); //less reason to think that mixture proportions are the same across elements
	site_sigma[3:5] ~ normal(0, 0.1);
	site_sigma[6:8] ~ normal(0, 0.05);
	to_vector(z_site) ~ normal(0, 1);
	L_Rho_site ~ lkj_corr_cholesky(2);
	//
	interaction_sigma[1:2] ~ normal(0, 0.25); //less reason to think that mixture proportions are the same across elements
	interaction_sigma[3:5] ~ normal(0, 0.05);
	interaction_sigma[6:8] ~ normal(0, 0.05);
	to_vector(z_interaction) ~ normal(0, 1);
	L_Rho_interaction ~ lkj_corr_cholesky(2);

  //measurement error
  Measurement_obs ~ normal(Measurement, Measurement_sd);
  Reference_obs ~ normal(Reference, Reference_sd);
  for(i in 1:N_Measurements) {
    LSI_measurement = log(Measurement[i]) - log(Reference[Dimension[i]]);
    LSI_measurement ~ normal(LSI[Specimen[i]], 0.02); //based on intra-element variation from Popkin, et al. (2012)
  }
  //observation of age and sex parameters
  for(i in 1:N_Immature_obs) {
    Immature_obs[i] ~ binomial(Immature_obs_n[i], site_p_immature[Immature_obs_site[i]]);
  }
  for(i in 1:N_Female_obs) {
    Female_obs[i] ~ binomial(Female_obs_n[i], site_theta_female[Female_obs_site[i]]);
  }
  
  //fitting the mixture model
  for(i in 1:N_Specimens) {
    if(Immature[i] == 1) {
      vector[3] lps;
      lps[1] = log(theta[1, Site[i], Element_Portion[i]] / Immature_Proportion[Site[i], Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_immature[Site[i], Element_Portion[i]], sigma_immature[Site[i], Element_Portion[i]]);
      lps[2] = log(theta[2, Site[i], Element_Portion[i]] * Immature_Proportion[Site[i], Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_female[Site[i], Element_Portion[i]], sigma_female[Site[i], Element_Portion[i]]);
      lps[3] = log(theta[3, Site[i], Element_Portion[i]] * Immature_Proportion[Site[i], Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_male[Site[i], Element_Portion[i]], sigma_male[Site[i], Element_Portion[i]]);
      target += log_sum_exp(lps);
    }
    if(Immature[i] == 0) {
      vector[2] lps;
      lps[1] = log(theta[2, Site[i], Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_female[Site[i], Element_Portion[i]], sigma_female[Site[i], Element_Portion[i]]);
      lps[2] = log(theta[3, Site[i], Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_male[Site[i], Element_Portion[i]], sigma_male[Site[i], Element_Portion[i]]);
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
  vector[3] grand_theta;
  real grand_p_immature;
  real grand_theta_female;
  real grand_mu_immature;
  real grand_mu_male;
  real grand_sigma_immature;
  real grand_sigma_female;
  real grand_sigma_male;
  //site-level parameters
  vector[N_Sites] site_mu_immature;
  vector[N_Sites] site_mu_female;
  vector[N_Sites] site_mu_male;
  vector[N_Sites] site_sigma_immature;
  vector[N_Sites] site_sigma_female;
  vector[N_Sites] site_sigma_male;
  //specimen-level probabilities of being young, female, or male
  array[N_Specimens] simplex[3] specimen_prob;

  //population-level parameter estimates
  grand_theta[1] = inv_logit(grand_theta_raw[1] + log(0.5));
  grand_theta[2] = (1 - grand_theta[1]) * inv_logit(grand_theta_raw[2] + log(1));
  grand_theta[3] = 1 - sum(grand_theta[1:2]);
  grand_p_immature = grand_theta[1];
  grand_theta_female = grand_theta[2] / sum(grand_theta[2:3]);
  grand_mu_immature = grand_mu_female - exp(grand_logdelta_immature);
  grand_mu_male = grand_mu_female + exp(grand_logdelta_male);
  grand_sigma_immature = exp(grand_logsigma_immature);
  grand_sigma_female = exp(grand_logsigma_female);
  grand_sigma_male = exp(grand_logsigma_male);
  //site-level parameter estimates
  site_mu_female = grand_mu_female + col(v_site, 3);
  site_mu_immature = site_mu_female - exp(grand_logdelta_immature + col(v_site, 4));
  site_mu_male = site_mu_female + exp(grand_logdelta_male + col(v_site, 5));
  site_sigma_immature = exp(grand_logsigma_immature + col(v_site, 6));
  site_sigma_female = exp(grand_logsigma_female + col(v_site, 7));
  site_sigma_male = exp(grand_logsigma_male + col(v_site, 8));

  //individual specimen mixture probabilities
  for(i in 1:N_Specimens) {
    if(Immature[i] == 1) {
      vector[3] lps;
      lps[1] = log(theta[1, Site[i], Element_Portion[i]] / Immature_Proportion[Site[i], Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_immature[Site[i], Element_Portion[i]], sigma_immature[Site[i], Element_Portion[i]]);
      lps[2] = log(theta[2, Site[i], Element_Portion[i]] * Immature_Proportion[Site[i], Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_female[Site[i], Element_Portion[i]], sigma_female[Site[i], Element_Portion[i]]);
      lps[3] = log(theta[3, Site[i], Element_Portion[i]] * Immature_Proportion[Site[i], Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_male[Site[i], Element_Portion[i]], sigma_male[Site[i], Element_Portion[i]]);
      specimen_prob[i][1] = exp(lps[1] - log_sum_exp(lps));
      specimen_prob[i][2] = exp(lps[2] - log_sum_exp(lps));
      specimen_prob[i][3] = exp(lps[3] - log_sum_exp(lps));
    }
    if(Immature[i] == 0) {
      vector[2] lps;
      lps[1] = log(theta[2, Site[i], Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_female[Site[i], Element_Portion[i]], sigma_female[Site[i], Element_Portion[i]]);
      lps[2] = log(theta[3, Site[i], Element_Portion[i]]) + normal_lpdf(LSI[i] | mu_male[Site[i], Element_Portion[i]], sigma_male[Site[i], Element_Portion[i]]);
      specimen_prob[i][1] = 0;
      specimen_prob[i][2] = exp(lps[1] - log_sum_exp(lps));
      specimen_prob[i][3] = exp(lps[2] - log_sum_exp(lps));
    }
  }
}
