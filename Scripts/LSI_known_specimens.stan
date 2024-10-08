/*
This Stan program models LSI values from sites as a multilevel model
with variation due to element portion (anatomical variation).
Three groups are identified (all specimens are known to be from a specific group)
Immature: specimens killed before reaching adult body size
Female: adult-sized females
Males: adult-sized males and/or castrates
LSI values are used in the model by calculating measurements and reference values using measurement error
LSI values are averaged across a specimen to create specimen-specific values for modeling

Modeling mu (mean) parameters:
Only the average body size of females (mu_female) is directly estimated, to prevent label-switching:
Adult male body size is defined as: mu_female + exp(logdelta_male)
Immature body size is defined as: mu_female - exp(logdelta_immature)

Modeling sigma (sd) parameters:
Sigma parameters (sigma_immature, sigma_female, sigma_male) are all allowed to vary independently

All model parameters vary based on the specimen's element portion.
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
N_Element_Portions: number of different element types (can separate proximal/distal), different measurements from the same element will be averaged per specimen
N_Dimensions: number of different measurement types

Specimen observations:
Variables that provide information about the specimens that help structure the mixture model and multilevel structure
Element_Portion: what element is this specimen
Group: what group does the specimen belong (1 = Immature, 2 = female, 3 = male)

Measurement observations:
These variables record the actual metrical observations used to calculate specimen-specific LSI values that are the base of the mixture models.
Includes code for measurement error on the observation of measured specimens (both analytical and reference specimens).
Measurement_obs, Reference_obs: the observed value of the measurement/reference value
Measurement_sd, Reference_sd: uncertainty statement on the value of the measurement/reference value
Dimension: What measurement type the measurement refers to
Specimen: what specimen the measurement belongs to

*/
data {
  //Analysis counts
  int<lower = 1> N_Specimens; //number of specimens
  int<lower = 1> N_Measurements; //number of total measurements
  int<lower = 1> N_Element_Portions; //number of elements
  int<lower = 1> N_Dimensions; //number of measurement types
  //Specimen observations
  array[N_Specimens] int<lower = 1, upper = N_Element_Portions> Element_Portion; //what element is the specimen?
  array[N_Specimens] int<lower = 1, upper = 3> Group; //1 = immature, 2 = female, 3 = male/castrate
	//Measurement observations
	array[N_Measurements] real<lower = 0> Measurement_obs; //observed measurements
	array[N_Measurements] real<lower = 0> Measurement_sd; //uncertainty around that measurement
	array[N_Dimensions] real<lower = 0> Reference_obs; //reference specimen measurements
	array[N_Dimensions] real<lower = 0> Reference_sd; //uncertainty around the measurement
  array[N_Measurements] int<lower = 1, upper = N_Dimensions> Dimension; //which measurement type is this?
	array[N_Measurements] int<lower = 1, upper = N_Specimens> Specimen; // which specimen is this measurement coming from?
}
/* The model parameters defined for the model:
These parameters define the overall relationship
between LSI parameters (mu, sigma) for different measurements.
The model uses a non-centered parameterization to provide
computational stability in the estimation of the posterior distribution
(Betancourt 2017).
*/
parameters {
  //observations (incorporating measurement error)
  vector<lower = 0>[N_Measurements] Measurement;
  vector<lower = 0>[N_Dimensions] Reference;
  vector[N_Specimens] LSI;
  //global hyper-parameters of the mixture model
  // real grand_mu;
  // real grand_logsigma;
  real grand_mu_female;
  real grand_logdelta_immature;
  real grand_logdelta_male;
  real grand_logsigma_immature;
  real grand_logsigma_female;
  real grand_logsigma_male;
  //non-centered parameterization of the model parameters
  //anatomic variation
  vector<lower = 0>[6] element_sigma;
  matrix[6, N_Element_Portions] z_element;
  cholesky_factor_corr[6] L_Rho_element;
}
/*
Transformations of the parameters:
These transformations use the non-centered
parameterization to calculate the element-specific
parameter values.
*/
transformed parameters {
  //element-level variation
  matrix[N_Element_Portions, 6] v_element;
  matrix[6, 6] Rho_element;
  //output variables (un-transformed)
  vector[N_Element_Portions] mu_female;
  vector[N_Element_Portions] logdelta_immature;
  vector[N_Element_Portions] logdelta_male;
  vector[N_Element_Portions] logsigma_immature;
  vector[N_Element_Portions] logsigma_female;
  vector[N_Element_Portions] logsigma_male;
  //output variables (transformed)
  matrix[3, N_Element_Portions] mu;
  matrix[3, N_Element_Portions] sigma;

  //calculating the offsets
  v_element = (diag_pre_multiply(element_sigma, L_Rho_element) * z_element)';
  Rho_element = L_Rho_element * L_Rho_element';

  //calculating the output variables (untransformed)
  mu_female = grand_mu_female + col(v_element, 1);
  logdelta_immature = grand_logdelta_immature + col(v_element, 2);
  logdelta_male = grand_logdelta_male + col(v_element, 3);
  logsigma_immature = grand_logsigma_immature + col(v_element, 4);
  logsigma_female = grand_logsigma_female + col(v_element, 5);
  logsigma_male = grand_logsigma_male + col(v_element, 6);
  
  //transformations
  for(i in 1:N_Element_Portions) {
    mu[1, i] = mu_female[i] - exp(logdelta_immature[i]);
    mu[2, i] = mu_female[i];
    mu[3, i] = mu_female[i] + exp(logdelta_male[i]);
    sigma[1, i] = exp(logsigma_immature[i]);
    sigma[2, i] = exp(logsigma_female[i]);
    sigma[3, i] = exp(logsigma_male[i]);
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
	grand_mu_female ~ normal(-0.1, 0.1);
  grand_logdelta_immature ~ normal(-3, 0.5);
	grand_logdelta_male ~ normal(-3, 0.5);
	grand_logsigma_immature ~ normal(-3, 0.1);
	grand_logsigma_female ~ normal(-3, 0.1);
	grand_logsigma_male ~ normal(-3, 0.1);
	//variability between analytical units
	element_sigma ~ normal(0, 0.05);
	// element_sigma[2:6] ~ normal(0, 0.5);
	to_vector(z_element) ~ normal(0, 1);
	L_Rho_element ~ lkj_corr_cholesky(2);
  //measurement error
  for(i in 1:N_Dimensions) {
    Reference_obs[i] ~ normal(Reference[i], Reference_sd[i]);
  }
  for(i in 1:N_Measurements) {
    Measurement_obs[i] ~ normal(Measurement[i], Measurement_sd[i]);
    LSI_measurement = log(Measurement[i] / Reference[Dimension[i]]);
    LSI_measurement ~ normal(LSI[Specimen[i]], 0.02); //based on intra-element variation from Popkin, et al. (2012)
  }

  //fitting the mixture model
  for(i in 1:N_Specimens) {
    LSI[i] ~ normal(mu[Group[i], Element_Portion[i]], sigma[Group[i], Element_Portion[i]]);
  }
}
/* Generated quantities:
This block can produce any derived quantities that may be of interest.
Currently, it produces transformed versions of model parameters (not accounting for element or interactions)
and specimen-level probabilities of being young, female, or male (specimen_probability)
*/
generated quantities {
  //names of the original mu and sigma values (so that plotting and summarizing functions don't need to be changed)
  vector[N_Element_Portions] mu_immature;
  vector[N_Element_Portions] mu_male;
  vector[N_Element_Portions] sigma_immature;
  vector[N_Element_Portions] sigma_female;
  vector[N_Element_Portions] sigma_male;
  //grand (population-level) parameters
  real grand_mu_immature;
  real grand_mu_male;
  real grand_sigma_immature;
  real grand_sigma_female;
  real grand_sigma_male;
  //renamed parameters
  for(i in 1:N_Element_Portions) {
    mu_immature[i] = mu[1, i];
    mu_male[i] = mu[3, i];
    sigma_immature[i] = sigma[1, i];
    sigma_female[i] = sigma[2, i];
    sigma_male[i] = sigma[3, i];
  }
  //population-level parameter estimates
  grand_mu_immature = grand_mu_female - exp(grand_logdelta_immature);
  grand_mu_male = grand_mu_female + exp(grand_logdelta_male);
  grand_sigma_immature = exp(grand_logsigma_immature);
  grand_sigma_female = exp(grand_logsigma_female);
  grand_sigma_male = exp(grand_logsigma_male);
}
