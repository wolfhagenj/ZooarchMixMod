####Single Assemblage Mixture Modeling####
#This script is part of the ZooarchMixMod GitHub project (https://github.com/wolfhagenj/ZooarchMixMod)
#It provides instructions to take standardized faunal measurements from a single assemblage and model
#them as a mixture of immature, adult-sized female, and adult-sized male specimens using a Bayesian
#multilevel mixture model.

#This code is meant to be copied and adapted for your own personal use. While it is designed for data
#set up like the OpenContext data used in the original manuscript, the same principles can be adapted
#to data set up in a different way. The principles can even be used to standardize the data yourself
#outside of R. See the ReadMe for the GitHub repository to see what the model expects in terms of
#data setup.

####Load the packages####
#These packages are necessary for the script. Use the following line (without the comment #) to install
#the packages if you have not already installed them:
#install.packages("data.table", "cmdstanr", "rstan")
#For "cmdstanr", see the following installation instructions: https://mc-stan.org/cmdstanr/articles/cmdstanr.html
#For "rstan", see the following installation instructions: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
#In general, you will need to install rtools to have a C++ toolchain; the above links help provide steps to do so.
library("cmdstanr")
library("data.table")
library("rstan")

####Import the data (if necessary)####
#Import mixture model data. This script uses the results of the "data_cleaning.R" as an example
#Alternatively, you could run the script in the same R environment earlier so that the objects exist
single_assemblage_mixmod_data <- fread("./Data/Site Mixture Model Data.csv")
single_assemblage_demographic_observations <- fread("./Data/Site Demographic Observations.csv")

####Set up the data for the Stan model####
#Stan requires the data set up as a list object
single_assemblage_mixmod_standata <- list(
  #Sample sizes
  N_Specimens = single_assemblage_mixmod_data[, .N, Specimen_No][, .N],
  N_Measurements = single_assemblage_mixmod_data[, .N],
  N_Element_Portions = single_assemblage_mixmod_data[, .N, Element_Portion][, .N],
  N_Measurement_Sets = single_assemblage_mixmod_data[, .N, Measurement_Set][, .N],
  #Specimen observations
  Element_Portion = single_assemblage_mixmod_data[, .N, .(Specimen_No, Element_Portion, Immature)][order(Specimen_No), Element_Portion],
  Immature = single_assemblage_mixmod_data[, .N, .(Specimen_No, Element_Portion, Immature)][order(Specimen_No), Immature],
  Immature_Proportion = as.matrix(dcast(single_assemblage_mixmod_data[, .(Immature_Proportion = mean(Immature)), .(Site_No, Element_Portion)], Site_No ~ Element_Portion, value.var = "Immature_Proportion", fill = 0))[, 2:(single_assemblage_mixmod_data[, .N, Element_Portion][, .N] + 1)],
  #Measurement observations
  Measurement_obs = single_assemblage_mixmod_data[, Measurement_value],
  Measurement_sd = single_assemblage_mixmod_data[, Measurement_value * 0.01], #Calculate measurement error for observed measurements and reference data (1% based on data from Breslawksi and Byers 2015)
  Reference_obs = single_assemblage_mixmod_data[, .N, .(Measurement_Set, Reference_value)][order(Measurement_Set), Reference_value],
  Reference_sd = single_assemblage_mixmod_data[, .N, .(Measurement_Set, Reference_value)][order(Measurement_Set), Reference_value * 0.01],
  Measurement_Set = single_assemblage_mixmod_data[, Measurement_Set],
  Specimen = single_assemblage_mixmod_data[, Specimen_No],
  #Demographic observations
  Immature_obs = single_assemblage_demographic_observations[, N_Unfused],
  Immature_obs_n = single_assemblage_demographic_observations[, N_Ageable],
  Female_obs = single_assemblage_demographic_observations[, N_Female],
  Female_obs_n = single_assemblage_demographic_observations[, N_Sexable],
  #Prior distributions for hyper-parameters
  #These can be changed to reflect different scenarios
  #They are largely based on the Shetland sheep population
  prior_theta_raw_1 = c(-0.5, 1.5),
  prior_theta_raw_2 = c(0, 1.5),
  prior_mu_female = c(0, 0.1),
  prior_logdelta_immature = c(-3.5, 0.4),
  prior_logdelta_male = c(-2.7, 0.1),
  prior_logsigma_immature = c(-3.05, 0.1),
  prior_logsigma_female = c(-3.1, 0.1),
  prior_logsigma_male = c(-3.1, 0.1)
)
singlesite_mixture_stanmodel <- cmdstan_model("./Scripts/LSI_mixture_model_singlesite.stan")
single_assemblage_samples <- singlesite_mixture_stanmodel$sample(
  data = single_assemblage_mixmod_standata,
  chains = 4,
  parallel_chains = 4,
  refresh = 250,
  adapt_delta = 0.80,
  max_treedepth = 15
)
#IMPORTANT: Divergent transitions require rerunning the model, increase the adapt_delta argument
#see https://mc-stan.org/misc/warnings for more details

single_assemblage_stanfit <- rstan::read_stan_csv(single_assemblage_samples$output_files())
single_assemblage_post <- extract(single_assemblage_stanfit)
