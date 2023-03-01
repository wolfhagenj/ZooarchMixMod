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
#If you do this, you need to use the following lines (uncommented)
#single_assemblage_mixmod_data <- site_mixmod_data
#single_assemblage_demographic_observations <- site_demographic_observations
single_assemblage_mixmod_data <- fread("./Data/Site Mixture Model Data.csv")
single_assemblage_demographic_observations <- fread("./Data/Site Demographic Observations.csv")

#Add Site_No category to mixmod data (used in defining the immature proportions)
#This variable isn't included in the data_cleaning script
single_assemblage_mixmod_data[, Site_No := 1]

####Set up the data for the Stan model####
#Stan requires the data set up as a list object
single_assemblage_mixmod_standata <- list(
  #Sample sizes
  N_Specimens = single_assemblage_mixmod_data[, .N, Specimen_No][, .N],
  N_Measurements = single_assemblage_mixmod_data[, .N],
  N_Element_Portions = single_assemblage_mixmod_data[, .N, Element_Portion][, .N],
  N_Dimensions = single_assemblage_mixmod_data[, .N, Dimension][, .N],
  #Specimen observations
  Element_Portion = single_assemblage_mixmod_data[, .N, .(Specimen_No, Element_Portion, Immature)][order(Specimen_No), Element_Portion],
  Immature = single_assemblage_mixmod_data[, .N, .(Specimen_No, Element_Portion, Immature)][order(Specimen_No), Immature],
  Immature_Proportion = as.matrix(dcast(single_assemblage_mixmod_data[, .(Immature_Proportion = mean(Immature)), .(Site_No, Element_Portion)], Site_No ~ Element_Portion, value.var = "Immature_Proportion", fill = 0))[, 2:(single_assemblage_mixmod_data[, .N, Element_Portion][, .N] + 1)],
  #Measurement observations
  Measurement_obs = single_assemblage_mixmod_data[, Measurement_value],
  Measurement_sd = single_assemblage_mixmod_data[, Measurement_value * 0.01], #Calculate measurement error for observed measurements and reference data (1% based on data from Breslawksi and Byers 2015)
  Reference_obs = single_assemblage_mixmod_data[, .N, .(Dimension, Reference_value)][order(Dimension), Reference_value],
  Reference_sd = single_assemblage_mixmod_data[, .N, .(Dimension, Reference_value)][order(Dimension), Reference_value * 0.01],
  Dimension = single_assemblage_mixmod_data[, Dimension],
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
#IMPORTANT: Warnings about divergent transitions typically requires rerunning the model
#Increase the adapt_delta argument to address the issue
#see https://mc-stan.org/misc/warnings for more details

####Evaluating model fit####
#There are no hard-and-fast rules to determine whether your model has fit correctly (or converged).
#You need multiple chains to evaluate this: these are separate runs of the model at different starting
#points. Ideally, you're getting the same results regardless of where you started.
#See https://mc-stan.org/rstan/reference/Rhat.html for more information
#The variables you can examine to argue that you have a correctly-fit model are:
#rhat: compares variation between and within chains; values should be near 1.00
#ess_bulk: estimates sample size (of the posterior draws) using rank normalized draws. Can be higher than
#the actual number of posterior draws; this roughly shows how many "effectively" independent samples you have
#to produce your mean and median estimates
#ess_tail: the same process, though focused on the quantile estimates (q5, q95) and variance.
#See also https://betanalpha.github.io/assets/case_studies/rstan_workflow.html for more information

#These are the names of the model hyper-parameters that were given user-defined prior distributions
single_assemblage_samples$summary(c("grand_theta_raw",
                                    "grand_mu_female", "grand_logdelta_immature", "grand_logdelta_male",
                                    "grand_logsigma_immature", "grand_logsigma_female", "grand_logsigma_male"))
#This displays an estimate of the parameter for each element portion
single_assemblage_samples$summary(c("mu_female"))

#Beyond summaries, it is also worth looking at trace plots that overlay each chain on top of one another
#Ideally, all of the chains should overlap each other and there shouldn't be any directionality across the
#length of the chains. This is sometimes called looking for "fuzzy caterpillars".
#To do this, first save the model fit into an rstan stanfit object
single_assemblage_stanfit <- rstan::read_stan_csv(single_assemblage_samples$output_files())

#Then use the convenient plotting function for a defined variable
traceplot(single_assemblage_stanfit, c("grand_theta_raw",
                                       "grand_mu_female", "grand_logdelta_immature", "grand_logdelta_male",
                                       "grand_logsigma_immature", "grand_logsigma_female", "grand_logsigma_male"))
#
traceplot(single_assemblage_stanfit, "mu_female")

####Exporting the posterior distributions####
#For posterior analysis, you'll want to have the exact posterior estimates of the parameter values.
#These can be extracted from the stanfit object
single_assemblage_post <- extract(single_assemblage_stanfit)

#This object can then be saved to be uploaded later or used directly in another script as an open R object.

#IMPORTANT: To ensure reproducibility, you can also extract the random seed used to generate your model
#and then place the seed in the samples call using the argument
#seed = [whatever number you were given]
get_seed(single_assemblage_stanfit)
#This would provide strict reproducibility, as you would get the exact same posterior samples.
#Overall results, however, shouldn't depend on the specific seed you use.