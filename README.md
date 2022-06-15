# Estimating the Age and Sex Composition of Faunal Assemblages with Bayesian Multilevel Mixture Models
This GitHub repository maintains the code, data, and RMarkdown files necessary to reproduce the pre-print of an article prepared for submission to the *Journal of Archaeological Method and Theory*, written by Jesse Wolfhagen. The method relies on [Stan](https://mc-stan.org/) for Bayesian analysis and is written for use in the R statistical computing framework.

In addition, this repository maintains code necessary to apply the method to other zooarchaeological measurement datasets, including code for data wrangling if faunal datasets are maintained following the [OpenContext EOL Computational Data Challenge standards](https://opencontext.org/tables/f07bce4fb08cfe926505c9e534d89a09).

## How to apply this method to your data

To apply the method to your own zooarchaeological data, you first need to set up your data in the following format. I have typically used the {data.table} package to prepare my data, which is what the "data_cleaning.R" script uses. However, the following steps can also be done using a different package or even outside of R, as long as the following requirements are met:

|Arch_ID|Site_No|Specimen_No|Element|Element_Portion|Measurement|Measurement_Set|Immature|Value_obs|Value_sd|Reference_obs|Reference_sd|
|---|---|---|---|---|---|---|---|---|---|---|---|
|SITE_059|1|1|Scapula|1|Sca_GLP|1|1|34.5|0.345|33.0|0.330|
|SITE_003|1|2|Scapula|1|Sca_GLP|1|1|30.0|0.300|33.0|0.330|
|SITE_846|1|3|Prox. Radius|2|Rad_Bp|2|1|28.0|0.280|33.5|0.335|
|SITE_841|1|4|Dist. Radius|3|Rad_Bp|2|0|35.0|0.350|33.5|0.335|
|SITE_841|1|4|Dist. Radius|3|Rad_Bd|3|0|32.1|0.321|32.0|0.320|

This is a dummy example using some fake data of sheep measurement data.
- The variable "ARCH_ID" is a clarifying label not necessary for the model, but is the specimen's original record number (and is helpful to relate the results back to other data). Other relevant specimen-level data can also be included in the table but will not be directly included in the model.
- "Site_No" is required but isn't called explicitly if only a single assemblage is analyzed. In such a case, a value of 1 is still necessary for every observation.
- "Specimen_No" is required; the model uses this variable as the baseline of the observations (i.e., every unique specimen is used as the observations). Specimens can have multiple observations (as in the radius with "Specimen_No" 4) that inform the modeled LSI value.
- "Element" is a clarifying label and thus not necessarily required; it is the label of the element portion and is helpful to interpret results
- "Element_Portion" is required; the model uses this numeric code to group observations into different categories for the multilevel model
  - It is important to ensure that all specimens have a single "Element_Portion" (even if, as in the example, both the proximal and distal parts of the element have observed measurements).
- "Measurement" and "Measurement_Set" have the same general structure as "Element" and "Element_Portion", though every type of measurement should have its own associated "Measurement_Set" (see the values for the two measurements in "Specimen_No" 4). Note that "Measurement" is a clarifying label and is thus not strictly required.
- "Immature" is required; this is a binary (0/1) variable about whether the specimen can be *potentially* immature based on its fusion status or biology.
  - Note that this is a specimen-level variable; all observations of the same specimen must have the same value
- "Value_obs" is the actual observed measurement value
  - "Value_sd" is the measurement error of that observation (in the paper, I define this as 1% of the observed value)
- "Reference_obs" is the relevant observed measurement from the standard animal (used to calculate LSI values)
  - "Reference_sd" is the meausrement error of that observation, following the same guidelines

These variables represent the minimal set of variables necessary to apply a Bayesian multilevel mixture model to a dataset (minus the three variables noted as clarifying labels). Note that the application code, as written, expects these particular variable names. Again, the "data_cleaning.R" script can help create the numeric codes (Specimen_No, Element_Portion, Measurement_Set) that the Stan model expects. These are done automatically by translating the clarifying labels into factors, which means the numeric codes follow alphabetical order. However, you can set up your own order if you want the numeric codes to be in a different order (e.g., one that mirrors anatomical position).

### Running the model (single assemblage)

Once you've cleaned your data, you can apply the model to a single assemblage using the "single_assemblage_mixmod.R" script. Note that if you are downloading single files from this repository, you will also need the "LSI_mixture_model_singlesite.stan" file (you can also download the "LSI_mixture_model_singlesite.exe" file but you can compile that on your own computer using the .stan file). The .R file lists the necessary packages. I recommend looking at the following guides to [install cmdstanr](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) and [install rstan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

Ensure that the .stan file is in the same place as your working directory when using the model. I recommend [creating an R project](https://r4ds.had.co.nz/workflow-projects.html) and placing the .stan file in a subfolder ("Scripts", for example) to standardize the relationship between your working directory and the location of the file across different applications.

### Running the model (multiple assemblages)

To apply mixture models to multiple assemblages at once (to, say, look at biometric variation across a region), use the "multisite_mixmod.R" script. As in the single assemblage example, you will also need the "LSI_mixture_model_multisite.stan" file. See the links and the discussion in the previous section for recommendations on installing [cmdstanr](https://mc-stan.org/cmdstanr/articles/cmdstanr.html), [rstan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started), and using [R projects](https://r4ds.had.co.nz/workflow-projects.html) for your analyses.

### You've applied the model, now what? (posterior analysis)

The file "posterior_analysis.R" contains the R functions used to create posterior simulations of the assemblage's composition (and/or fusion rates), as shown in the manuscript.

This also includes functions to produce the (archaeological) composition plots as seen in the manuscript. The necessary packages ({data.table}, {ggplot2}, {ggdist}, {ggpubr}) are included at the top of the script; include those in your own scripts if copy/pasting the functions into your own scripts (which is definitely allowed!).

Note that the script for simulating the composition assemblages uses multi-core processing to speed things up. This is not necessary but will make your life easier. If using those, include the {doParallel} and {parallel} packages.
