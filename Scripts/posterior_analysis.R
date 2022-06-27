####Posterior Analysis####
#This script is part of the ZooarchMixMod GitHub project (https://github.com/wolfhagenj/ZooarchMixMod)
#It provides instructions to take Bayesian multilevel mixture model results (posterior samples) and use
#them to simulate assemblages of known-identity specimens (immature, adult-sized female, and adult-sized males).
#These simulated assemblages can then be used to estimate assemblage composition, fusion rates, or other derived
#quantities. Further, these methods work for measured, modeled, or full assemblages.

#This code is meant to be copied and adapted for your own personal use. While it is designed for data
#set up like the OpenContext data used in the original manuscript, the same principles can be adapted
#to data set up in a different way. The principles can even be used to standardize the data yourself
#outside of R. See the ReadMe for the GitHub repository to see what the model expects in terms of
#data setup.

####Load the packages####
#These packages are necessary for the script. Use the following line (without the comment #) to install
#the packages if you have not already installed them:
#install.packages("data.table", "parallel", "doParallel", "ggdist", "ggplot2", "ggpubr")
library("data.table")
library("parallel")
library("doParallel")
library("ggdist")
library("ggplot2")
library("ggpubr")

####Functions####
#These functions are used in this script to simulate true identities of an assemblage and then estimate
#composition (or fusion rates) in each posterior sample.
#The ifm_sample() function is central to all of these, as it does the actual sampling of identities.
#The theta_finder() function is used to identify the correct set of probabilities for a specimen.

#Simulate a specimen as Immature/Female/Male based on probabilities
ifm_sample <- function(p_immature, p_female, p_male) {
  sample(c("Immature", "Female", "Male"), 1, replace = T, prob = c(p_immature, p_female, p_male))
}

#Identifies the correct set of Immature/Female/Male probabilities for a specimen based on whether it was directly modeled
#and the element portion that it is
theta_finder <- function(Specimen_No, Element_Portion, Immature, Immature_Proportions, element_thetas, specimen_probs) {
  thetas <- rep(NA, 3)
  immature_proportion <- Immature_Proportions[Portion %in% Element_Portion, Immature_Proportion]
  thetas[1:3] <- ifelse(rep(!is.na(Specimen_No), 3),
                        specimen_probs[Specimen_No, 1:3],
                        ifelse(rep(Immature == 1, 3),
                               c(element_thetas[1, Element_Portion] / immature_proportion,
                                 element_thetas[2, Element_Portion] * immature_proportion,
                                 element_thetas[3, Element_Portion] * immature_proportion),
                               c(0,
                                 element_thetas[2, Element_Portion] / (1 - min(c(element_thetas[1, Element_Portion], 0.99999))), #fail-safe for extreme values
                                 element_thetas[3, Element_Portion] / (1 - min(c(element_thetas[1, Element_Portion], 0.99999)))))
  )
  list(p_immature = thetas[1] / sum(thetas), p_female = thetas[2] / sum(thetas), p_male = thetas[3] / sum(thetas))
}

#Simulates the composition of an assemblage with Immature/Female/Male probabilities
ifm_composition <- function(measured_specimens) {
  #simulate sex assignments
  measured_specimens[, Simulated_Group := sapply(1:.N, function(x) ifm_sample(p_immature[x], p_female[x], p_male[x]))]
  #counts of the different elements by group (make sure there's a 0 if a value is missing)
  melt(measured_specimens[, .(`Immature` = sum(Simulated_Group %in% "Immature"),
                              `Female` = sum(Simulated_Group %in% "Female"),
                              `Male` = sum(Simulated_Group %in% "Male")), .(Site_No, Element_Portion, Element)], id.vars = c("Site_No", "Element_Portion", "Element"), variable.name = "Simulated_Group", value.name = "N")[order(Site_No, Element_Portion, Simulated_Group)]
}

#Simulates the fusion rate of an assemblage with Immature/Female/Male probabilities
ifm_fusion_rate <- function(measured_specimens) {
  #simulate sex assignments
  measured_specimens[, Simulated_Group := sapply(1:.N, function(x) ifm_sample(p_immature[x], p_female[x], p_male[x]))]
  #counts of the different elements by group and fusion status
  
  melt(measured_specimens[, .(`Immature` = sum(Simulated_Group %in% "Immature"),
                              `Female` = sum(Simulated_Group %in% "Female"),
                              `Male` = sum(Simulated_Group %in% "Male")), .(Site_No, Fusion_Element, Stage, Fusion)], id.vars = c("Site_No", "Fusion_Element", "Stage", "Fusion"), variable.name = "Simulated_Group", value.name = "N")[order(Site_No, Stage, Fusion_Element, Fusion, Simulated_Group)]
}

#function to get integer y-axis values only (source: https://gist.github.com/jhrcook/eb7b63cc57c683a6eb4986c4107a88ec)
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

#Plots the distributions of group-specific composition estimates
composition_analysis <- function(composition_data, site_no = 1, site_name = "Site Measured Assemblage", full_element_name_key = site_element_name_key) {
  ggplot(composition_data[Site_No %in% site_no, .(N = sum(N)), .(Iteration, Simulated_Group, Element, Element_Portion)][, .(Element, Element_Portion, group = Simulated_Group, Plot_Group = paste(Element, as.numeric(Simulated_Group), sep = "."), N)]) + aes(y = N, x = Plot_Group) +
    stat_slab(normalize = "groups", slab_type = "histogram", breaks = 0:(composition_data[Site_No %in% site_no, .(N = sum(N)), .(Iteration, Simulated_Group, Element, Element_Portion)][, max(N)]), aes(fill = group, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(0.80, 0.95, 1), labels = scales::percent_format())))) +
    stat_pointinterval(.width = c(0.8, 0.95), aes(color = group), position = position_dodge(width = 0.2, preserve = "single")) +
    scale_fill_ramp_discrete(name = "Interval", range = c(0.5, 0.15), na.translate = F) +
    scale_color_manual(name = "Group", values = c("black", "blue", "red"), na.translate = F) +
    scale_fill_manual(name = "Group", values = c("black", "blue", "red"), na.translate = F) +
    scale_x_discrete(name = "Element", breaks = composition_data[, .N, .(Element, Element_Portion)][order(Element_Portion), paste(Element, ".2", sep = "")], labels = composition_data[, .N, .(Element, Element_Portion)][full_element_name_key, on = "Element"][!is.na(Element_Portion)][order(Element_Portion), Full_Element_Name]) +
    scale_y_continuous(name = "Estimated Count", breaks = integer_breaks()) + coord_cartesian(expand = FALSE) +
    labs(title = site_name) +
    theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(size = 8, angle = -90, hjust = 0), axis.title = element_text(size = 12), axis.text.y = element_text(size = 8), axis.title.x = element_blank())
}

####Prepare the assemblages for composition analysis####
#Uses the results of the single_assemblage_mixmod script (and data_cleaning script) to estimate
#composition for the measured, modeled, and full assemblages.

#First, calculate new sets of element_portion-specific mixing probabilities for unmodeled element portions (for full assemblage)
site_unmodeled_element_portion_theta <- data.table(Iteration = rep(1:4000, full_assemblage[, .N, .(Element, Element_Portion)][Element_Portion > measured_assemblage[, .N, Element_Portion][, .N]][order(Element_Portion), .N]),
                                                        Element_Portion = rep(full_assemblage[, .N, .(Element, Element_Portion)][Element_Portion > measured_assemblage[, .N, Element_Portion][, .N]][order(Element_Portion), Element_Portion], each = 4000),
                                                        Element = rep(full_assemblage[, .N, .(Element, Element_Portion)][Element_Portion > measured_assemblage[, .N, Element_Portion][, .N]][order(Element_Portion), Element], each = 4000))
#specify hyper-parameter values for the model
for(i in 1:full_assemblage[, .N, .(Element, Element_Portion)][Element_Portion > measured_assemblage[, .N, Element_Portion][, .N]][order(Element_Portion), .N]) {
  element_theta_raw <- t(sapply(1:4000, function(x) MASS::mvrnorm(1, mu = single_assemblage_post$grand_theta_raw[x, 1:2], Sigma = diag(single_assemblage_post$element_sigma[x, 1:2]) %*% single_assemblage_post$Rho_element[x, 1:2, 1:2] %*% diag(single_assemblage_post$element_sigma[x, 1:2]))))
  site_unmodeled_element_portion_theta[(4000 * (i - 1) + 1):(4000 * i), c("theta_raw_1", "theta_raw_2") := .(element_theta_raw[, 1], element_theta_raw[, 2])]
}
#perform stick-breaking procedure to turn theta_raw into theta values
site_unmodeled_element_portion_theta[, theta_1 := boot::inv.logit(theta_raw_1 + log(0.5))]
site_unmodeled_element_portion_theta[, theta_2 := (1 - theta_1) * boot::inv.logit(theta_raw_2 + log(1))]
site_unmodeled_element_portion_theta[, theta_3 := (1 - (theta_1 + theta_2))]

#Second, calculate the immature proportions for the three datasets (affects how heavily potentially-immature specimens are weighted towards that category)
measured_immature_proportions <- measured_assemblage[, .(Immature_Proportion = mean(Immature)), .(Portion = Element_Portion)][order(Portion)]
modeled_immature_proportions <- modeled_assemblage[, .(Immature_Proportion = mean(Immature)), .(Portion = Element_Portion)][order(Portion)]
full_immature_proportions <- full_assemblage[, .(Immature_Proportion = mean(Immature)), .(Portion = Element_Portion)][order(Portion)]

#Third, this provides display names for the elements in the assemblages. Adjust as necessary:
#"Element" is the label for each Element_Portion
#"Full_Element_Name" is the name displayed in the composition plot
site_element_name_key <- data.table(
  Element = c("Ast", "Cal", "Sca", "Hum", "Mtc_prox", "Mtc_dist", "Mtt_dist", "Mtt_prox", "PH1", "Rad_dist", "Fem_prox", "Fem_dist", "Tib_dist", "Tib_prox", "Rad_prox", "Pel", "PH2", "Uln"),
  Full_Element_Name = c("Astragalus", "Calcaneus", "Scapula", "Humerus", "P. Metacarpal", "D. Metacarpal", "D. Metatarsal", "P. Metatarsal", "Proximal Phalanx", "D. Radius", "P. Femur", "D. Femur", "D. Tibia", "P. Tibia", "P. Radius", "Pelvis", "Middle Phalanx", "Ulna")
)

####Simulate the composition####
#Ideally, you will use parallel computing to do this much faster (using the doParallel and parallel packages)
#To do this, you need to assign different cores on your computer to a cluster (this is necessary for Windows, but should work on all platforms)
#setting resources to do simulations with multicore processing (on Windows)
num_cores <- detectCores(logical = T)
cl <- makeCluster(num_cores - 4) #creating a virtual cluster
registerDoParallel(cl) #registering the cluster

#Then export the objects, packages, and functions to the virtual cluster
clusterExport(cl, list('data.table', 'melt', 'ifm_sample', 'theta_finder', 'ifm_composition', 'measured_assemblage', 'measured_immature_proportions', 'modeled_assemblage', 'modeled_immature_proportions', 'full_assemblage', 'full_immature_proportions', 'single_assemblage_post', 'site_unmodeled_element_portion_theta'))

#Finally, run the simulations. This may take time, you can use a package like beepr::beep() to have the computer play a chime when it is done
measured_assemblage_simulated_composition <- rbindlist(parLapply(cl = cl, 1:4000, fun = function(x) ifm_composition(measured_assemblage[, theta_finder(Specimen_No, Element_Portion, Immature, measured_immature_proportions, element_thetas = single_assemblage_post$theta[x, , ], specimen_probs = single_assemblage_post$specimen_prob[x, , ]), .(ID, Specimen_No, Element, Element_Portion, Immature)][, .(ID, Specimen_No, Site_No = 1, Element, Element_Portion, Immature, p_immature, p_female, p_male)])[, .(Iteration = x, Site_No, Element_Portion, Element, Simulated_Group, N)]))
modeled_assemblage_simulated_composition <- rbindlist(parLapply(cl = cl, 1:4000, fun = function(x) ifm_composition(modeled_assemblage[, theta_finder(Specimen_No, Element_Portion, Immature, modeled_immature_proportions, element_thetas = single_assemblage_post$theta[x, , ], specimen_probs = single_assemblage_post$specimen_prob[x, , ]), .(ID, Specimen_No, Element, Element_Portion, Immature)][, .(ID, Specimen_No, Site_No = 1, Element, Element_Portion, Immature, p_immature, p_female, p_male)])[, .(Iteration = x, Site_No, Element_Portion, Element, Simulated_Group, N)]))
full_assemblage_simulated_composition <- rbindlist(parLapply(cl = cl, 1:4000, fun = function(x) ifm_composition(full_assemblage[, theta_finder(Specimen_No, Element_Portion, Immature, full_immature_proportions, element_thetas = cbind(single_assemblage_post$theta[x, , ], t(site_unmodeled_element_portion_theta[Iteration %in% x, .(theta_1, theta_2, theta_3)])), specimen_probs = single_assemblage_post$specimen_prob[x, , ]), .(ID, Specimen_No, Element, Element_Portion, Immature)][, .(ID, Specimen_No, Site_No = 1, Element, Element_Portion, Immature, p_immature, p_female, p_male)])[, .(Iteration = x, Site_No, Element_Portion, Element, Simulated_Group, N)]))

####Run the composition analyses####
#This saves each plot as an object, which can then be called to display the plot
#or combined together using the ggarrange() function
site_measured_composition <- composition_analysis(measured_assemblage_simulated_composition, site_no = 1, site_name = "Measured Assemblage")
site_modeled_composition <- composition_analysis(modeled_assemblage_simulated_composition, site_no = 1, site_name = "Modeled Assemblage")
site_full_composition <- composition_analysis(full_assemblage_simulated_composition, site_no = 1, site_name = "Full Assemblage")
#
site_measured_composition

#Putting the figures together into a single figure
site_composition_figure <- ggarrange(site_measured_composition,
                                     site_modeled_composition,
                                     site_full_composition,
                                     nrow = 3, ncol = 1, common.legend = T, legend = "bottom")
#
final_figure <- annotate_figure(site_composition_figure, top = text_grob("Estimated Composition of Site Cattle", face = "bold", size = 14))
#
final_figure
