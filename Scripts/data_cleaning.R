####Data Cleaning Script####
#This script is part of the ZooarchMixMod GitHub project (https://github.com/wolfhagenj/ZooarchMixMod)
#It provides instructions to standardize zooarchaeological measurement data so that it can be used by
#the Bayesian multilevel mixture models described in the project.

#This code is meant to be copied and adapted for your own personal use. While it is designed for data
#set up like the OpenContext data used in the original manuscript, the same principles can be adapted
#to data set up in a different way. The principles can even be used to standardize the data yourself
#outside of R. This script helps with creating the necessary numeric codes that the Stan model expects.

####Load the packages####
#These packages are necessary for the script. Use the following line (without the comment #) to install
#the packages if you have not already installed them:
#install.packages("data.table", "zoolog")
library("data.table") #this loads the package so you can use the functions
library("zoolog")

####Bring in your data with fread()####
#These scripts are structured around you setting up your workflow in an R project.
#See https://r4ds.had.co.nz/workflow-projects.html for more information about R projects and their use.
#Importantly, the script assumes that your working directory has a sub-folder called "Data" where your
#dataset is located.
dataset <- fread("./Data/Barcin Hoyuk Zooarchaeology Data (OpenContext - DOI 10.678-M7MS3QN7).csv")
#this is the Barcin Hoyuk dataset from OpenContext; replace with your own file name

#List the Element labels you want to find in your data
modeled_elements <- c("scapula", "humerus", "radius bone", "fused metacarpal bones 3 and 4",
                      "femur", "tibia", "talus", "calcaneus", "fused metatarsal bones 3 and 4",
                      "proximal phalanx", "middle phalanx")
#these names follow the Encyclopedia of Life as part of the EOL Computational Challenge

#Sometimes other element names may be included in your specific dataset, include those too
extra_element_names <- c("fused metacarpal bones 3 and 4; metacarpal bone of digit 3; metacarpal bone of digit 4; metacarpal bone of digit 5",
                         "metatarsal bone of digit 5; metatarsal bone of digit 2; fused metatarsal bones 3 and 4")

#Include the taxonomic names in your dataset; again, these follow the EOL standards in this example
modeled_taxon <- c("Bos taurus Linnaeus, 1758", "Bos primigenius")

#Collect the relevant observations with all relevant measurements
#data.table uses the following syntax
#datatablename[filter, .(columns)]
#See https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html for an introduciton to data.table syntax
measurement_data <- dataset[`Has Biological Taxonomy [Label]` %in% c(modeled_taxon) & `Has anatomical identification [Label]` %in% c(modeled_elements, extra_element_names),
                                  .(ID = paste("Barcin", `Label`, sep = " "),
                                    #I have separated these columns into separate lines to make it easier to read, but this is not necessary
                                    #NOTE: if running a multisite analysis, I find it helpful to append the site name to whatever label the
                                    #specimen has to ensure that there are no repeated labels across different sites
                                    Site = "Barcin",
                                    Taxon = `Has Biological Taxonomy [Label]`,
                                    Anatomy = `Has anatomical identification [Label]`,
                                    `Proximal Fusion` = `Has fusion character [Proximal Label]`,
                                    `Distal Fusion` = `Has fusion character [Distal Label]`,
                                    GLP = `Scap_glp`,
                                    #some datasets may have conditional variables you need to collapse together
                                    Bd = ifelse(`Has anatomical identification [Label]` %in% "talus", Tal_bd, Bd),
                                    #if one of the measurement types you're interested in is missing from the assemblage completely, use NA
                                    #this provides consistency across your datasets (and is necessary when doing multisite analysis, so that
                                    #you have the same column names in each dataset)
                                    BT = NA,
                                    Bp = Bp,
                                    BFp = Bfp,
                                    DC = NA,
                                    GB = Gb)][!is.na(GLP) | !is.na(Bd) | !is.na(BT) | !is.na(Bp) | !is.na(BFp) | !is.na(DC) | !is.na(GB)]
#this last line is another filter, keeping only those specimens that have at least one of the potential measurements

#You may need to troubleshoot your dataset, removing specimens with anomalous measurements or correcting typos
#It's best practice to record those steps within your analytical script rather than changing the data directly unless you are the data producer
#The := command is unique to data.table, allowing you to create new variables (which we'll do in other cases) or redefine values (as in here)
measurement_data[ID %in% c("Barcin Bone 448", "Barcin Bone 358"), c("Bp") := NA] #anomalously small Bp measurements (17.2, 18.1); error in taxon coding?


#Create a key between the "Anatomy" label from the original dataset and the "Element" label that will be used in the model
site_element_key <- data.table(`Anatomy Label` = c("scapula", "humerus", "radius bone", "fused metacarpal bones 3 and 4; metacarpal bone of digit 3; metacarpal bone of digit 4; metacarpal bone of digit 5", "femur", "tibia", "talus", "calcaneus", "metatarsal bone of digit 5; metatarsal bone of digit 2; fused metatarsal bones 3 and 4", "proximal phalanx", "middle phalanx"),
                                 `Element Label` = c("Sca", "Hum", "Rad", "Mtc", "Fem", "Tib", "Ast", "Cal", "Mtt", "PH1", "PH2"))

####Get the reference measurements####
#The "zoolog" package has a lot of standard measurements to calculate LSI values. In this example, I use the "Bos primigenius"
#standard animal, though I also add a few measurement sets for elements that aren't included
bos_standard_animal <- data.table(referencesDatabase$`Bos primigenius`$Degerbol)
#add in the measurement of the Scapula GLP (89.0 mm), Calcaneus GB (46.0 mm); Degerbol 1970: Table 13, Table 18
bos_standard_animal <- rbind(bos_standard_animal, data.table(TAX = rep("Bos primigenius", 2), EL = c("Scapula", "Calcaneus"), Measure = c("GLP", "GB"), Standard = c(89.0, 46)))

#Create "Element" and "Measurement" variables--specify these for your own analyses if you're not using these elements
bos_standard_animal[EL %in% "Scapula", c("Element", "Measurement") := .("Sca", paste("Sca", Measure, sep = "_"))]
bos_standard_animal[EL %in% "Humerus", c("Element", "Measurement") := .("Hum", paste("Hum", Measure, sep = "_"))]
bos_standard_animal[EL %in% "Radius", c("Element", "Measurement") := .("Rad", paste("Rad", Measure, sep = "_"))]
bos_standard_animal[EL %in% "Metacarpus", c("Element", "Measurement") := .("Mtc", paste("Mtc", Measure, sep = "_"))]
bos_standard_animal[EL %in% "Femur", c("Element", "Measurement") := .("Fem", paste("Fem", Measure, sep = "_"))]
bos_standard_animal[EL %in% "Tibia", c("Element", "Measurement") := .("Tib", paste("Tib", Measure, sep = "_"))]
bos_standard_animal[EL %in% "Astragalus", c("Element", "Measurement") := .("Ast", paste("Ast", Measure, sep = "_"))]
bos_standard_animal[EL %in% "Calcaneus", c("Element", "Measurement") := .("Cal", paste("Cal", Measure, sep = "_"))]
bos_standard_animal[EL %in% "Metatarsus", c("Element", "Measurement") := .("Mtt", paste("Mtt", Measure, sep = "_"))]
bos_standard_animal[EL %in% "Phalanx 1 ant.", c("Element", "Measurement") := .("PH1", paste("PH1", Measure, sep = "_"))]
bos_standard_animal[EL %in% "Phalanx 2 ant.", c("Element", "Measurement") := .("PH1", paste("PH2", Measure, sep = "_"))]

####Use a function to restructure the data####
#The dataset as currently constructed has every specimen as a single row (with potentially multiple measurements)
#We want every observed measurement to be a single row while maintaining labels that tie measurements from the same specimen together
#This function restructures the data into the format expected for the Bayesian multilevel mixture model
#(Sidenote: this same logic would underpin tying all measurements together from a single individual, in the case of multiple articulated elements)

assemblage_restructure <- function(dataset, reference_dataset, element_key, included_measurements = c("Sca_GLP", "Hum_Bd", "Hum_BT", "Rad_Bp", "Rad_BFp", "Rad_Bd", "Mtc_Bp", "Mtc_Bd", "Fem_DC", "Fem_Bd", "Tib_Bp", "Tib_Bd", "Ast_Bd", "Cal_GB", "Mtt_Bp", "Mtt_Bd", "PH1_Bp", "PH2_Bp")) {
  #Element: identify what end of the bone (proximal/distal) is included in the measurements (for element portions).
  #goes through element_key to change anatomical part names to the new labels
  #Note: key may have multiple Anatomy Labels that fit onto the same `Element Label`
  for(i in 1:nrow(element_key)) {
    dataset[Anatomy %in% element_key[i, `Anatomy Label`], Element := element_key[i, `Element Label`]]
    #Separate out proximal and distal element portions based on presence of later-fusing measurements (e.g., radius bone is proximal/distal based on presence of Bd measurement)
    if(element_key[i, `Element Label`] %in% c("Rad", "Mtc", "Fem", "Mtt")) {
      dataset[Anatomy %in% element_key[i, `Anatomy Label`], Element := ifelse(!is.na(Bd), paste(element_key[i, `Element Label`], "dist", sep = "_"), paste(element_key[i, `Element Label`], "prox", sep = "_"))]
    }
    if(element_key[i, `Element Label`] %in% c("Tib")) {
      dataset[Anatomy %in% element_key[i, `Anatomy Label`], Element := ifelse(!is.na(Bp), paste(element_key[i, `Element Label`], "prox", sep = "_"), paste(element_key[i, `Element Label`], "dist", sep = "_"))]
    }
  }
  #Immature: identify whether specimens COULD be immature based on element portion and fusion status, if relevant
  dataset[Element %in% c("Sca", "Ast"), Immature := 1] #does not fuse / subject to post-fusion growth and has no relevant fusion data to exclude immature status
  dataset[Element %in% c("Hum", "Cal", "PH1", "PH2"), Immature := as.numeric(`Proximal Fusion` %in% c("Proximal epiphysis fused", "Proximal epiphysis fusing") == F)]
  dataset[Element %in% c("Rad_prox", "Rad_dist", "Mtc_prox", "Mtc_dist", "Mtt_prox", "Mtt_dist"), Immature := as.numeric(`Distal Fusion` %in% c("Distal epiphysis fused", "Distal epiphysis fusing") == F)]
  dataset[Element %in% c("Fem_prox", "Fem_dist", "Tib_prox", "Tib_dist"), Immature := as.numeric((`Proximal Fusion` %in% c("Proximal epiphysis fused", "Proximal epiphysis fusing") | `Distal Fusion` %in% c("Distal epiphysis fused", "Distal epiphysis fusing")) == F)]
  #transform into long-form data
  dataset_longdata <- melt(dataset, id.vars = c("ID", "Site", "Anatomy", "Element", "Proximal Fusion", "Distal Fusion", "Immature"), measure.vars = c("GLP", "Bd", "BT", "Bp", "BFp", "DC", "GB"), variable.name = "Measure", value.name = "Measurement_value")[!is.na(Measurement_value)]
  #Measurement: combine element and measurement names (for measurement sets)
  #goes through element_key to change anatomical part names to the new labels
  #Note: key may have mlutiple Anatomy Labels that fit onto the same `Element Label`
  for(i in 1:nrow(element_key)) {
    dataset_longdata[Anatomy %in% element_key[i, `Anatomy Label`], Measurement := paste(element_key[i, `Element Label`], Measure, sep = "_")]
  }
  #bring in the reference measurement for each value (joined by measurement) and limiting the assemblage to the measurements that we're interested in modeling
  mixmod_data <- dataset_longdata[reference_dataset[, .(Reference_value = Standard, Measurement)], on = c("Measurement")][!is.na(ID) & Measurement %in% included_measurements]
  #Create numeric labels for measurement sets, element portions, and specimens (for Stan modeling)
  mixmod_data[, Dimension := as.numeric(as.factor(Measurement))]
  mixmod_data[, Element_Portion := as.numeric(as.factor(Element))]
  mixmod_data[, Specimen_No := as.numeric(as.factor(ID))]
  #return new dataset
  mixmod_data
}

#List of measurements we want to use in our analysis (what to keep from the dataset)
#This can be amended to focus on your particular analysis
measurement_list <- c("Sca_GLP", "Hum_Bd", "Hum_BT", "Rad_Bp", "Rad_BFp", "Rad_Bd", "Mtc_Bp", "Mtc_Bd", "Fem_DC", "Fem_Bd", "Tib_Bp", "Tib_Bd", "Ast_Bd", "Cal_GB", "Mtt_Bp", "Mtt_Bd", "PH1_Bp", "PH2_Bp")

#Running the function to create the dataset standardized for the Bayesian multilevel mixture model
#Creates the numeric codes necessary--these are, by default, based on alphabetization, but can be
#changed in the function or afterwards to reflect different orderings of the categories (they're arbitrary)
site_mixmod_data <- assemblage_restructure(dataset = measurement_data,
                                             reference_dataset = bos_standard_animal,
                                             element_key = site_element_key,
                                             included_measurements = c(measurement_list))

#You may get a warning that not all of your 'measure.vars' are the same type.
#This happens when you have completely missing measurement sets (e.g., all of the "DC" variable is NA).
#Double-check that the mixmod_data looks correct, then the warning can be safely ignored

####Collecting demographic observations####
#The model uses demographic observations to estimate the proportion of immature animals and the adult sex ratio
#In the paper, this is focused on phalanx fusion (for immature animals) and pelvis morphology (for sex ratios)
#NOTE: for sheep or goat models, I tend to combine both the taxon and "sheep/goat" to account for less identifiability
#among unfused elements
#These are modeling choices, however, so can be changed depending on the context.

site_demographic_observations <- dataset[`Has Biological Taxonomy [Label]` %in% c(modeled_taxon), .(Site = "Barcin Hoyuk",
                                                                                                     N_Unfused = sum(`Has anatomical identification [Label]` %in% c("proximal phalanx", "middle phalanx") & `Has fusion character [Proximal Label]` %in% "Proximal epiphysis unfused"),
                                                                                                     N_Ageable = sum(`Has anatomical identification [Label]` %in% c("proximal phalanx", "middle phalanx") & `Has fusion character [Proximal Label]` %in% "" == F),
                                                                                                     N_Female = sum(`Has anatomical identification [Label]` %in% "innominate bone" & `Has physiological sex determination [Label]` %in% "female organism"),
                                                                                                     N_Sexable = sum(`Has anatomical identification [Label]` %in% "innominate bone" & `Has physiological sex determination [Label]` %in% c("female organism", "male organism")))]
#the "Site" variable isn't necessary, but helps keep the data clear, especially when bringing multiple sites together

####Collecting the "modeled" assemblage####
#As noted in Section 2.5 of the manuscript, a "measurement assemblage" represents a sample from the "modeled" assemblage
#(all recovered specimens from those element portions), assuming that "measurability" is random. We are interested in
#the composition of the measured assemblage as a way of informing us about the composition of the modeled assemblage.
#The model parameters thus relate to the modeled assemblage; we can use those parameters to estimate the composition
#but first need to collect the entire "modeled assemblage" to make composition (or fusion) plots

#IMPORTANT: the key thing is to make sure that the measurable remains get their original Specimen_No values
#used by the model, this ensures that you get the correct prob_specimen estimates.
#It is also necessary to specify the relevant Element_Portion and Immature variables for measured and non-measured specimens

#For convenience, we can also create a specific "measured_assemblage" that has only a single row per specimen
#but also preserves the Specimen_No, Element_Portion, and Immature codes from the site_mixmod_data
measured_assemblage <- site_mixmod_data[, .N, .(ID, Site, Specimen_No, Element, Element_Portion, Immature)]
#here, the "N" variable is the number of observed measurements the specimen has

#Create the modeled assemblage
#if data are set up in the EOL format, this is simply done by removing measurement variables from the original
#line that created the measurement assemblage
modeled_assemblage <- dataset[`Has Biological Taxonomy [Label]` %in% c(modeled_taxon) & `Has anatomical identification [Label]` %in% c(modeled_elements, extra_element_names),
                              .(ID = paste("Barcin", `Label`, sep = " "),
                                #I have separated these columns into separate lines to make it easier to read, but this is not necessary
                                #NOTE: if running a multisite analysis, I find it helpful to append the site name to whatever label the
                                #specimen has to ensure that there are no repeated labels across different sites
                                Site = "Barcin",
                                Taxon = `Has Biological Taxonomy [Label]`,
                                Anatomy = `Has anatomical identification [Label]`,
                                `Proximal Fusion` = `Has fusion character [Proximal Label]`,
                                `Distal Fusion` = `Has fusion character [Distal Label]`)]

#The next step is to align the "Specimen_No" variable for the specimens with measurements
#and leave it as NA (missing value) for non-measured specimens.
#This is done by using a "join"
modeled_assemblage <- measured_assemblage[, .(ID, Specimen_No, Element_Portion, Element, Immature)][modeled_assemblage, on = "ID"]

#The final step is to assign the Element_Portion and Immature variables for all non-measured specimens.
#These must be the same codes as used in the original model (as they impact what parameter values are used)
#These codes are part of the assemblage_restructure() function
for(i in 1:site_element_key[, .N]) {
  modeled_assemblage[is.na(Element_Portion) & Anatomy %in% site_element_key[i, `Anatomy Label`], Element := site_element_key[i, `Element Label`]]
  #Separate out proximal and distal element portions based on presence of later-fusing measurements (e.g., radius bone is proximal/distal based on presence of Bd measurement)
  if(site_element_key[i, `Element Label`] %in% c("Rad", "Mtc", "Fem", "Mtt")) {
    modeled_assemblage[is.na(Element_Portion) & Anatomy %in% site_element_key[i, `Anatomy Label`], Element := ifelse(`Distal Fusion` %in% "" == F, paste(site_element_key[i, `Element Label`], "dist", sep = "_"), paste(site_element_key[i, `Element Label`], "prox", sep = "_"))]
  }
  if(site_element_key[i, `Element Label`] %in% c("Tib")) {
    modeled_assemblage[is.na(Element_Portion) & Anatomy %in% site_element_key[i, `Anatomy Label`], Element := ifelse(`Proximal Fusion` %in% "" == F, paste(site_element_key[i, `Element Label`], "prox", sep = "_"), paste(site_element_key[i, `Element Label`], "dist", sep = "_"))]
  }
}

#Troubleshooting: remove unobserved element portions from modeled assemblage
#This is in case there were element portions in the original modeled set
#that had no actual measurements (e.g., if there were no measured Distal Radius specimens).
#Since we don't have element portion-specific parameter estimates for that element,
#we cannot considered it part of the "modeled" assemblage

#This is again done by a join
modeled_assemblage <- measured_assemblage[, .N, .(Element_Portion, Element)][, .(Element_Portion, Element)][modeled_assemblage, on = "Element"][!is.na(Element_Portion)]

#Get the Immature variable (potentially immature status based on fusion/tooth eruption)
#This is only assigning Immature status for non-measured specimens
#Measured specimens had their Immature status drawn from the site_mixmod_data records
#(in case Immature values were changed through troubleshooting)
modeled_assemblage[is.na(Immature) & Element %in% c("Sca", "Ast"), Immature := 1] #does not fuse / subject to post-fusion growth and has no relevant fusion data to exclude immature status
modeled_assemblage[is.na(Immature) & Element %in% c("Hum", "Cal", "PH1", "PH2"), Immature := as.numeric(`Proximal Fusion` %in% c("Proximal epiphysis fused", "Proximal epiphysis fusing") == F)]
modeled_assemblage[is.na(Immature) & Element %in% c("Rad_prox", "Rad_dist", "Mtc_prox", "Mtc_dist", "Mtt_prox", "Mtt_dist"), Immature := as.numeric(`Distal Fusion` %in% c("Distal epiphysis fused", "Distal epiphysis fusing") == F)]
modeled_assemblage[is.na(Immature) & Element %in% c("Fem_prox", "Fem_dist", "Tib_prox", "Tib_dist"), Immature := as.numeric((`Proximal Fusion` %in% c("Proximal epiphysis fused", "Proximal epiphysis fusing") | `Distal Fusion` %in% c("Distal epiphysis fused", "Distal epiphysis fusing")) == F)]

####Collecting the "full" assemblage####
#The same logic that underpins the estimation of a "modeled" assemblage can be extended to element portions
#that had no observed measurements at all (a "full" assemblage). This is useful for extending an analysis
#to element portions (parts of the body that you want to study) that either do not have standard breadth
#measurements or did not have any observed measurements but you want to model them directly anyway.

#First, define the additional elements (and the associated "element portion" labels you want to use)
additional_elements <- c("innominate bone", "ulna") #same naming convention as the original elements
additional_element_labels <- c("Pel", "Uln") #IMPORTANT! Make sure that the order is the same across both

#Expand the element key to create new element labels
site_full_element_key <- rbind(site_element_key,
                                    data.table(`Anatomy Label` = additional_elements, `Element Label` = additional_element_labels))

#Second, create the full assemblage
#if data are set up in the EOL format, this is simply done by removing measurement variables from the original
#line that created the measurement assemblage and
#adding additional_elements to the filter on `Has anatomical identification [Label]`
full_assemblage <- dataset[`Has Biological Taxonomy [Label]` %in% c(modeled_taxon) & `Has anatomical identification [Label]` %in% c(modeled_elements, extra_element_names, additional_elements),
                           .(ID = paste("Barcin", `Label`, sep = " "),
                             #I have separated these columns into separate lines to make it easier to read, but this is not necessary
                             #NOTE: if running a multisite analysis, I find it helpful to append the site name to whatever label the
                             #specimen has to ensure that there are no repeated labels across different sites
                             Site = "Barcin",
                             Taxon = `Has Biological Taxonomy [Label]`,
                             Anatomy = `Has anatomical identification [Label]`,
                             `Proximal Fusion` = `Has fusion character [Proximal Label]`,
                             `Distal Fusion` = `Has fusion character [Distal Label]`)]

#The rest of the formatting is largely the same as with the modeled assemblage, though including the additional elements
#Some troubleshooting may be necessary to include the right mix of additional elements you want to include

#The next step is to align the "Specimen_No" variable for the specimens with measurements
#and leave it as NA (missing value) for non-measured specimens.
#This is done by using a "join"
full_assemblage <- measured_assemblage[, .(ID, Specimen_No, Element_Portion, Element, Immature)][full_assemblage, on = "ID"]

#The final step is to assign the Element_Portion and Immature variables for all non-measured specimens.
#These must be the same codes as used in the original model (as they impact what parameter values are used)
#These codes are part of the assemblage_restructure() function
for(i in 1:site_full_element_key[, .N]) {
  full_assemblage[is.na(Element_Portion) & Anatomy %in% site_full_element_key[i, `Anatomy Label`], Element := site_full_element_key[i, `Element Label`]]
  #Separate out proximal and distal element portions based on presence of later-fusing measurements (e.g., radius bone is proximal/distal based on presence of Bd measurement)
  if(site_full_element_key[i, `Element Label`] %in% c("Rad", "Mtc", "Fem", "Mtt")) {
    full_assemblage[is.na(Element_Portion) & Anatomy %in% site_full_element_key[i, `Anatomy Label`], Element := ifelse(`Distal Fusion` %in% "" == F, paste(site_full_element_key[i, `Element Label`], "dist", sep = "_"), paste(site_full_element_key[i, `Element Label`], "prox", sep = "_"))]
  }
  if(site_full_element_key[i, `Element Label`] %in% c("Tib")) {
    full_assemblage[is.na(Element_Portion) & Anatomy %in% site_full_element_key[i, `Anatomy Label`], Element := ifelse(`Proximal Fusion` %in% "" == F, paste(site_full_element_key[i, `Element Label`], "prox", sep = "_"), paste(site_full_element_key[i, `Element Label`], "dist", sep = "_"))]
  }
}

#Assign element portion numbers (matching original measured assemblage)
#Unlike the modeled assemblage, we don't want to remove element portions without any observed measurements.
#We do, however, want to keep the same Element_Portion labels for those element portions WITH observed measurements,
#then add new labels for the unobserved element portions

#This is again done by a join
full_assemblage <- measured_assemblage[, .N, .(Element_Portion, Element)][, .(Element_Portion, Element)][full_assemblage, on = "Element"]

#Defining the new Element_Portion labels will be done here alphabetically by "Element" label (among unobserved element portions)
unobserved_element_portions <- full_assemblage[is.na(Element_Portion), .N, Element][order(Element), Element]
for(i in 1:length(unobserved_element_portions)) {
  full_assemblage[Element %in% unobserved_element_portions[i], Element_Portion := (measured_assemblage[, .N, Element_Portion][, .N] + i)]
}

#Get the Immature variable (potentially immature status based on fusion/tooth eruption)
#This is only assigning Immature status for non-measured specimens
#Measured specimens had their Immature status drawn from the site_mixmod_data records
#(in case Immature values were changed through troubleshooting)
full_assemblage[is.na(Immature) & Element %in% c("Sca", "Ast"), Immature := 1] #does not fuse / subject to post-fusion growth and has no relevant fusion data to exclude immature status
full_assemblage[is.na(Immature) & Element %in% c("Hum", "Cal", "PH1", "PH2"), Immature := as.numeric(`Proximal Fusion` %in% c("Proximal epiphysis fused", "Proximal epiphysis fusing") == F)]
full_assemblage[is.na(Immature) & Element %in% c("Rad_prox", "Rad_dist", "Mtc_prox", "Mtc_dist", "Mtt_prox", "Mtt_dist"), Immature := as.numeric(`Distal Fusion` %in% c("Distal epiphysis fused", "Distal epiphysis fusing") == F)]
full_assemblage[is.na(Immature) & Element %in% c("Fem_prox", "Fem_dist", "Tib_prox", "Tib_dist"), Immature := as.numeric((`Proximal Fusion` %in% c("Proximal epiphysis fused", "Proximal epiphysis fusing") | `Distal Fusion` %in% c("Distal epiphysis fused", "Distal epiphysis fusing")) == F)]

#Codes for the additional elements (update for other lists of additional elements)
full_assemblage[is.na(Immature) & Element %in% c("Pel"), Immature := as.numeric((`Proximal Fusion` %in% c("Proximal epiphysis fused", "Proximal epiphysis fusing") | `Distal Fusion` %in% c("Distal epiphysis fused", "Distal epiphysis fusing")) == F)]
full_assemblage[is.na(Immature) & Element %in% c("Uln"), Immature := as.numeric((`Proximal Fusion` %in% c("Proximal epiphysis fused", "Proximal epiphysis fusing") | `Distal Fusion` %in% c("Distal epiphysis fused", "Distal epiphysis fusing")) == F)]

####Exporting the files####
#While the code has created objects that can be used within an R session, you may also want to save the files to upload later
#These lines of code will export the data (to a sub-folder called "Data") as .CSV files that can be loaded by other scripts or
#used in other programs (like Excel)
write.csv(site_mixmod_data, file = "./Data/Site Mixture Model Data.csv", row.names = FALSE)
write.csv(site_demographic_observations, file = "./Data/Site Demographic Observations.csv", row.names = FALSE)
write.csv(measured_assemblage, file = "./Data/Site Measured Assemblage Data.csv", row.names = FALSE)
write.csv(modeled_assemblage, file = "./Data/Site Modeled Assemblage Data.csv", row.names = FALSE)
write.csv(full_assemblage, file = "./Data/Site Full Assemblage Data.csv", row.names = FALSE)
