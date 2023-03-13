###Text mining algorithm version 2

######### load libraries and install limma package ################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
install.packages("pdftools")
install.packages("dplyr")
install.packages("readxl")
install.packages("rvest")
install.packages("string")

library(pdftools)
library(limma)
library(dplyr)
library(readxl)
library(rvest)
library(stringr)
library(htmltools)
library(htmlwidgets)
library(tidyverse)
library(writexl)
library(formattable)

#################################################
############### PDF import function #############
#################################################
read.pdf.EE <- function(src)
{
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  
  QTD_COLUMNS <- 2
  read_text <- function(text) {
    result <- ''
    #Get all index of " " from page.
    lstops <- gregexpr(pattern =" ",text)
    #Puts the index of the most frequents ' ' in a vector.
    stops <- as.integer(names(sort(table(unlist(lstops)),decreasing=TRUE)[1:2]))
    #Slice based in the specified number of columns (this can be improved)
    for(i in seq(1, QTD_COLUMNS, by=1))
    {
      temp_result <- sapply(text, function(x){
        start <- 1
        stop <-stops[i] 
        if(i > 1)            
          start <- stops[i-1] + 1
        if(i == QTD_COLUMNS)#last column, read until end.
          stop <- nchar(x)+1
        substr(x, start=start, stop=stop)
      }, USE.NAMES=FALSE)
      temp_result <- trim(temp_result)
      result <- append(result, temp_result)
    }
    result
  }
  
  txt <- pdf_text(src)
  result <- ''
  for (i in 1:length(txt)) { 
    page <- txt[i]
    t1 <- unlist(strsplit(page, "\n"))      
    maxSize <- max(nchar(t1))
    
    t1 <- paste0(t1,strrep(" ", maxSize-nchar(t1)))
    result = append(result,read_text(t1))
  }
  result <- result[!result %in% ""]
  return(result)
}

read.pdf.EE(src)

#################################################
###### RoB mining tool preclinical studies ######
#################################################

### Regex ###
ALSmodel_regex <- c("G93A", "G86R", "G85R", "G37R", "SOD186R",
                    "Tg-SOD1", "SOD1 Tg", "SOD1/rag2", "\\brag2",
                    "Q331K", "\\bFUS\\b",
                    "A315T", "rNLS8", "\\bNLS8\\b", 
                    "wobbler", "\\bwr", "\\bMeHg\\b", "\\bcycad")

MSmodel_regex <- c("cuprizone", "\\bEAE\\b", "experimental autoimmune encephalomyelitis",
                   "experimental allergic encephalomyelitis", "lysolecithin",
                   "delayed type hypersensitivity", "japanese macaque encephalomyelitis",
                   "JME", "lipopolysaccharide", "TMEV", "cytokine injection", "hyponatremia", "NMO",
                   "optic neuritis", "ethidium bromide", "EtBr")

species_list <- c("\\brat[^ihe]", "\\bmarmoset", "\\bguinea pig", "\\bmacaque", "\\brhesus monkey",
                  "\\bswine", "\\bdog", "\\bcat(s|\\b)[^.]", "\\brabbit", "\\bprimate",
                  "\\bhamster", "\\bsheep", "\\bmouse", "\\bmice", "\\bswine", "minipig")

female_regex <- "\\bfemale(s|)\\b"
only_female_regex <- "only.{1,10}\\bfemale(s|)\\b"
male_regex <- "\\bmale(s|)\\b"
only_male_regex <- "only.{1,10}\\bmale(s|)\\b"
sex_regex <- c("either sex", "both genders")

histology_list <- c("histology",  "histological", "histopathology", "histopathologic",
                    "staining", "electron microscopy","Confocal Microscopy",
                    "immunohistochemistry", "immunehistochemistry", 
                    "immunostaining(s|)", "fluorescent staining(s|)", 
                    "immunofluorescence", "immunohistochemical",
                    "Nissl", "Hematoxylin-eosin", "Hematoxylin", "Eosin", "H&E")
behaviour_list <- c("DigiGait", "Treadmill", "kinematics", "grip strength", 
                    "rota(-|)rod", "motor performance", "grip test", "motor coordination",
                    "motor dysfunction", "muscular strength", "gripmeter", "hang wire", 
                    "behavio(u|)r score", "open field", "runtime", "pole test",
                    "catwalk", "gait analysis", "grid hanging", "pole climbing",
                    "climbing pole", "hanging endurance", "neurological score", 
                    "hanging wire", "swim test", "Morris water maze", "Morris maze",
                    "SHIRPA", "reflex test", "walking speed", "stair( |)case", "clinical score", "clinical scores",
                    "whisker behavio(u|)r", "clinical assessment", "clinical signs",
                    "clinical symptoms", "clinical grading", "clinical scoring", "disability score",
                    "examined clinically", "clinical evaluation", "neurological status",
                    "neurological deficit", "visual acuity", "ambulation", 
                    "paralysis", "flaccid", "weakness", "posture", "moribund", "ataxia", 
                    "paraplegia", "paresis", "disability")
imaging_list <- c("MR imaging", "magnetic resonance imaging", "\\bMRI\\b",
                  "\\bNODDI\\b", "diffusion tensor imaging", "\\bDTI\\b",
                  "diffusion( |-)weighted imaging", "\\bDWI\\b", "\\bMTR\b",
                  "diffusion kurtosis imaging")


###Create functions for quantification
Specifity_fun <- function(TP, TN, FP, FN) {
  output <- TN/(TN+FP)
  return(output)
}
Sensitivity_fun <- function(TP, TN, FP, FN) {
  output <- TP/(TP+FN)
  return(output)
}
Precision_fun <- function(TP, TN, FP, FN) {
  output <- TP/(TP+FP)
  return(output)
}
F1_fun <- function(TP, TN, FP, FN) {
  output <- 2*TP/(2*TP+FP+FN)
  return(output)
}
Accuracy_fun <- function(TP, TN, FP, FN) {
  output <- (TP+TN)/(TP+TN+FP+FN)
  return(output)
}

#Rob function
src = pdf.rob.vector[pdf.i]
pdf.extraction <- function(src)
{
  temp.x <- read.pdf.EE(src = src)
  
  #this creates a vector with 1 line as a element
  
  ############################################
  ########## Define paper sections ###########
  ############################################
  
  #Abstract
  
  if(length(grep("Abstract|a b s t r a c t", temp.x, ignore.case = T)) == 0)
  {
    Abstract.start <- 1
  }else{
    Abstract.start <- min(grep("Abstract|a b s t r a c t", temp.x, ignore.case = T))
  }
  
  #Introduction
  
  if(length(grep("Introduction", temp.x, ignore.case = T)) == 0)
  {
    Introduction.start <- 1
  }else{
    Introduction.start <- min(grep("Introduction", temp.x, ignore.case = T))
  }
  
  #Methods
  
  #if(length(grep("materials(\\s+)and(\\s+)methods", temp.x, ignore.case = T)) == 0)
  #{
  #  if(length(grep("experimental(\\s+)procedures", temp.x, ignore.case = T)[grep("experimental(\\s+)procedures", temp.x, ignore.case = T) > 33]) == 0)
  #  {
  #    if(length(grep("methods", temp.x, ignore.case = T)[grep("methods", temp.x, ignore.case = T) > 33]) == 0)
  #    {Method.start <- 33}
  #    else{
  #      Method.start <- min(grep("methods", temp.x, ignore.case = T)[grep("methods", temp.x, ignore.case = T) > 33])
  #    }
  #  }else{
  #    Method.start <- min(grep("experimental(\\s+)procedures", temp.x, ignore.case = T)[grep("experimental(\\s+)procedures", temp.x, ignore.case = T) > 33])
  #  }    
  #}else{
  #  Method.start <- min(grep("materials(\\s+)and(\\s+)methods", temp.x, ignore.case = T))
  #}
  
  ##Removing "experimental procedures" as methods defining feature because
  #"experimental procedures" does commonly not correspond to methods start
  #but is used in a different context
  if(length(grep("materials(\\s+)and(\\s+)methods", temp.x, ignore.case = T)) == 0)
  {
    if(length(grep("methods", temp.x, ignore.case = T)[grep("methods", temp.x, ignore.case = T) > 33]) == 0)
    {
      Method.start <- 33
    }else{
      Method.start <- min(grep("methods", temp.x, ignore.case = T)[grep("methods", temp.x, ignore.case = T) > 33])
    }
  }else{
    Method.start <- min(grep("materials(\\s+)and(\\s+)methods", temp.x, ignore.case = T))
  }
  
  #Results --> most difficult to define results start since "result" appears commonly in scientific texts
  
  if(length(grep("Results", temp.x, ignore.case = T)[grep("Results", temp.x, ignore.case = T) > 33]) == 0)
  {
    Results.start <- 33
  }
  if(length(grep("Results", temp.x, ignore.case = T)[grep("Results", temp.x, ignore.case = T) > 33]) > 0)
  {
    Results.start <- min(grep("Results", temp.x, ignore.case = T)[grep("Results", temp.x, ignore.case = T) > 33])
  }
  
  #Discussion
  
  if(length(grep("Discussion", temp.x, ignore.case = T)[grep("Discussion", temp.x, ignore.case = T) > 33]) == 0)
  {
    if(length(grep("Results", temp.x, ignore.case = T)[grep("Results", temp.x, ignore.case = T) > 33]) == 0)
    {
      Discussion.start <- 33
    }else{
      Discussion.start <-  min(grep("Results", temp.x, ignore.case = T)[grep("Results", temp.x, ignore.case = T) > 33])
    }
  }else{
    Discussion.start <- min(grep("Discussion", temp.x, ignore.case = T)) 
  }
  
  #References
  
  if(length(grep("References", temp.x, ignore.case = T)[grep("References", temp.x, ignore.case = T) > 33]) == 0)
  {
    References.start <- length(temp.x)
  }else{
    References.start <- min(grep("References", temp.x, ignore.case = T)[grep("References", temp.x, ignore.case = T) > 33])
  }
  
  
  ############################################
  ############# Paper mapping ################
  ############################################ 
  
  #Abstract
  
  if(Abstract.start == 1)
  {
    if(length(grep("Introduction", temp.x, ignore.case = T)) == 0){temp.x.abstract <- temp.x[1:80]}
    else{temp.x.abstract <- temp.x[1:Introduction.start]}}
  if(Abstract.start > 1)
  {
    if(length(grep("Introduction", temp.x, ignore.case = T)) == 0){temp.x.abstract <- temp.x[Abstract.start:80]}
    else{temp.x.abstract <- temp.x[Abstract.start:Introduction.start]}}
  
  temp.x.abstract <- paste(unlist(temp.x.abstract), sep = " ")
  temp.x.abstract <- paste(temp.x.abstract, collapse = " ")
  
  #Introduction
  
  if(Introduction.start < Method.start){intro.range <- Introduction.start:Method.start
  }else{intro.range <- Introduction.start+80}#which value covers 90% of all Intros? 80 is a raw estimate
  
  temp.x.intro <- temp.x[intro.range]
  temp.x.intro <- paste(unlist(temp.x.intro), sep = " ")
  temp.x.intro <- paste(temp.x.intro, collapse = " ")
  
  #Experimental (in most cases, you want to search methods and results for our terms, experimental = methods + results)
  
  if(Discussion.start < Method.start)
  {
    experimental.range <- Method.start:References.start
  }else{
    experimental.range <- Method.start:Discussion.start
  }
  
  temp.x.method <- temp.x[experimental.range]
  temp.x.method <- paste(unlist(temp.x.method), sep = " ")
  temp.x.method <- paste(temp.x.method, collapse = " ")
  
  #Whole paper w/o References
  paper.range <- Abstract.start:References.start
  paper.range.full <- 1:length(temp.x)
  temp.x.paper <- temp.x[paper.range]
  temp.x.paper <- paste(unlist(temp.x.paper), sep = " ")
  temp.x.paper <- paste(temp.x.paper, collapse = " ")
  temp.x.paper.full <- temp.x[paper.range.full]
  temp.x.paper.full <- paste(unlist(temp.x.paper.full), sep = " ")
  temp.x.paper.full <- paste(temp.x.paper.full, collapse = " ")
  
  ############################################
  ################ Metadata ##################
  ############################################
  
  #file title
  
  file.title.temp <- str_extract(pdf.rob.vector[pdf.i], regex("/.+", ignore_case = T))
  file.title <- str_extract(file.title.temp, regex("[^/].+"))
  
  #author
  if(str_detect(pdf.rob.vector[pdf.i], regex("/[a-z| ]+-", ignore_case = T)) == F){
    author <- NA
  }else{
    author.temp <- str_extract(pdf.rob.vector[pdf.i], regex("/[a-z| ]+-", ignore_case = T))
    author <- str_extract(author.temp, regex("[a-z]+", ignore_case = T))
  }
  
  #year
  if(str_detect(pdf.rob.vector[pdf.i], regex("[^\\d| ]\\d\\d\\d\\d[^\\d| ]", ignore_case = T)) == T){
    year.temp <- str_extract(pdf.rob.vector[pdf.i], regex("[^\\d| ]\\d\\d\\d\\d[^\\d| ]", ignore_case = T))
    year.temp <- str_extract(year.temp, regex("\\d\\d\\d\\d", ignore_case = T))
    year <- ifelse(year.temp < 2019 & year.temp > 1900, year.temp, NA)
  }else{
    year.temp <- t(strsplit2(temp.x.abstract, " "))
    year.temp <- gsub(";", "", year.temp)
    year.temp <- gsub(",", "", year.temp)
    
    year.temp <- as.numeric(year.temp)
    year.temp <- year.temp[!is.na(year.temp)]
    year.temp <- year.temp[nchar(year.temp) == 4]
    year.temp <- year.temp[ year.temp < 2019]
    year.temp <- year.temp[ year.temp > 1900]
    year.temp <- max(year.temp, na.rm = T)
    year <- ifelse(year.temp < 2019 & year.temp > 1900, year.temp, NA)
  }
  
  #paper title
  if(str_detect(pdf.rob.vector[pdf.i], regex("-[a-z| ]+\\.", ignore_case = T)) == F){
    paper.title <- NA
  }else{
    paper.title.temp <- str_extract(pdf.rob.vector[pdf.i], regex("-[a-z| ]+\\.", ignore_case = T))
    paper.title <- str_extract(paper.title.temp, regex("[a-z| ]+", ignore_case = T))
  }
  
  #doi
  if(str_detect(temp.x.method, regex("10.\\d{4,9}/[-._;()/:a-z0-9A-Z]+", ignore_case = T)) == T){
    doi <- str_extract_all(temp.x.method, regex("10.\\d{4,9}/[-._;()/:a-z0-9A-Z]+", ignore_case = T))
  }else{
    doi <- NA
  }
  
  #E-mail
  if(str_detect(toString(temp.x.paper), regex("([_a-z0-9-]+(\\.[_a-z0-9-]+)*@[a-z0-9-]+(\\.[a-z0-9-]+)*(\\.[a-z]{2,4}))", ignore_case = T)) == F){
    email <- NA
  }else{
    email.temp <- str_extract_all(temp.x.paper, regex("([_a-z0-9-]+(\\.[_a-z0-9-]+)*@[a-z0-9-]+(\\.[a-z0-9-]+)*(\\.[a-z]{2,4}))", ignore_case = T))
    email.temp <- Filter(length, email.temp)
    email <- paste(email.temp, collapse = "; ")}
  
  ############################################
  ########### Experimental data ##############
  ############################################
  
  #Experimental model
  model.count.temp <- vector()
  model.temp <- vector()
  for(numCount in seq_along(ALSmodel_regex)) 
  {
    model.count.temp[[numCount]] <- sum(str_count(temp.x.method, regex(ALSmodel_regex[[numCount]], ignore_case = F)))#ignore case!
    model.temp[[numCount]] <- str_extract(paste(temp.x.method, collapse = " "), regex(ALSmodel_regex[[numCount]], ignore_case = F))#ignore case!
  }
  model.temp.df <- as.data.frame(t(rbind(model.temp, model.count.temp)), stringsAsFactors = F)
  model.temp.df$model.count.temp <- as.numeric(model.count.temp)
  model.temp.df$model.temp <- as.character(model.temp)
  model.temp.df <- model.temp.df[order(desc(model.temp.df$model.count.temp)),]
  
  model1 <- paste(model.temp.df$model.temp[1])
  model_count1 <- paste(model.temp.df$model.count.temp[1])
  model2 <- paste(model.temp.df$model.temp[2])
  model_count2 <- paste(model.temp.df$model.count.temp[2])
  
  if(model1 == "NA"){model1 <- "not reported"}
  if(model2 == "NA"){model2 <- "not reported"}
  
  #Species
  species.count.temp <- vector()
  species.temp <- vector()
  for(numCount in seq_along(species_list)) 
  {
    species.count.temp[[numCount]] <- sum(str_count(temp.x.method, regex(species_list[[numCount]])))
    species.temp[[numCount]] <- str_extract(temp.x.method, regex(species_list[[numCount]]))
  }
  species.temp.df <- as.data.frame(t(rbind(species.temp, species.count.temp)), stringsAsFactors = F)
  species.temp.df$species.count.temp <- as.numeric(species.count.temp)
  species.temp.df$species.temp <- as.character(species.temp)
  species.temp.df <- species.temp.df[order(desc(species.temp.df$species.count.temp)),]
  species.temp.df <- species.temp.df[!(species.temp.df$species.count.temp == 0),]
  
  species1 <- paste(species.temp.df$species.temp[1])
  species_count1 <- paste(species.temp.df$species.count.temp[1])
  species2 <- paste(species.temp.df$species.temp[2])
  species_count2 <- paste(species.temp.df$species.count.temp[2])
  
  if(species1 == "NA"){species1 <- "not reported"}
  if(species2 == "NA"){species2 <- "not reported"}
  
  #Sex
  has_female <- str_subset(temp.x.method, regex(female_regex, ignore_case = T))
  has_female <- str_extract(temp.x.method, regex(female_regex, ignore_case = T))
  has_female <- paste(unlist(has_female), collapse = ";")
  ifelse(str_detect(has_female, regex(female_regex, ignore_case = T)) == T,
         female <- "female", female <- 0)
  
  has_male <- str_subset(temp.x.method, regex(male_regex, ignore_case = T))
  has_male <- str_extract(temp.x.method, regex(male_regex, ignore_case = T))
  has_male <- paste(unlist(has_male), collapse = ";")
  ifelse(str_detect(has_male, regex(male_regex, ignore_case = T)) == T,
         male <- "male", male <- 0)
  
  if(sum(str_detect(female, regex(female_regex, ignore_case = T))) > 0 &
     sum(str_detect(male, regex(male_regex, ignore_case = T))) == 0){
    sex.test <- "female"
  }
  if(sum(str_detect(female, regex(female_regex, ignore_case = T))) == 0 &
     sum(str_detect(male, regex(male_regex, ignore_case = T))) > 0){
    sex.test <- "male"
  }
  if(sum(str_detect(female, regex(female_regex, ignore_case = T))) > 0 &
     sum(str_detect(male, regex(male_regex, ignore_case = T))) > 0){
    sex.test <- "female and male"
  }
  if(sum(str_detect(female, regex(female_regex, ignore_case = T))) == 0 &
     sum(str_detect(male, regex(male_regex, ignore_case = T))) == 0){
    sex.test <- "not reported"
  }
  if(sum(str_detect(temp.x.method, regex(sex_regex, ignore_case = T))) > 0){
    sex.test <- "female and male"
  }
  if(sum(str_detect(temp.x.method, regex(only_female_regex, ignore_case = T))) > 0){
    sex.test <- "female"
  }
  if(sum(str_detect(temp.x.method, regex(only_male_regex, ignore_case = T))) > 0){
    sex.test <- "male"
  }
  
  #Outcome
  has_histology <- str_extract(temp.x.method, regex(histology_list, ignore_case = T))
  has_histology <- paste(unlist(has_histology), collapse = ";")
  ifelse(str_detect(has_histology, regex(histology_list, ignore_case = T)) == F, histology_outcome <- "not reported", histology_outcome <- "histology")
  
  has_behaviour <- str_extract(temp.x.method, regex(behaviour_list, ignore_case = T))
  has_behaviour <- paste(unlist(has_behaviour), collapse = ";")
  ifelse(str_detect(has_behaviour, regex(behaviour_list, ignore_case = T)) == F, behaviour_outcome <- "not reported", behaviour_outcome <- "behaviour")
  
  has_imaging <- str_extract(temp.x.method, regex(imaging_list, ignore_case = T))
  has_imaging <- paste(unlist(has_imaging), collapse = ";")
  ifelse(str_detect(has_imaging, regex(imaging_list, ignore_case = T)) == F, imaging_outcome <- "not reported", imaging_outcome <- "imaging")
  
  #outcome <- paste(c(histology,behaviour,imaging), collapse=", ")
  
  #if(sum(str_detect(outcome, regex("histology, behaviour, imaging", ignore_case = T))) > 0){
  #  outcome.test <- "histology, behaviour, imaging"
  #}
  #if(sum(str_detect(outcome, regex("histology, behaviour, 0", ignore_case = T))) > 0){
  #  outcome.test <- "histology, behaviour"
  #}
  #if(sum(str_detect(outcome, regex("histology, 0, imaging", ignore_case = T))) > 0){
  #  outcome.test <- "histology, imaging"
  #}
  #if(sum(str_detect(outcome, regex("0, behaviour, imaging", ignore_case = T))) > 0){
  #  outcome.test <- "behaviour, imaging"
  #}
  #if(sum(str_detect(outcome, regex("histology, 0, 0", ignore_case = T))) > 0){
  #  outcome.test <- "histology"
  #}
  #if(sum(str_detect(outcome, regex("0, behaviour, 0", ignore_case = T))) > 0){
  #  outcome.test <- "behaviour"
  #}
  #if(sum(str_detect(outcome, regex("0, 0, imaging", ignore_case = T))) > 0){
  #  outcome.test <- "imaging"
  #}
  #if(sum(str_detect(outcome, regex("0, 0, 0", ignore_case = T))) > 0){
  #  outcome.test <- NA
  #}
  
  ############################################
  ############## Risk of Bias ################
  ############################################
  
  #Welfare
  if(str_detect(temp.x.method, regex("approved|approval|carried out in accordance|conducted in accordance|in accordance|animal experiments were performed|conformed to guidelines|complied with relevant|in compliance with the|animal handling conformed|met the state regulations|in agreement with|in strict adherence|in adherence", ignore_case = T)) == T){
    welfare <- "yes"
  }else{
    welfare <- "not reported"
  }
  
  if(sum(str_count(temp.x.method, regex("approved|approval|carried out in accordance|conducted in accordance|performed in strict concordance|in accordance|according to guidelines|conformed to the guidelines", ignore_case = T))) > 1){
    welfare_signal <- "Manual check recommended"
  }else{
    welfare_signal <- "No manual check required"
  }
  
  #Blinding
  if(str_detect(temp.x.method, regex("blinded|were blind|blind", ignore_case = T)) == T){
    blinding <- "yes"
  }else{
    blinding <- "not reported"
  }
  ##CAVE; NOT blindness
  if(sum(str_count(temp.x.method, regex("blinded|were blind|blind", ignore_case = T))) > 1){
    blinding_signal <- "Manual check recommended"
  }else{
    blinding_signal <- "No manual check required"
  }
  
  #Randomization
  if(str_detect(temp.x.method, regex("random", ignore_case = T)) == T){
    randomization <- "yes"
  }else{
    randomization <- "not reported"
  }
   ##CAVE: NOT RANDOM EFFECTS OR RANDOM PRIMER(S), RANDOM HEXAMER(S); RANDOM Coefficient(s), random variable(s), random variation, random network, randomized-controlled, randomized controlled, randomized\, double
  
  if(sum(str_count(temp.x.method, regex("random", ignore_case = T))) > 1){
    randomization_signal <- "Manual check recommended"
  }else{
    randomization_signal <- "No manual check required"
  }
  
  #sample size calculation
  if(str_detect(temp.x.paper.full, regex("achieve.{1,10}power|sample size calculation|power analysis|study power|sample size|power calculation", ignore_case = T)) == T){
    powertest <- "yes"
  }else{
    powertest <- "not reported"
  }
  ##NOT powerful, power level, wave power
  
  if(sum(str_count(temp.x.paper.full, regex("achieve.{1,10}power|sample size calculation|power analysis|study power|sample size", ignore_case = T))) > 1){
    powertest_signal <- "Manual check recommended"
  }else{
    powertest_signal <- "No manual check required"
  }
  
  #ARRIVE guidelines
  if(str_detect(temp.x.paper, regex("ARRIVE", ignore_case = F)) == T){#ignore case!
    ARRIVE <- "yes"
  }else{
    ARRIVE <- "not reported"
  }
  ##CAVE: AFTER 2010
  if(sum(str_count(temp.x.paper, regex("ARRIVE", ignore_case = F))) > 1){
    ARRIVE_signal <- "Manual check recommended"
  }else{
    ARRIVE_signal <- "No manual check required"
  }
  
  #conflict of interest
  if(str_detect(temp.x.paper.full, regex("conflict[s| ] ?of interest|disclosures|disclosure|competing interests|nothing to disclose|received grant[s| ]|disclosure statement|financial disclosure|financial interests", ignore_case = T)) == T){
    disclosure <- "yes"
  }else{
    disclosure <- "not reported"
  }
  
  if(sum(str_count(temp.x.paper.full, regex("conflict[s| ] ?of interest|disclosures|competing interests|nothing to disclose|received grant[s| ]|disclosure statement", ignore_case = T))) > 1){
    disclosure_signal <- "Manual check recommended"
  }else{
    disclosure_signal <- "No manual check required"
  }
  
  #Data availability
  if(str_detect(temp.x.paper.full, regex("data availability|availability of data", ignore_case = T)) == T){
    availability <- "yes"
  }else{
    availability <- "not reported"
  }
  
  ############################################
  ################# OUTPUT ###################
  ############################################
  
  return.x <- c(file.title, author, year, paper.title, 
                model1, model_count1, model2, model_count2,
                species1, species_count1, species2, species_count2,
                sex.test, histology_outcome, behaviour_outcome, imaging_outcome,
                randomization, randomization_signal,
                blinding, blinding_signal,
                welfare, welfare_signal,
                disclosure, disclosure_signal,
                powertest, powertest_signal,
                ARRIVE, ARRIVE_signal,
                availability,
                doi[[1]][1], email,
                Abstract.start, Introduction.start, Method.start, Results.start, Discussion.start,
                paste(c(min(experimental.range), max(experimental.range)), collapse=" - "),
                paste(c(min(paper.range.full), max(paper.range.full)), collapse=" - "))
  names(return.x) <- c("File title", "First author", "Year","Paper title", 
                       "Model 1", "Model count 1", "Model 2", "Model count 2",
                       "Species 1", "Species count 1", "Species 2", "Species count 2",
                       "Sex", "Outcome histology", "Outcome behaviour", "Outcome imaging",
                       "Randomization", "Randomization_QC",
                       "Blinding", "Blinding_QC",
                       "Welfare", "Welfare_QC",
                       "Conflict", "Conflict_QC",
                       "Samplesize", "Sample_QC",
                       "ARRIVE", "ARRIVE_QC",
                       "Data_availability_statement",
                       "doi", "email", 
                       "Abstract start", "Introduction start", "Methods start", "Results start", "Discussion start",
                       "Methods range", "Paper range")
  return(return.x)
}


#################################################
####### Apply function to reference sets ########
#################################################
#test set 1: define folder with PDFs
pdf.rob.vector <- list.files("C:/Users/benja/Dropbox/Research/Experiments/Karolinska/ALS_SystematicReview/Animal_study/Validation_sample2/", full.names = T)
pdf.rob.vector <- pdf.rob.vector[grep("pdf", pdf.rob.vector, ignore.case = T)]
pdf.rob.vector <- pdf.rob.vector[!grepl("pdftools", pdf.rob.vector)]
pdf.rob.vector

#Extract: apply function to PDF test set
mined.rob <- vector()

for(pdf.i in 1:length(pdf.rob.vector)) 
{
  print(pdf.i)
  print(pdf.rob.vector[pdf.i])
  mined.rob <-  rbind(mined.rob, pdf.extraction(src = pdf.rob.vector[pdf.i]))
}


pdf.extraction(src = pdf.rob.vector[1])

ret <- lapply(pdf.rob.vector, pdf.extraction)
df <- as.data.frame(do.call("rbind", ret))
str(df)

mined.rob
write_xlsx(as.data.frame(mined.rob), "C:/Users/benja/Dropbox/Research/Experiments/Karolinska/ALS_SystematicReview/Animal_study/rob-pilot.xlsx")

read.pdf.EE(pdf.rob.vector[14])
options(max.print = 2000)

#################################################
############# Accuracy assessment ###############
#################################################

#evaluate accuracy of mining
#load extracted data set
rob.extracted <- read_excel("C:/Users/benja/Dropbox/Research/Experiments/Karolinska/ALS_SystematicReview/Animal_study/Rob_test1.xlsx", sheet = 5)
drops <- c("Title")
rob.extracted <- rob.extracted[ , !(names(rob.extracted) %in% drops)]

#compared with mined data set
keeps <- c("File title", "Randomization", "Blinding", "Welfare", "Conflict",
           "Samplesize", "ARRIVE", "Data_availability_statement",
           "Species 1", "Sex", "Model 1", 
           "Outcome histology", "Outcome behaviour", "Outcome imaging")
rob.mined <- mined.rob[ , keeps, drop = FALSE]
colnames(rob.mined) <- c("Reference", "Randomization_m", "Blinding_m", "Welfare_m",
                         "Conflict_m", "Samplesize_m", "ARRIVE_m", "Data_availability_m",
                         "Species_m", "Sex_m", "Model_m",
                         "Outcome_histology_m", "Outcome_behaviour_m", "Outcome_imaging_m")
rob.mined <- gsub("Users/benja/Dropbox/Research/Experiments/Karolinska/ALS_SystematicReview/Animal_study/Training_sample/", "", rob.mined) 

rob.mined == rob.extracted

#Primary test for rob function
rob.mined <- subset(rob.mined, select = c("Randomization_m", "Blinding_m", "Welfare_m", "Conflict_m", "Samplesize_m", "ARRIVE_m",
                                          "Data_availability_m", "Species_m", "Sex_m", "Model_m",
                                          "Outcome_histology_m", "Outcome_behaviour_m", "Outcome_imaging_m"))
rob.extracted <- subset(rob.extracted, select = c("Randomization_e", "Blinding_e", "Welfare_e", "Conflict_e", "Samplesize_e", "ARRIVE_e",
                                                  "Data_availability_e", "Species_e", "Sex_e", "Model_e",
                                                  "Outcome_histology_e", "Outcome_behaviour_e", "Outcome_imaging_e"))
rob.mined.df <- as.data.frame(rob.mined)
rob.extracted.df <- as.data.frame(rob.extracted)
rob.mined.df == rob.extracted.df
rob.comparison <- table(rob.mined.df == rob.extracted.df)
rob.comparison[names(rob.comparison)==F]

quantification_df <- cbind(rob.extracted, rob.mined)
quantification_df$Species_m <- gsub("mice", "mouse", quantification_df$Species_m)
quantification_df$Species_m <- gsub("rats", "rat", quantification_df$Species_m)

quantification_df$Randomization <- ifelse(quantification_df$Randomization_e == "not reported" &
                                     quantification_df$Randomization_m == "not reported", "TN",
                                   ifelse(quantification_df$Randomization_e == "yes" &
                                            quantification_df$Randomization_m == "yes", "TP",
                                          ifelse(quantification_df$Randomization_e == "not reported" &
                                                   quantification_df$Randomization_m == "yes", "FP", "FN")))
quantification_df$Blinding <- ifelse(quantification_df$Blinding_e == "not reported" &
                                     quantification_df$Blinding_m == "not reported", "TN",
                                   ifelse(quantification_df$Blinding_e == "yes" &
                                            quantification_df$Blinding_m == "yes", "TP",
                                          ifelse(quantification_df$Blinding_e == "not reported" &
                                                   quantification_df$Blinding_m == "yes", "FP", "FN")))
quantification_df$Welfare <- ifelse(quantification_df$Welfare_e == "not reported" &
                                     quantification_df$Welfare_m == "not reported", "TN",
                                   ifelse(quantification_df$Welfare_e == "yes" &
                                            quantification_df$Welfare_m == "yes", "TP",
                                          ifelse(quantification_df$Welfare_e == "not reported" &
                                                   quantification_df$Welfare_m == "yes", "FP", "FN")))
quantification_df$Conflict <- ifelse(quantification_df$Conflict_e == "not reported" &
                                     quantification_df$Conflict_m == "not reported", "TN",
                                   ifelse(quantification_df$Conflict_e == "yes" &
                                            quantification_df$Conflict_m == "yes", "TP",
                                          ifelse(quantification_df$Conflict_e == "not reported" &
                                                   quantification_df$Conflict_m == "yes", "FP", "FN")))
quantification_df$Samplesize <- ifelse(quantification_df$Samplesize_e == "not reported" &
                                     quantification_df$Samplesize_m == "not reported", "TN",
                                   ifelse(quantification_df$Samplesize_e == "yes" &
                                            quantification_df$Samplesize_m == "yes", "TP",
                                          ifelse(quantification_df$Samplesize_e == "not reported" &
                                                   quantification_df$Samplesize_m == "yes", "FP", "FN")))
quantification_df$ARRIVE <- ifelse(quantification_df$ARRIVE_e == "not reported" &
                                     quantification_df$ARRIVE_m == "not reported", "TN",
                                   ifelse(quantification_df$ARRIVE_e == "yes" &
                                            quantification_df$ARRIVE_m == "yes", "TP",
                                          ifelse(quantification_df$ARRIVE_e == "not reported" &
                                                   quantification_df$ARRIVE_m == "yes", "FP", "FN")))
quantification_df$Data_availability <- ifelse(quantification_df$Data_availability_e == "not reported" &
                                     quantification_df$Data_availability_m == "not reported", "TN",
                                   ifelse(quantification_df$Data_availability_e == "yes" &
                                            quantification_df$Data_availability_m == "yes", "TP",
                                          ifelse(quantification_df$Data_availability_e == "not reported" &
                                                   quantification_df$Data_availability_m == "yes", "FP", "FN")))
quantification_df$Species <- ifelse(quantification_df$Species_e == "not reported" & quantification_df$Species_m == "not reported", "TN",
                                ifelse(quantification_df$Species_e == quantification_df$Species_m, "TP",
                                       ifelse(quantification_df$Species_e == "not reported" &
                                                quantification_df$Species_m != "not reported", "FP", "FN")))
quantification_df$Sex <- ifelse(quantification_df$Sex_e == "not reported" & quantification_df$Sex_m == "not reported", "TN",
                                              ifelse(quantification_df$Sex_e == quantification_df$Sex_m, "TP",
                                                ifelse(quantification_df$Sex_e == "not reported" &
                                                        quantification_df$Sex_m != "not reported", "FP", "FN")))
quantification_df$Model <- ifelse(quantification_df$Model_e == "not reported" & quantification_df$Model_m == "not reported", "TN",
                                              ifelse(quantification_df$Model_e == quantification_df$Model_m, "TP",
                                                     ifelse(quantification_df$Model_e == "not reported" &
                                                              quantification_df$Model_m != "not reported", "FP", "FN")))
quantification_df$Outcome_histology <- ifelse(quantification_df$Outcome_histology_e == "not reported" & quantification_df$Outcome_histology_m == "not reported", "TN",
                                ifelse(quantification_df$Outcome_histology_e == quantification_df$Outcome_histology_m, "TP",
                                       ifelse(quantification_df$Outcome_histology_e == "not reported" &
                                                quantification_df$Outcome_histology_m != "not reported", "FP", "FN")))
quantification_df$Outcome_behaviour <- ifelse(quantification_df$Outcome_behaviour_e == "not reported" & quantification_df$Outcome_behaviour_m == "not reported", "TN",
                                              ifelse(quantification_df$Outcome_behaviour_e == quantification_df$Outcome_behaviour_m, "TP",
                                                     ifelse(quantification_df$Outcome_behaviour_e == "not reported" &
                                                              quantification_df$Outcome_behaviour_m != "not reported", "FP", "FN")))
quantification_df$Outcome_imaging <- ifelse(quantification_df$Outcome_imaging_e == "not reported" & quantification_df$Outcome_imaging_m == "not reported", "TN",
                                              ifelse(quantification_df$Outcome_imaging_e == quantification_df$Outcome_imaging_m, "TP",
                                                     ifelse(quantification_df$Outcome_imaging_e == "not reported" &
                                                              quantification_df$Outcome_imaging_m != "not reported", "FP", "FN")))

quantification_df <- quantification_df[, c(1, 14, 27, 2, 15, 28,
                                           3, 16, 29, 4, 17, 30,
                                           5, 18, 31, 6, 19, 32,
                                           7, 20, 33, 8, 21, 34,
                                           9, 22, 35, 10, 23, 36,
                                           11, 24, 37, 12, 25, 38,
                                           13, 26, 39)]
table(quantification_df$Samplesize)

formattable(quantification_df, list(
  Sex = formatter("span", 
                         style = ~style(display = "block",
                                        font.weight = "bold", 
                                        color = "white",
                                        "border-radius" = "4px",
                                        "padding-right" = "4px",
                                        "background-color" =  
                                          ifelse(Sex == "FP" | Sex == "FN", "red", "blue")))
))


###Calculate precision measures
Specificity_Random <- Specifity_fun(length(which(quantification_df$Randomization == "TP")), 
              length(which(quantification_df$Randomization == "TN")),
              length(which(quantification_df$Randomization == "FP")),
              length(which(quantification_df$Randomization == "FN")))
Sensitivity_Random <- Sensitivity_fun(length(which(quantification_df$Randomization == "TP")), 
              length(which(quantification_df$Randomization == "TN")),
              length(which(quantification_df$Randomization == "FP")),
              length(which(quantification_df$Randomization == "FN")))
Precision_Random <- Precision_fun(length(which(quantification_df$Randomization == "TP")), 
              length(which(quantification_df$Randomization == "TN")),
              length(which(quantification_df$Randomization == "FP")),
              length(which(quantification_df$Randomization == "FN")))
F1_Random <- F1_fun(length(which(quantification_df$Randomization == "TP")), 
              length(which(quantification_df$Randomization == "TN")),
              length(which(quantification_df$Randomization == "FP")),
              length(which(quantification_df$Randomization == "FN")))
Accuracy_Random <- Accuracy_fun(length(which(quantification_df$Randomization == "TP")), 
                    length(which(quantification_df$Randomization == "TN")),
                    length(which(quantification_df$Randomization == "FP")),
                    length(which(quantification_df$Randomization == "FN")))

Specificity_Blinding <- Specifity_fun(length(which(quantification_df$Blinding == "TP")), 
                                    length(which(quantification_df$Blinding == "TN")),
                                    length(which(quantification_df$Blinding == "FP")),
                                    length(which(quantification_df$Blinding == "FN")))
Sensitivity_Blinding <- Sensitivity_fun(length(which(quantification_df$Blinding == "TP")), 
                                      length(which(quantification_df$Blinding == "TN")),
                                      length(which(quantification_df$Blinding == "FP")),
                                      length(which(quantification_df$Blinding == "FN")))
Precision_Blinding <- Precision_fun(length(which(quantification_df$Blinding == "TP")), 
                                      length(which(quantification_df$Blinding == "TN")),
                                      length(which(quantification_df$Blinding == "FP")),
                                      length(which(quantification_df$Blinding == "FN")))
F1_Blinding <- F1_fun(length(which(quantification_df$Blinding == "TP")), 
                    length(which(quantification_df$Blinding == "TN")),
                    length(which(quantification_df$Blinding == "FP")),
                    length(which(quantification_df$Blinding == "FN")))
Accuracy_Blinding <- Accuracy_fun(length(which(quantification_df$Blinding == "TP")), 
                      length(which(quantification_df$Blinding == "TN")),
                      length(which(quantification_df$Blinding == "FP")),
                      length(which(quantification_df$Blinding == "FN")))

Specificity_Welfare <- Specifity_fun(length(which(quantification_df$Welfare == "TP")), 
                                      length(which(quantification_df$Welfare == "TN")),
                                      length(which(quantification_df$Welfare == "FP")),
                                      length(which(quantification_df$Welfare == "FN")))
Sensitivity_Welfare <- Sensitivity_fun(length(which(quantification_df$Welfare == "TP")), 
                                        length(which(quantification_df$Welfare == "TN")),
                                        length(which(quantification_df$Welfare == "FP")),
                                        length(which(quantification_df$Welfare == "FN")))
Precision_Welfare <- Precision_fun(length(which(quantification_df$Welfare == "TP")), 
                                        length(which(quantification_df$Welfare == "TN")),
                                        length(which(quantification_df$Welfare == "FP")),
                                        length(which(quantification_df$Welfare == "FN")))
F1_Welfare <- F1_fun(length(which(quantification_df$Welfare == "TP")), 
                      length(which(quantification_df$Welfare == "TN")),
                      length(which(quantification_df$Welfare == "FP")),
                      length(which(quantification_df$Welfare == "FN")))
Accuracy_Welfare <- Accuracy_fun(length(which(quantification_df$Welfare == "TP")), 
                     length(which(quantification_df$Welfare == "TN")),
                     length(which(quantification_df$Welfare == "FP")),
                     length(which(quantification_df$Welfare == "FN")))

Specificity_Conflict <- Specifity_fun(length(which(quantification_df$Conflict == "TP")), 
                                     length(which(quantification_df$Conflict == "TN")),
                                     length(which(quantification_df$Conflict == "FP")),
                                     length(which(quantification_df$Conflict == "FN")))
Sensitivity_Conflict <- Sensitivity_fun(length(which(quantification_df$Conflict == "TP")), 
                                       length(which(quantification_df$Conflict == "TN")),
                                       length(which(quantification_df$Conflict == "FP")),
                                       length(which(quantification_df$Conflict == "FN")))
Precision_Conflict <- Precision_fun(length(which(quantification_df$Conflict == "TP")), 
                                       length(which(quantification_df$Conflict == "TN")),
                                       length(which(quantification_df$Conflict == "FP")),
                                       length(which(quantification_df$Conflict == "FN")))
F1_Conflict <- F1_fun(length(which(quantification_df$Conflict == "TP")), 
                     length(which(quantification_df$Conflict == "TN")),
                     length(which(quantification_df$Conflict == "FP")),
                     length(which(quantification_df$Conflict == "FN")))
Accuracy_Conflict <- Accuracy_fun(length(which(quantification_df$Conflict == "TP")), 
                      length(which(quantification_df$Conflict == "TN")),
                      length(which(quantification_df$Conflict == "FP")),
                      length(which(quantification_df$Conflict == "FN")))

Specificity_Samplesize <- Specifity_fun(length(which(quantification_df$Samplesize == "TP")), 
                                      length(which(quantification_df$Samplesize == "TN")),
                                      length(which(quantification_df$Samplesize == "FP")),
                                      length(which(quantification_df$Samplesize == "FN")))
Sensitivity_Samplesize <- Sensitivity_fun(length(which(quantification_df$Samplesize == "TP")), 
                                        length(which(quantification_df$Samplesize == "TN")),
                                        length(which(quantification_df$Samplesize == "FP")),
                                        length(which(quantification_df$Samplesize == "FN")))
Precision_Samplesize <- Precision_fun(length(which(quantification_df$Samplesize == "TP")), 
                                        length(which(quantification_df$Samplesize == "TN")),
                                        length(which(quantification_df$Samplesize == "FP")),
                                        length(which(quantification_df$Samplesize == "FN")))
F1_Samplesize <- F1_fun(length(which(quantification_df$Samplesize == "TP")), 
                      length(which(quantification_df$Samplesize == "TN")),
                      length(which(quantification_df$Samplesize == "FP")),
                      length(which(quantification_df$Samplesize == "FN")))
Accuracy_Samplesize <- Accuracy_fun(length(which(quantification_df$Samplesize == "TP")), 
                        length(which(quantification_df$Samplesize == "TN")),
                        length(which(quantification_df$Samplesize == "FP")),
                        length(which(quantification_df$Samplesize == "FN")))

Specificity_ARRIVE <- Specifity_fun(length(which(quantification_df$ARRIVE == "TP")), 
                                        length(which(quantification_df$ARRIVE == "TN")),
                                        length(which(quantification_df$ARRIVE == "FP")),
                                        length(which(quantification_df$ARRIVE == "FN")))
Sensitivity_ARRIVE <- Sensitivity_fun(length(which(quantification_df$ARRIVE == "TP")), 
                                          length(which(quantification_df$ARRIVE == "TN")),
                                          length(which(quantification_df$ARRIVE == "FP")),
                                          length(which(quantification_df$ARRIVE == "FN")))
Precision_ARRIVE <- Precision_fun(length(which(quantification_df$ARRIVE == "TP")), 
                                          length(which(quantification_df$ARRIVE == "TN")),
                                          length(which(quantification_df$ARRIVE == "FP")),
                                          length(which(quantification_df$ARRIVE == "FN")))
F1_ARRIVE <- F1_fun(length(which(quantification_df$ARRIVE == "TP")), 
                        length(which(quantification_df$ARRIVE == "TN")),
                        length(which(quantification_df$ARRIVE == "FP")),
                        length(which(quantification_df$ARRIVE == "FN")))
Accuracy_ARRIVE <- Accuracy_fun(length(which(quantification_df$ARRIVE == "TP")), 
                    length(which(quantification_df$ARRIVE == "TN")),
                    length(which(quantification_df$ARRIVE == "FP")),
                    length(which(quantification_df$ARRIVE == "FN")))

Specificity_Data_availability <- Specifity_fun(length(which(quantification_df$Data_availability == "TP")), 
                                    length(which(quantification_df$Data_availability == "TN")),
                                    length(which(quantification_df$Data_availability == "FP")),
                                    length(which(quantification_df$Data_availability == "FN")))
Sensitivity_Data_availability <- Sensitivity_fun(length(which(quantification_df$Data_availability == "TP")), 
                                      length(which(quantification_df$Data_availability == "TN")),
                                      length(which(quantification_df$Data_availability == "FP")),
                                      length(which(quantification_df$Data_availability == "FN")))
Precision_Data_availability <- Precision_fun(length(which(quantification_df$Data_availability == "TP")), 
                                      length(which(quantification_df$Data_availability == "TN")),
                                      length(which(quantification_df$Data_availability == "FP")),
                                      length(which(quantification_df$Data_availability == "FN")))
F1_Data_availability <- F1_fun(length(which(quantification_df$Data_availability == "TP")), 
                    length(which(quantification_df$Data_availability == "TN")),
                    length(which(quantification_df$Data_availability == "FP")),
                    length(which(quantification_df$Data_availability == "FN")))
Accuracy_Data_availability <- Accuracy_fun(length(which(quantification_df$Data_availability == "TP")), 
                               length(which(quantification_df$Data_availability == "TN")),
                               length(which(quantification_df$Data_availability == "FP")),
                               length(which(quantification_df$Data_availability == "FN")))

Specificity_Species <- Specifity_fun(length(which(quantification_df$Species == "TP")), 
                                               length(which(quantification_df$Species == "TN")),
                                               length(which(quantification_df$Species == "FP")),
                                               length(which(quantification_df$Species == "FN")))
Sensitivity_Species <- Sensitivity_fun(length(which(quantification_df$Species == "TP")), 
                                                 length(which(quantification_df$Species == "TN")),
                                                 length(which(quantification_df$Species == "FP")),
                                                 length(which(quantification_df$Species == "FN")))
Precision_Species <- Precision_fun(length(which(quantification_df$Species == "TP")), 
                                                 length(which(quantification_df$Species == "TN")),
                                                 length(which(quantification_df$Species == "FP")),
                                                 length(which(quantification_df$Species == "FN")))
F1_Species <- F1_fun(length(which(quantification_df$Species == "TP")), 
                               length(which(quantification_df$Species == "TN")),
                               length(which(quantification_df$Species == "FP")),
                               length(which(quantification_df$Species == "FN")))
Accuracy_Species <- Accuracy_fun(length(which(quantification_df$Species == "TP")), 
                     length(which(quantification_df$Species == "TN")),
                     length(which(quantification_df$Species == "FP")),
                     length(which(quantification_df$Species == "FN")))

Specificity_Sex <- Specifity_fun(length(which(quantification_df$Sex == "TP")), 
                                     length(which(quantification_df$Sex == "TN")),
                                     length(which(quantification_df$Sex == "FP")),
                                     length(which(quantification_df$Sex == "FN")))
Sensitivity_Sex <- Sensitivity_fun(length(which(quantification_df$Sex == "TP")), 
                                       length(which(quantification_df$Sex == "TN")),
                                       length(which(quantification_df$Sex == "FP")),
                                       length(which(quantification_df$Sex == "FN")))
Precision_Sex <- Precision_fun(length(which(quantification_df$Sex == "TP")), 
                                       length(which(quantification_df$Sex == "TN")),
                                       length(which(quantification_df$Sex == "FP")),
                                       length(which(quantification_df$Sex == "FN")))
F1_Sex <- F1_fun(length(which(quantification_df$Sex == "TP")), 
                     length(which(quantification_df$Sex == "TN")),
                     length(which(quantification_df$Sex == "FP")),
                     length(which(quantification_df$Sex == "FN")))
Accuracy_Sex <- Accuracy_fun(length(which(quantification_df$Sex == "TP")), 
                 length(which(quantification_df$Sex == "TN")),
                 length(which(quantification_df$Sex == "FP")),
                 length(which(quantification_df$Sex == "FN")))

Specificity_Model <- Specifity_fun(length(which(quantification_df$Model == "TP")), 
                                 length(which(quantification_df$Model == "TN")),
                                 length(which(quantification_df$Model == "FP")),
                                 length(which(quantification_df$Model == "FN")))
Sensitivity_Model <- Sensitivity_fun(length(which(quantification_df$Model == "TP")), 
                                   length(which(quantification_df$Model == "TN")),
                                   length(which(quantification_df$Model == "FP")),
                                   length(which(quantification_df$Model == "FN")))
Precision_Model <- Precision_fun(length(which(quantification_df$Model == "TP")), 
                                   length(which(quantification_df$Model == "TN")),
                                   length(which(quantification_df$Model == "FP")),
                                   length(which(quantification_df$Model == "FN")))
F1_Model <- F1_fun(length(which(quantification_df$Model == "TP")), 
                 length(which(quantification_df$Model == "TN")),
                 length(which(quantification_df$Model == "FP")),
                 length(which(quantification_df$Model == "FN")))
Accuracy_Model <- Accuracy_fun(length(which(quantification_df$Model == "TP")), 
                   length(which(quantification_df$Model == "TN")),
                   length(which(quantification_df$Model == "FP")),
                   length(which(quantification_df$Model == "FN")))

Specificity_Outcome_histology <- Specifity_fun(length(which(quantification_df$Outcome_histology == "TP")), 
                                   length(which(quantification_df$Outcome_histology == "TN")),
                                   length(which(quantification_df$Outcome_histology == "FP")),
                                   length(which(quantification_df$Outcome_histology == "FN")))
Sensitivity_Outcome_histology <- Sensitivity_fun(length(which(quantification_df$Outcome_histology == "TP")), 
                                     length(which(quantification_df$Outcome_histology == "TN")),
                                     length(which(quantification_df$Outcome_histology == "FP")),
                                     length(which(quantification_df$Outcome_histology == "FN")))
Precision_Outcome_histology <- Precision_fun(length(which(quantification_df$Outcome_histology == "TP")), 
                                     length(which(quantification_df$Outcome_histology == "TN")),
                                     length(which(quantification_df$Outcome_histology == "FP")),
                                     length(which(quantification_df$Outcome_histology == "FN")))
F1_Outcome_histology <- F1_fun(length(which(quantification_df$Outcome_histology == "TP")), 
                   length(which(quantification_df$Outcome_histology == "TN")),
                   length(which(quantification_df$Outcome_histology == "FP")),
                   length(which(quantification_df$Outcome_histology == "FN")))
Accuracy_Outcome_histology <- Accuracy_fun(length(which(quantification_df$Outcome_histology == "TP")), 
                               length(which(quantification_df$Outcome_histology == "TN")),
                               length(which(quantification_df$Outcome_histology == "FP")),
                               length(which(quantification_df$Outcome_histology == "FN")))

Specificity_Outcome_behaviour <- Specifity_fun(length(which(quantification_df$Outcome_behaviour == "TP")), 
                                               length(which(quantification_df$Outcome_behaviour == "TN")),
                                               length(which(quantification_df$Outcome_behaviour == "FP")),
                                               length(which(quantification_df$Outcome_behaviour == "FN")))
Sensitivity_Outcome_behaviour <- Sensitivity_fun(length(which(quantification_df$Outcome_behaviour == "TP")), 
                                                 length(which(quantification_df$Outcome_behaviour == "TN")),
                                                 length(which(quantification_df$Outcome_behaviour == "FP")),
                                                 length(which(quantification_df$Outcome_behaviour == "FN")))
Precision_Outcome_behaviour <- Precision_fun(length(which(quantification_df$Outcome_behaviour == "TP")), 
                                             length(which(quantification_df$Outcome_behaviour == "TN")),
                                             length(which(quantification_df$Outcome_behaviour == "FP")),
                                             length(which(quantification_df$Outcome_behaviour == "FN")))
F1_Outcome_behaviour <- F1_fun(length(which(quantification_df$Outcome_behaviour == "TP")), 
                               length(which(quantification_df$Outcome_behaviour == "TN")),
                               length(which(quantification_df$Outcome_behaviour == "FP")),
                               length(which(quantification_df$Outcome_behaviour == "FN")))
Accuracy_Outcome_behaviour <- Accuracy_fun(length(which(quantification_df$Outcome_behaviour == "TP")), 
                               length(which(quantification_df$Outcome_behaviour == "TN")),
                               length(which(quantification_df$Outcome_behaviour == "FP")),
                               length(which(quantification_df$Outcome_behaviour == "FN")))

Specificity_Outcome_imaging <- Specifity_fun(length(which(quantification_df$Outcome_imaging == "TP")), 
                                               length(which(quantification_df$Outcome_imaging == "TN")),
                                               length(which(quantification_df$Outcome_imaging == "FP")),
                                               length(which(quantification_df$Outcome_imaging == "FN")))
Sensitivity_Outcome_imaging <- Sensitivity_fun(length(which(quantification_df$Outcome_imaging == "TP")), 
                                                 length(which(quantification_df$Outcome_imaging == "TN")),
                                                 length(which(quantification_df$Outcome_imaging == "FP")),
                                                 length(which(quantification_df$Outcome_imaging == "FN")))
Precision_Outcome_imaging <- Precision_fun(length(which(quantification_df$Outcome_imaging == "TP")), 
                                             length(which(quantification_df$Outcome_imaging == "TN")),
                                             length(which(quantification_df$Outcome_imaging == "FP")),
                                             length(which(quantification_df$Outcome_imaging == "FN")))
F1_Outcome_imaging <- F1_fun(length(which(quantification_df$Outcome_imaging == "TP")), 
                               length(which(quantification_df$Outcome_imaging == "TN")),
                               length(which(quantification_df$Outcome_imaging == "FP")),
                               length(which(quantification_df$Outcome_imaging == "FN")))
Accuracy_Outcome_imaging <- Accuracy_fun(length(which(quantification_df$Outcome_imaging == "TP")), 
                             length(which(quantification_df$Outcome_imaging == "TN")),
                             length(which(quantification_df$Outcome_imaging == "FP")),
                             length(which(quantification_df$Outcome_imaging == "FN")))

###plot precision measures
Table_quantification <- matrix(c(Specificity_Random, Sensitivity_Random, Precision_Random, F1_Random, Accuracy_Random,
                                Specificity_Blinding, Sensitivity_Blinding, Precision_Blinding, F1_Blinding, Accuracy_Blinding,
                                Specificity_Welfare, Sensitivity_Welfare, Precision_Welfare, F1_Welfare, Accuracy_Welfare,
                                Specificity_Conflict, Sensitivity_Conflict, Precision_Conflict, F1_Conflict, Accuracy_Conflict,
                                Specificity_Samplesize, Sensitivity_Samplesize, Precision_Samplesize, F1_Samplesize, Accuracy_Samplesize,
                                Specificity_ARRIVE, Sensitivity_ARRIVE, Precision_ARRIVE, F1_ARRIVE, Accuracy_ARRIVE,
                                Specificity_Data_availability, Sensitivity_Data_availability, Precision_Data_availability, F1_Data_availability, Accuracy_Data_availability,
                                Specificity_Species, Sensitivity_Species, Precision_Species, F1_Species, Accuracy_Species,
                                Specificity_Sex, Sensitivity_Sex, Precision_Sex, F1_Sex, Accuracy_Sex,
                                Specificity_Model, Sensitivity_Model, Precision_Model, F1_Model, Accuracy_Model,
                                Specificity_Outcome_histology, Sensitivity_Outcome_histology, Precision_Outcome_histology, F1_Outcome_histology, Accuracy_Outcome_histology,
                                Specificity_Outcome_behaviour, Sensitivity_Outcome_behaviour, Precision_Outcome_behaviour, F1_Outcome_behaviour, Accuracy_Outcome_behaviour,
                                Specificity_Outcome_imaging, Sensitivity_Outcome_imaging, Precision_Outcome_imaging, F1_Outcome_imaging, Accuracy_Outcome_imaging),
                               ncol=13)

rownames(Table_quantification) <- c("Specificity", "Sensitivity", "Precision", "F1", "Accuracy")
colnames(Table_quantification) <- c("Randomization", "Blinding", "Welfare",
                                    "Conflict", "Samplesize", "ARRIVE",
                                    "Data availability", "Species", "Sex",
                                    "Model", "Outcome histology", "Outcome behaviour", "Outcome imaging")
Table_quantification <- as.table(Table_quantification)
Table_quantification



#accuracy in percent
as.numeric(rob.comparison[names(rob.comparison)==T])/(nrow(rob.mined)*ncol(rob.mined))*100

test <- read.pdf.EE(pdf.rob.vector[14])
test1 <- "rats"

if (test1 == "mice"){test1 <- "mouse"}
if (test1 == "rats"){test1 <- "rat"}

str_view(test, regex("\\bcat(s|\\b)[^.]", ignore_case = T))
test[1:500]

test <- paste(species.temp.df$species.temp[1])




#per column
#Randomization
rob.comparison.randomization <- table(rob.mined.df$Randomization == rob.extracted.df$Randomization)
print("Randomization")
as.numeric(rob.comparison.randomization[names(rob.comparison.randomization)==T])/nrow(rob.mined.df)*100

#Blinding
rob.comparison.blinding <- table(rob.mined.df$Blinding == rob.extracted.df$Blinding)
print("Blinding")
as.numeric(rob.comparison.blinding[names(rob.comparison.blinding)==T])/nrow(rob.mined.df)*100

#Welfare
rob.comparison.welfare <- table(rob.mined.df$Welfare == rob.extracted.df$Welfare)
print("Welfare")
as.numeric(rob.comparison.welfare[names(rob.comparison.welfare)==T])/nrow(rob.mined.df)*100

#Conflict
rob.comparison.conflict <- table(rob.mined.df$Conflict == rob.extracted.df$Conflict)
print("Conflict")
as.numeric(rob.comparison.conflict[names(rob.comparison.conflict)==T])/nrow(rob.mined.df)*100

#Samplesize
rob.comparison.samplesize <- table(rob.mined.df$Samplesize == rob.extracted.df$Samplesize)
print("Samplesize")
as.numeric(rob.comparison.samplesize[names(rob.comparison.samplesize)==T])/nrow(rob.mined.df)*100

#ARRIVE
rob.comparison.arrive <- table(rob.mined.df$ARRIVE == rob.extracted.df$ARRIVE)
print("ARRIVE")
as.numeric(rob.comparison.arrive[names(rob.comparison.arrive)==T])/nrow(rob.mined.df)*100



temp.x.method
test <- read.pdf.EE(pdf.rob.vector[11])
length(test)




#################################################
############# function development ##############
#################################################

##Disease model
test <- read.pdf.EE(pdf.rob.vector[1])
test <- c("we used the following model for our study G93A-SOD1. ALso, we used cycad cycadcycad and cycad. In addition, we used G86R andalsoG86R")

ALSmodel_regex <- c("G93A", "G86R", "G85R", "G37R", "SOD186R",
                    "Tg-SOD1", "SOD1 Tg", "SOD1/rag2", "\\brag2", "SOD1 mice",
                    "TDP-43Q331K", "Q331K", "\\bFUS\\b",
                    "A315T", "Q331K", "rNLS8", "\\bNLS8\\b", 
                    "wobbler", "\\bwr", "\\bMeHg\\b", "\\bcycad")

#mine for the model and the number of appearances and show the first two most common hits
model.count.temp <- vector()
model.temp <- vector()
for(numCount in seq_along(ALSmodel_regex)) 
{
  print(str_count(test, regex(numCount)))
  model.count.temp[[numCount]] <- sum(str_count(test, regex(ALSmodel_regex[[numCount]], ignore_case = F)))#ignore case!
  model.temp[[numCount]] <- str_extract(paste(test, collapse = " "), regex(ALSmodel_regex[[numCount]], ignore_case = F))#ignore case!
}
model.temp.df <- as.data.frame(t(rbind(model.temp, model.count.temp)), stringsAsFactors = F)
model.temp.df$model.count.temp <- as.numeric(model.count.temp)
model.temp.df$model.temp <- as.character(model.temp)
model.temp.df <- model.temp.df[order(desc(model.temp.df$model.count.temp)),]

model1 <- paste(model.temp.df$model.temp[1])
model_count1 <- paste(model.temp.df$model.count.temp[1])
model2 <- paste(model.temp.df$model.temp[2])
model_count2 <- paste(model.temp.df$model.count.temp[2])

if(model1 == "NA"){model1 <- NA}
if(model2 == "NA"){model2 <- NA}

model1
model_count1
model2
model_count2

## Total numbers of used terms for ALS animal models
sum(str_count(test, regex(ALSmodel_regex)))

###Species
#alternative approach: order species according to most common appearance, seems better
test <- read.pdf.EE(pdf.rob.vector[1])
test <- c("we used the following model for our study G93A-SOD1 mice. ALso, we used cycad cycadcycad and cycad in mice. In addition, we used G86R andalsoG86R")


species_list <- c("\\brat[^ihe]", "\\bmarmoset", "\\bguinea pig", "\\bmacaque", "\\brhesus monkey",
                  "\\bswine", "\\bdog", "\\bcat", "\\brabbit", "\\bprimate",
                  "\\bhamster", "\\bsheep", "\\bmouse", "\\bmice")

species.count.temp <- vector()
species.temp <- vector()
for(numCount in seq_along(species_list)) 
{
  print(str_count(test, regex(numCount)))
  species.count.temp[[numCount]] <- sum(str_count(test, regex(species_list[[numCount]])))
  species.temp[[numCount]] <- str_extract(test, regex(species_list[[numCount]]))
}
species.temp.df <- as.data.frame(t(rbind(species.temp, species.count.temp)), stringsAsFactors = F)
species.temp.df$species.count.temp <- as.numeric(species.count.temp)
species.temp.df$species.temp <- as.character(species.temp)
species.temp.df <- species.temp.df[order(desc(species.temp.df$species.count.temp)),]
species.temp.df <- species.temp.df[!(species.temp.df$species.count.temp == 0),]

species1 <- paste(species.temp.df$species.temp[1])
species_count1 <- paste(species.temp.df$species.count.temp[1])
species2 <- paste(species.temp.df$species.temp[2])
species_count2 <- paste(species.temp.df$species.count.temp[2])

if(species1 == "NA"){species1 <- NA}
if(species2 == "NA"){species2 <- NA}

###Sex
test <- c("they used female animals and either species")

female_regex <- "\\bfemale(s|)\\b"
male_regex <- "\\bmale(s|)\\b"
sex_regex <- c("either sex", "both genders")

has_female <- str_subset(test, regex(female_regex, ignore_case = T))
has_female <- str_extract(test, regex(female_regex, ignore_case = T))
has_female <- paste(unlist(has_female), collapse = ";")
ifelse(str_detect(has_female, regex(female_regex, ignore_case = T)) == F, female <- 0, female <- "female")

has_male <- str_subset(test, regex(male_regex, ignore_case = T))
has_male <- str_extract(test, regex(male_regex, ignore_case = T))
has_male <- paste(unlist(has_male), collapse = ";")
ifelse(str_detect(has_male, regex(male_regex, ignore_case = T)) == F, male <- 0, male <- "male")


if(sum(str_detect(female, regex(female_regex, ignore_case = T))) > 0 &
   sum(str_detect(male, regex(male_regex, ignore_case = T))) == 0){
  sex.test <- "female"
}
if(sum(str_detect(female, regex(female_regex, ignore_case = T))) == 0 &
   sum(str_detect(male, regex(male_regex, ignore_case = T))) > 0){
  sex.test <- "male"
}
if(sum(str_detect(female, regex(female_regex, ignore_case = T))) > 0 &
   sum(str_detect(male, regex(male_regex, ignore_case = T))) > 0){
  sex.test <- "female and male"
}
if(sum(str_detect(female, regex(female_regex, ignore_case = T))) == 0 &
   sum(str_detect(male, regex(male_regex, ignore_case = T))) == 0){
  sex.test <- NA
}
if(sum(str_detect(test, regex(sex_regex, ignore_case = T))) > 0){
  sex.test <- "female and male"
}
sex.test

###Outcomes
test <- read.pdf.EE(pdf.rob.vector[1])
test <- c("we used magnetic resonance imaging and DWI and gait analysis")

histology_list <- c("histology",  "histological", "histopathology", "histopathologic",
                    "staining", "electron microscopy","Confocal Microscopy",
                    "immunohistochemistry", "immunehistochemistry", 
                    "immunostaining(s|)", "fluorescent staining(s|)", 
                    "immunofluorescence", "immunohistochemical",
                    "Nissl", "Hematoxylin-eosin", "Hematoxylin", "Eosin", "H&E")
behaviour_list <- c("DigiGait", "Treadmill", "kinematics", "grip strength", 
                    "rota(-|)rod", "motor performance", "grip test", "motor coordination",
                    "motor dysfunction", "muscular strength", "gripmeter", "hang wire", 
                    "behavio(u|)r score", "open field", "runtime", "pole test",
                    "catwalk", "gait analysis", "grid hanging", "pole climbing",
                    "climbing pole", "hanging endurance", "neurological score", 
                    "hanging wire", "swim test", "Morris water maze", "Morris maze",
                    "SHIRPA", "reflex test", "walking speed", "stair( |)case", 
                    "whisker behavio(u|)r")
imaging_list <- c("MR imaging", "magnetic resonance imaging", "\\bMRI\\b",
                  "\\bNODDI\\b", "diffusion tensor imaging", "\\bDTI\\b",
                  "diffusion( |-)weighted imaging", "\\bDWI\\b",
                  "diffusion kurtosis imaging")

has_histology <- str_extract(test, regex(histology_list, ignore_case = T))
has_histology <- paste(unlist(has_histology), collapse = ";")
ifelse(str_detect(has_histology, regex(histology_list, ignore_case = T)) == F, histology <- 0, histology <- "histology")

has_behaviour <- str_extract(test, regex(behaviour_list, ignore_case = T))
has_behaviour <- paste(unlist(has_behaviour), collapse = ";")
ifelse(str_detect(has_behaviour, regex(behaviour_list, ignore_case = T)) == F, behaviour <- 0, behaviour <- "behaviour")

has_imaging <- str_extract(test, regex(imaging_list, ignore_case = T))
has_imaging <- paste(unlist(has_imaging), collapse = ";")
ifelse(str_detect(has_imaging, regex(imaging_list, ignore_case = T)) == F, imaging <- 0, imaging <- "imaging")

outcome <- paste(c(histology,behaviour,imaging),collapse=", ")

if(sum(str_detect(outcome, regex("histology, behaviour, imaging", ignore_case = T))) > 0){
  outcome.test <- "histology, behaviour, imaging"
}
if(sum(str_detect(outcome, regex("histology, behaviour, 0", ignore_case = T))) > 0){
  outcome.test <- "histology, behaviour"
}
if(sum(str_detect(outcome, regex("histology, 0, imaging", ignore_case = T))) > 0){
  outcome.test <- "histology, imaging"
}
if(sum(str_detect(outcome, regex("0, behaviour, imaging", ignore_case = T))) > 0){
  outcome.test <- "behaviour, imaging"
}
if(sum(str_detect(outcome, regex("histology, 0, 0", ignore_case = T))) > 0){
  outcome.test <- "histology"
}
if(sum(str_detect(outcome, regex("0, behaviour, 0", ignore_case = T))) > 0){
  outcome.test <- "behaviour"
}
if(sum(str_detect(outcome, regex("0, 0, imaging", ignore_case = T))) > 0){
  outcome.test <- "imaging"
}
if(sum(str_detect(outcome, regex("0, 0, 0", ignore_case = T))) > 0){
  outcome.test <- "NA"
}
outcome.test

###DOI
test <- "A systematic review with meta-analysis on the antihypertensive efficacy of Nigerian medicinal plants,M. A. Abdulazeez S. A. Muhammad Y. Saidu A. B. Sallau A. A. Arzai M. A. Tabari A. Hafiz M. Y. Gwarzo J. Manosroi A. Idi M. Bashir S. L. Pedro,Journal of Ethnopharmacology,Journal of Ethnopharmacology,\"ETHNOPHARMACOLOGICAL RELEVANCE: Despite the promising effects of herbal preparations in lowering blood pressure (BP), hypertension remains a major clinical challenge in Nigeria. The BP-lowering effects of medicinal plants are due to the presence of bioactive compounds. AIM OF THE STUDY: This meta-analysis presents a precise estimate of the therapeutic benefits of medicinal plants utilized in Nigeria for the management of hypertension in animals and humans. METHODS: A systematic literature search was performed through Cochrane, PubMed, Science Direct and Scopus databases from inception until February 28, 2021 using search terms related to randomized controlled trials of Nigerian medicinal plants for hypertension. Additional studies were identified through manual search. BP was the main outcome that was measured after the intervention. Meta-analysis was performed using the Review Manager and Meta-Essential. RESULTS: Nineteen trials comprising of 16 preclinical and 3 clinical studies were enrolled for the meta-analysis. A total number of 16 plants was identified of which H. sabdariffa was the highest reported plant. The plant extracts significantly lowered the systolic blood pressure (SBP) and diastolic blood pressure (DBP) of the hypertensive subjects compared to control. Weighted mean difference (WMD) for SBP (-43.60 mmHg, 95% CI: -63.18, -24.01 p<0.0001) and DBP (-29.50 mmHg, 95 CI: -43.66, -15.34 p<0.0001) was observed for the preclinical studies. For clinical trials, the WMD was -13.98 mmHg, 95 CI: -19.08, -8.88 p<0.00001 for SBP and -10.00 mmHg, 95 CI: -12.22, -7.78 p<0.00001 for DBP. High heterogeneity was observed for the outcome measures of preclinical studies, but not for the clinical studies. The observed substantial heterogeneity in preclinical studies may be linked to methodological shortcomings as evidenced by the results of the risk of bias assessment. There was no evidence of publication bias in animal trials for BP using the funnel plot and Egger's regression test (SBP, p=0.239 and DBP, p=0.112). CONCLUSIONS: This study provides evidence of medicinal preparations for the treatment of hypertension. A well-conducted trial with methodological rigour and a longer duration of follow-up is required for their effective clinical utilization.\",https://ovidsp.ovid.com/ovidweb.cgi?T=JS&CSC=Y&NEWS=N&PAGE=fulltext&D=med20&AN=34157327,,2021,10.1016/j.jep.2021.114342,Journal Article,\"Animals Antihypertensive Agents/ip [Isolation & Purification] *Antihypertensive Agents/pd [Pharmacology] Blood Pressure/de [Drug Effects] Humans *Hypertension/dt [Drug Therapy] Nigeria *Plant Extracts/pd [Pharmacology] Plants, Medicinal/ch [Chemistry] Randomized Controlled Trials as Topic Research Design 0 (Antihypertensive Agents) 0 (Plant Extracts)\",, NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA"
test <- "10.1016/j.jep.2021.114342 and 10.1016/j.jclinepi.2022.05.019 test experiment"
test <- " test experiment does not contain a doi classifier"

str_detect(test, regex("10.\\d{4,9}/[-._;()/:a-z0-9A-Z]+", ignore_case = T))


if(str_detect(test, regex("10.\\d{4,9}/[-._;()/:a-z0-9A-Z]+", ignore_case = T)) == T){
  doi <- str_extract_all(test, regex("10.\\d{4,9}/[-._;()/:a-z0-9A-Z]+", ignore_case = T))
}else{
  doi <- NA
}
doi[[1]][1]






#converts spelled number words to actual numbers
word2num <- function(word){
  wsplit <- strsplit(tolower(word)," ")[[1]]
  one_digits <- list(zero=0, one=1, two=2, three=3, four=4, five=5,
                     six=6, seven=7, eight=8, nine=9)
  teens <- list(eleven=11, twelve=12, thirteen=13, fourteen=14, fifteen=15,
                sixteen=16, seventeen=17, eighteen=18, nineteen=19)
  ten_digits <- list(ten=10, twenty=20, thirty=30, forty=40, fifty=50,
                     sixty=60, seventy=70, eighty=80, ninety=90)
  doubles <- c(teens,ten_digits)
  out <- 0
  i <- 1
  while(i <= length(wsplit)){
    j <- 1
    if(i==1 && wsplit[i]=="hundred")
      temp <- 100
    else if(i==1 && wsplit[i]=="thousand")
      temp <- 1000
    else if(wsplit[i] %in% names(one_digits))
      temp <- as.numeric(one_digits[wsplit[i]])
    else if(wsplit[i] %in% names(teens))
      temp <- as.numeric(teens[wsplit[i]])
    else if(wsplit[i] %in% names(ten_digits))
      temp <- (as.numeric(ten_digits[wsplit[i]]))
    if(i < length(wsplit) && wsplit[i+1]=="hundred"){
      if(i>1 && wsplit[i-1] %in% c("hundred","thousand"))
        out <- out + 100*temp
      else
        out <- 100*(out + temp)
      j <- 2
    }
    else if(i < length(wsplit) && wsplit[i+1]=="thousand"){
      if(i>1 && wsplit[i-1] %in% c("hundred","thousand"))
        out <- out + 1000*temp
      else
        out <- 1000*(out + temp)
      j <- 2
    }
    else if(i < length(wsplit) && wsplit[i+1] %in% names(doubles)){
      temp <- temp*100
      out <- out + temp
    }
    else{
      out <- out + temp
    }
    i <- i + j
  }
  return(list(word,out))
}



######### define path where PDFs are ##################

pdf.test.vector <- list.files("C:/Users/benja/Dropbox/Data/MS/Experiments/Karolinska/ALS_SystematicReview/Animal_study/Rob_test1_sample", full.names = T)
pdf.test.vector <- pdf.test.vector[grep("pdf", pdf.test.vector, ignore.case = T)]
pdf.test.vector <- pdf.test.vector[!grepl("pdftools", pdf.test.vector)]
pdf.test.vector


######## predefine strings for searches ############
#Species
#careful - add at the end otherwise you always need to change all possible entries in the script

######## data mining function ######################

pdf.extraction(pdf.test.vector[pdf.i])

src = pdf.test.vector[pdf.i]
pdf.extraction <- function(src)
{
  temp.x <- read.pdf.EE(src = src)
  
  #this creates a vector with 1 line as a element
  
  ############################################
  ###################### Define paper sections
  
  #Abstract
  
  if(length(grep("Abstract|a b s t r a c t", temp.x, ignore.case = T)) == 0)
  {
    Abstract.start <- 1
  }else{
    Abstract.start <- grep("Abstract|a b s t r a c t", temp.x, ignore.case = T)  
  }
  
  #Introduction
  
  if(length(grep("Introduction", temp.x, ignore.case = T)) == 0)
  {
    Introduction.start <- 1
  }else{
    Introduction.start <- grep("Introduction", temp.x, ignore.case = T)  
  }
  
  #if(length(grep("Introduction", temp.x, ignore.case = T)) == 0)
  #{
  #  Introduction.start <- 1
  #}
  #if(length(grep("Introduction", temp.x, ignore.case = T)) > 0)
  #{
  #  Introduction.start <- min(grep("Introduction", temp.x, ignore.case = T))
  #  
  #}
  #if(length(grep("Background", temp.x, ignore.case = T)) > 0)#paper [123] has 5 times "background" ??? but length output = 1
  #{
  #  Introduction.start <- min(grep("Background", temp.x, ignore.case = T)  )
  #}
  
  #Methods
  
  #Method.start <- min(c(grep("methods", temp.x, ignore.case = T)[grep("methods", temp.x, ignore.case = T) > 4]),
  #                    grep("Experimental procedures", temp.x, ignore.case = T)[grep("Experimental procedures", temp.x, ignore.case = T) > 4])
  
  if(length(grep("materials and methods", temp.x, ignore.case = T)) == 0)
  {
    if(length(grep("experimental procedures", temp.x, ignore.case = T)[grep("experimental procedures", temp.x, ignore.case = T) > 33]) == 0)
    {
      if(length(grep("methods", temp.x, ignore.case = T)[grep("methods", temp.x, ignore.case = T) > 33]) == 0)
      {Method.start <- 33}
      else{
        Method.start <- min(grep("methods", temp.x, ignore.case = T)[grep("methods", temp.x, ignore.case = T) > 33])
      }
    }else{
      Method.start <- min(grep("experimental procedures", temp.x, ignore.case = T)[grep("experimental procedures", temp.x, ignore.case = T) > 33])
    }    
  }else{
    Method.start <- min(grep("materials and methods", temp.x, ignore.case = T))
  }
  
  
  #if(length(grep("materials and methods", temp.x, ignore.case = T)[grep("materials and methods", temp.x, ignore.case = T) > 33]) == 0)
  #{
  #  Method.start <- 33#according to median start of Introduction from 240 sample papers (an earlier start would risk that methods appears in a structured abstract)
  #}
  #if(length(grep("methods", temp.x, ignore.case = T)[grep("methods", temp.x, ignore.case = T) > 33]) > 0)
  #{
  #  Method.start <- min(grep("methods", temp.x, ignore.case = T)[grep("methods", temp.x, ignore.case = T) > 33])
  #}
  #if(length(grep("Experimental procedures", temp.x, ignore.case = T)[grep("Experimental procedures", temp.x, ignore.case = T) > 33]) > 0)
  #{
  #  Method.start <- min(grep("Experimental procedures", temp.x, ignore.case = T)[grep("Experimental procedures", temp.x, ignore.case = T) > 33])
  #}
  #works well but you need to set a minimum term to exclude the word methods in structured abstracts
  
  #Results --> most difficult to define results start since "result" appears commonly in scientific texts
  if(length(grep("Results", temp.x, ignore.case = T)[grep("Results", temp.x, ignore.case = T) > 33]) == 0)
  {
    Results.start <- 33
  }
  if(length(grep("Results", temp.x, ignore.case = T)[grep("Results", temp.x, ignore.case = T) > 33]) > 0)
  {
    Results.start <- min(grep("Results", temp.x, ignore.case = T)[grep("Results", temp.x, ignore.case = T) > 33])
  }
  
  #Results.start <- (grep("Results", temp.x, ignore.case = T)[grep("Results", temp.x, ignore.case = T) > 33])#any appearance of "results" in the abstract should be avoided
  #Results.start <- min(Results.start[Results.start>Method.start])
  #get the first "results" after the start of the methods
  
  #Discussion
  
  if(length(grep("Discussion", temp.x, ignore.case = T)[grep("Discussion", temp.x, ignore.case = T) > 33]) == 0)
  {
    if(length(grep("Results", temp.x, ignore.case = T)[grep("Results", temp.x, ignore.case = T) > 33]) == 0)
    {
      discussion.start <- 33
    }else{
      discussion.start <-  min(grep("Results", temp.x, ignore.case = T)[grep("Results", temp.x, ignore.case = T) > 33])
    }
  }else{
    discussion.start <- min(grep("Discussion", temp.x, ignore.case = T)) 
  }
  
  #if(length(grep("Conclusion" OR "Conclusions", temp.x, ignore.case = T)) == 0)
  #{
  #  if(sum(grep("Discussion", temp.x, ignore.case = T)) == 0)
  #  {
  #    discussion.start <- min(grep("Results", temp.x, ignore.case = T))
  #  }else{
  #    discussion.start <-  min(grep("Discussion", temp.x, ignore.case = T))
  #  }
  #}else {
  #  discussion.start <- min(c(grep("Discussion", temp.x, ignore.case = T), grep("Conclusion", temp.x, ignore.case = T)))   
  #}
  
  #References
  
  if(length(grep("References", temp.x, ignore.case = T)[grep("References", temp.x, ignore.case = T) > 33]) == 0)
  {
    References.start <- length(temp.x)
  }else{
    References.start <- min(grep("References", temp.x, ignore.case = T)[grep("References", temp.x, ignore.case = T) > 33])
  }
  
  
  ############################################
  ######################### Map paper sections 
  
  #Abstract
  
  if(Abstract.start == 1)
  {
    if(length(grep("Introduction", temp.x, ignore.case = T)) == 0){temp.x.abstract <- temp.x[1:80]}
    else{temp.x.abstract <- temp.x[1:Introduction.start]}}
  if(Abstract.start > 1)
  {
    if(length(grep("Introduction", temp.x, ignore.case = T)) == 0){temp.x.abstract <- temp.x[Abstract.start:80]}
    else{temp.x.abstract <- temp.x[Abstract.start:Introduction.start]}}
  
  temp.x.abstract <- paste(unlist(temp.x.abstract), sep = " ")
  temp.x.abstract <- paste(temp.x.abstract, collapse = " ")
  
  #if(Introduction.start == 1)
  #{
  #  temp.x.abstract <- temp.x[1:Method.start]
  #}
  #if(Introduction.start > 1)
  #{
  #  temp.x.abstract <- temp.x[Introduction.start]  
  #}
  
  #}else{
  #  if(Introduction.start < Method.start){intro.range <- Introduction.start:Method.start}
  #  if(Results.start < Method.start){intro.range <- Introduction.start:Results.start}
  #  temp.x.abstract <- temp.x[1:Introduction.start]
  #}
  
  #Introduction
  
  if(Introduction.start < Method.start){intro.range <- Introduction.start:Method.start
  }else{intro.range <- Introduction.start+80}#which value covers 90% of all Intros? 80 is a raw estimate
  
  temp.x.intro <- temp.x[intro.range]
  temp.x.intro <- paste(unlist(temp.x.intro), sep = " ")
  temp.x.intro <- paste(temp.x.intro, collapse = " ")
  
  #if(Results.start < Method.start){intro.range <- Introduction.start:Results.start}#difficult since results start is hard to define using simple text mining
  
  
  #Experimental (in most cases, you want to search methods and results for our terms, experimental = methods + results)
  
  if(discussion.start < Method.start)
  {
    experimental.range <- Method.start:References.start
  }else{
    experimental.range <- Method.start:discussion.start
  }
  
  temp.x.method <- temp.x[experimental.range]
  temp.x.method <- paste(unlist(temp.x.method), sep = " ")
  temp.x.method <- paste(temp.x.method, collapse = " ")
  
  #if(Method.start < Results.start){range <- Method.start:Results.start}#difficult since results start is hard to define using simple text mining
  
  #Whole paper w/o References
  paper.range <- Abstract.start:References.start
  temp.x.paper <- temp.x[paper.range]
  temp.x.paper <- paste(unlist(temp.x.paper), sep = " ")
  temp.x.paper <- paste(temp.x.paper, collapse = " ")
  
  ##############################################
  ########################## Regular expressions
  #file title
  
  file.title.temp <- str_extract(pdf.test.vector[pdf.i], regex("/.+", ignore_case = T))
  file.title <- str_extract(file.title.temp, regex("[^/].+"))
  
  
  #journal
  #journal.temp <- unique (grep(paste(Journal_vector, collapse="|"),
  #                             temp.x.abstract, ignore.case = T, value = T))##should work but it does not...
  
  #if(sum(grepl("journal", temp.x.abstract, ignore.case = T))==0)
  #{
  #  journal <- NA
  #}else{
  #  journal.temp <- temp.x.abstract[min(grep("journal",temp.x.abstract, ignore.case = T))]
  #  journal.temp <- strsplit2(journal.temp, split=(","))
  #  journal <- journal.temp[1,grep("journal", journal.temp, ignore.case = T)]
  #  
  #}
  
  #journal.temp <- strsplit2(journal.temp, split=(". "))
  
  
  #author
  
  if(str_detect(pdf.test.vector[pdf.i], regex("/[a-z| ]+-", ignore_case = T)) == F){
    author <- NA
  }else{
    author.temp <- str_extract(pdf.test.vector[pdf.i], regex("/[a-z| ]+-", ignore_case = T))
    author <- str_extract(author.temp, regex("[a-z]+", ignore_case = T))
  }
  
  
  #year
  #if(sum(c(grepl("submitted",temp.x, ignore.case = T),
  #         grepl("Published",temp.x, ignore.case = T),  
  #         grepl("accepted",temp.x, ignore.case = T))) == 0)
  #{
  #  year.temp <- t(strsplit2(temp.x.abstract, " "))
  #  year.temp <- gsub(";", "", year.temp)
  #  year.temp <- gsub(",", "", year.temp)
  
  #  year.temp <- as.numeric(year.temp)
  #  year.temp <- year.temp[!is.na(year.temp)]
  #  year.temp <- year.temp[nchar(year.temp) == 4]
  
  # }else{
  #  year.temp <- temp.x[min(c(grep("submitted",temp.x, ignore.case = T),
  #                           grep("Published",temp.x, ignore.case = T),  
  #                            grep("accepted",temp.x, ignore.case = T)))]
  #  year.temp <- t(strsplit2(year.temp, " "))
  #  year.temp <- gsub(";", "", year.temp)
  #  year.temp <- as.numeric(year.temp)
  #  year.temp <- year.temp[!is.na(year.temp)]
  #  year <- max(year.temp)
  #}
  
  if(str_detect(pdf.test.vector[pdf.i], regex("[^\\d| ]\\d\\d\\d\\d[^\\d| ]", ignore_case = T)) == T){
    year.temp <- str_extract(pdf.test.vector[pdf.i], regex("[^\\d| ]\\d\\d\\d\\d[^\\d| ]", ignore_case = T))
    year.temp <- str_extract(year.temp, regex("\\d\\d\\d\\d", ignore_case = T))
    year <- ifelse(year.temp < 2019 & year.temp > 1900, year.temp, NA)
  }else{
    year.temp <- t(strsplit2(temp.x.abstract, " "))
    year.temp <- gsub(";", "", year.temp)
    year.temp <- gsub(",", "", year.temp)
    
    year.temp <- as.numeric(year.temp)
    year.temp <- year.temp[!is.na(year.temp)]
    year.temp <- year.temp[nchar(year.temp) == 4]
    year.temp <- year.temp[ year.temp < 2019]
    year.temp <- year.temp[ year.temp > 1900]
    year.temp <- max(year.temp, na.rm = T)
    year <- ifelse(year.temp < 2019 & year.temp > 1900, year.temp, NA)
  }
  
  #paper title
  if(str_detect(pdf.test.vector[pdf.i], regex("-[a-z| ]+\\.", ignore_case = T)) == F){
    paper.title <- NA
  }else{
    paper.title.temp <- str_extract(pdf.test.vector[pdf.i], regex("-[a-z| ]+\\.", ignore_case = T))
    paper.title <- str_extract(paper.title.temp, regex("[a-z| ]+", ignore_case = T))
  }
  
  
  
  #doi
  if(sum(grepl("doi:", temp.x.paper))==0){
    doi <- NA
  }else{
    doi.temp <- temp.x.paper[grep("doi:",temp.x.paper)]
    doi.temp <- t(strsplit2(doi.temp, " "))
    doi <- paste(doi.temp[grep("doi:", doi.temp)],collapse=";")
  }
  
  
  #E-mail
  #if(sum(grepl("@", email.temp))==0){
  #  email <- NA
  #}else{
  #  email.temp <- temp.x[grep("@",temp.x)]
  #  email.temp <- t(strsplit2(email.temp, " "))
  #  email <- paste(email.temp[grep("@", email.temp)],collapse=";")
  #}
  
  if(str_detect(toString(temp.x.paper), regex("([_a-z0-9-]+(\\.[_a-z0-9-]+)*@[a-z0-9-]+(\\.[a-z0-9-]+)*(\\.[a-z]{2,4}))", ignore_case = T)) == F){
    email <- NA
  }else{
    email.temp <- str_extract_all(temp.x.paper, regex("([_a-z0-9-]+(\\.[_a-z0-9-]+)*@[a-z0-9-]+(\\.[a-z0-9-]+)*(\\.[a-z]{2,4}))", ignore_case = T))
    email.temp <- Filter(length, email.temp)
    email <- paste(email.temp, collapse = "; ")}
  
  
  #Species
  species.count.temp <- vector()
  species.temp <- vector()
  for(numCount in seq_along(species_list)) 
  {
    species.count.temp[[numCount]] <- sum(str_count(temp.x.method, regex(species_list[[numCount]])))
    species.temp[[numCount]] <- str_extract(temp.x.method, regex(species_list[[numCount]]))
  }
  species.temp.df <- as.data.frame(t(rbind(species.temp, species.count.temp)), stringsAsFactors = F)
  species.temp.df$species.count.temp <- as.numeric(species.count.temp)
  species.temp.df$species.temp <- as.character(species.temp)
  species.temp.df <- species.temp.df[order(desc(species.temp.df$species.count.temp)),]
  species.temp.df <- species.temp.df[!(species.temp.df$species.count.temp == 0),]
  
  species1 <- paste(species.temp.df$species.temp[1])
  species_count1 <- paste(species.temp.df$species.count.temp[1])
  species2 <- paste(species.temp.df$species.temp[2])
  species_count2 <- paste(species.temp.df$species.count.temp[2])
  
  if(species1 == "NA"){species1 <- NA}
  if(species2 == "NA"){species2 <- NA}
  
  #species.temp <- str_subset(temp.x.method, species_match)
  #species.temp <- str_extract_all(temp.x.method, regex(species_match, ignore_case = T))
  #species.temp <- Filter(length, species.temp)
  #species.temp <- paste(unlist(species.temp), collapse = ";")
  #species.temp <- gsub(" ", "", species.temp)
  #species.temp <- str_to_lower(species.temp)
  #species.temp <- str_extract_all(species.temp, regex(species_match2, ignore_case = T))
  #species.temp <- unique(unlist(species.temp))
  #if(sum(str_detect(species.temp, regex("mouse", ignore_case = T))) > 0 &
  #   sum(str_detect(species.temp, regex("mice", ignore_case = T))) > 0){
  #  species.temp <- gsub("mice", "mouse", species.temp)
  #  species.temp <- unique(unlist(species.temp))
  #}
  #if(sum(str_detect(species.temp, regex("mice", ignore_case = T))) > 0){
  #  species.temp <- gsub("mice", "mouse", species.temp)
  #}
  #species <- paste(unlist(species.temp), collapse = ";")
  
  #Sex
  has_female <- str_subset(temp.x.method, regex(female_regex, ignore_case = T))
  has_female <- str_extract(temp.x.method, regex(female_regex, ignore_case = T))
  has_female <- paste(unlist(has_female), collapse = ";")
  ifelse(str_detect(has_female, regex(female_regex, ignore_case = T)) == F, female <- 0, female <- "female")
  
  has_male <- str_subset(temp.x.method, regex(male_regex, ignore_case = T))
  has_male <- str_extract(temp.x.method, regex(male_regex, ignore_case = T))
  has_male <- paste(unlist(has_male), collapse = ";")
  ifelse(str_detect(has_male, regex(male_regex, ignore_case = T)) == F, male <- 0, male <- "male")
  
  if(sum(str_detect(female, regex("female", ignore_case = T))) > 0 &
     sum(str_detect(male, regex("male", ignore_case = T))) == 0){
    sex.test <- "female"
  }
  if(sum(str_detect(female, regex("female", ignore_case = T))) == 0 &
     sum(str_detect(male, regex("male", ignore_case = T))) > 0){
    sex.test <- "male"
  }
  if(sum(str_detect(female, regex("female", ignore_case = T))) > 0 &
     sum(str_detect(male, regex("male", ignore_case = T))) > 0){
    sex.test <- "female and male"
  }
  if(sum(str_detect(female, regex("female", ignore_case = T))) == 0 &
     sum(str_detect(male, regex("male", ignore_case = T))) == 0){
    sex.test <- NA
  }
  
  
  #Model disease
  model.count.temp <- vector()
  model.temp <- vector()
  for(numCount in seq_along(model_regex)) 
  {
    model.count.temp[[numCount]] <- sum(str_count(temp.x.method, regex(model_regex[[numCount]], ignore_case = F)))#ignore case!
    model.temp[[numCount]] <- str_extract(temp.x.method, regex(model_regex[[numCount]], ignore_case = F))#ignore case!
  }
  model.temp.df <- as.data.frame(t(rbind(model.temp, model.count.temp)), stringsAsFactors = F)
  model.temp.df$model.count.temp <- as.numeric(model.count.temp)
  model.temp.df$model.temp <- as.character(model.temp)
  model.temp.df <- model.temp.df[order(desc(model.temp.df$model.count.temp)),]
  
  model1 <- paste(model.temp.df$model.temp[1])
  model_count1 <- paste(model.temp.df$model.count.temp[1])
  model2 <- paste(model.temp.df$model.temp[2])
  model_count2 <- paste(model.temp.df$model.count.temp[2])
  
  if(model1 == "NA"){model1 <- NA}
  if(model2 == "NA"){model2 <- NA}
  
  
  #Disease
  disease.count.temp <- vector()
  disease.temp <- vector()
  for(numCount in seq_along(disease_regex)) 
  {
    disease.count.temp[[numCount]] <- sum(str_count(temp.x.intro, regex(disease_regex[[numCount]], ignore_case = F)))
    disease.temp[[numCount]] <- str_extract(temp.x.intro, regex(disease_regex[[numCount]], ignore_case = F))
  }
  disease.temp.df <- as.data.frame(t(rbind(disease.temp, disease.count.temp)), stringsAsFactors = F)
  disease.temp.df$disease.count.temp <- as.numeric(disease.count.temp)
  disease.temp.df$disease.temp <- as.character(disease.temp)
  disease.temp.df <- disease.temp.df[order(desc(disease.temp.df$disease.count.temp)),]
  
  disease1 <- paste(disease.temp.df$disease.temp[1])
  disease_count1 <- paste(disease.temp.df$disease.count.temp[1])
  disease2 <- paste(disease.temp.df$disease.temp[2])
  disease_count2 <- paste(disease.temp.df$disease.count.temp[2])
  
  if(disease1 == "NA"){disease1 <- NA}
  if(disease2 == "NA"){disease2 <- NA}
  
  
  #Animal numbers
  if(str_detect(toString(temp.x.method), regex(animal_number_regex_total, ignore_case = T)) == T){
    if(str_detect(toString(temp.x.method), regex("\\d+", ignore_case = T)) == T){
      animal_number.temp <- str_extract_all(temp.x.method, regex(animal_number_regex_total, ignore_case = T))
      animal_number.temp <- Filter(length, animal_number.temp)
      animal_number.temp <- toString(animal_number.temp)
      animal_number.temp <- str_extract_all(animal_number.temp, regex("\\d+", ignore_case = T))
      animal_number <- sapply(animal_number.temp, function(x) sum(as.numeric(x)))
    }else{
      spelled_match.temp <- c()
      animal_number.temp <- c()
      for(numCount.numbers in seq_along(numbers_spelled)) 
      {for(numCount.animals in seq_along(animals))
      {spelled_match.temp <- str_c(numbers_spelled[[numCount.numbers]], animals[[numCount.animals]])
      print(spelled_match.temp)
      }
        animal_number.temp[[numCount.numbers]] <- str_extract_all(temp.x.method, regex(spelled_match.temp, ignore_case = T))
      }
      animal_number.temp <- Filter(length, animal_number.temp)
      animal_number.temp <- toString(animal_number.temp)
      animal_number.temp <- str_extract_all(animal_number.temp, regex(animal_number_regex_spelled, ignore_case = T))
      animal_number.temp <- unlist(animal_number.temp, use.names=FALSE)
      
      animal_number.temp2 <- c()
      for(numCount in seq_along(animal_number.temp)) 
      {
        animal_number.temp2[[numCount]] <- animal_number.temp[[numCount]] %>% word2num %>% .[[2]]
      }
      animal_number <- sum(animal_number.temp2)
    }
  }else{
    if(str_detect(toString(temp.x.method), regex(animal_number_regex_n, ignore_case = T)) == T){
      animal_number.temp <- str_extract_all(temp.x.method, regex(animal_number_regex_n, ignore_case = T))
      animal_number.temp <- Filter(length, animal_number.temp)
      animal_number.temp <- toString(animal_number.temp)
      animal_number.temp <- str_extract_all(animal_number.temp, "\\d+")
      animal_number <- sapply(animal_number.temp, function(x) sum(as.numeric(x)))
      
    }else{animal_number <- NA}}
  
  
  #Risk of bias assessment
  #Welfare
  if(str_detect(temp.x.method, regex("approved|approval|carried out in accordance|conducted in accordance", ignore_case = T)) == T){
    welfare <- "yes"
  }else{
    welfare <- "not reported"
  }

  #Blinding
  if(str_detect(temp.x.method, regex("blinded|were blind", ignore_case = T)) == T){
    blinding <- "yes"
  }else{
    blinding <- "not reported"
  }
  
  #Randomization
  if(str_detect(temp.x.method, regex("random", ignore_case = T)) == T){
    randomization <- "yes"
  }else{
    randomization <- "not reported"
  }
  
  #sample size calculation
  if(str_detect(temp.x.method, regex("achieve.{1,10}power|sample size calculation|power analysis", ignore_case = T)) == T){
    powertest <- "yes"
  }else{
    powertest <- "not reported"
  }
  
  #ARRIVE guidelines
  if(str_detect(temp.x.paper, regex("ARRIVE", ignore_case = F)) == T){#ignore case!
    ARRIVE <- "yes"
  }else{
    ARRIVE <- "not reported"
  }
  
  #conflict of interest
  if(str_detect(temp.x.paper, regex("conflict[s| ] ?of interest|disclosures|competing interests|nothing to disclose", ignore_case = T)) == T){
    disclosure <- "yes"
  }else{
    disclosure <- "not reported"
  }
  

  #OUTPUT
  return.x <- c(file.title, author, year, paper.title, 
                sex.test, species1, species_count1,
                animal_number, disease1, disease_count1, disease2, disease_count2,
                model1, model_count1, model2, model_count2,
                randomization, blinding, welfare, disclosure, powertest, ARRIVE,
                doi, email)
  names(return.x) <- c("File title", "First author", "Year","Paper title", 
                       "Sex", "Species1", "SCount1",
                       "Total animal count", "Disease1", "DCount1", "Disease2", "DCount2",
                       "Model1", "Mcount1", "Model2", "Mcount2",
                       "Random", "Blind", "Welfare statement", "Conflict of interest", "Sample size", "ARRIVE",
                       "doi", "email")
  return(return.x)
  
  
}




extraction.ee <- vector()

for(pdf.i in 1:length(pdf.test.vector)) 
{
  print(pdf.i)
  print(pdf.test.vector[pdf.i])
  extraction.ee <-  rbind(extraction.ee, pdf.extraction(src = pdf.test.vector[pdf.i]))
}

write_xlsx(as.data.frame(extraction.ee), "F:\\MS\\Experiments\\Neuroscience Translation\\extraction.xlsx")


extraction.ee

read.pdf.EE(pdf.test.vector[225])

pdf.extraction(src = pdf.test.vector[1])
pdf.extraction(src = pdf.test.vector[2])
pdf.extraction(src = pdf.test.vector[3])
pdf.extraction(src = pdf.test.vector[4])
pdf.extraction(src = pdf.test.vector[5])


######### Reference table for comparison ##############

Extraction <- read_excel("data_mining_1/Sample.xlsx")
Extraction <- as.data.frame(Extraction)
head(Extraction)



#### paper structure statistics

#which median line number starts introduction
src = pdf.test.vector[pdf.i]
pdf.format <- function(src)
{
  temp.x <- read.pdf.EE(src = src)
  if(length(grep("Introduction", temp.x, ignore.case = T)) == 0)
  {
    NA
  }else{
    min(grep("Introduction", temp.x, ignore.case = T))
  }}

paper_format <- vector()
for(pdf.i in 1:length(pdf.test.vector)) 
{
  print(pdf.i)
  print(pdf.test.vector[pdf.i])
  paper_format <-  rbind(paper_format, pdf.format(src = pdf.test.vector[pdf.i]))
}

median(paper_format, na.rm = T)
plot(paper_format)








#####definition of search strings

#paper
Journal_list <- "Search_files/journal_list/J_Medline.txt"#this approach works better
Journal_list <- Journal_list %>%  readChar(file.info(Journal_list)$size)
Journal_list <- str_split(Journal_list, "\n")
Journal_list <- paste(unlist(Journal_list), sep = " ")
Journal_list <- paste(Journal_list, collapse = " ")

#full journal title
#str_view_all(Journal_list, regex("JournalTitle.{0,150}MedAbbr?", ignore_case = T))
Journal_list_fulltitle <- str_extract_all(Journal_list, regex("JournalTitle.{0,150}MedAbbr?", ignore_case = T))
Journal_list_fulltitle <- paste(unlist(Journal_list_fulltitle), sep = " ")
Journal_list_fulltitle <- paste(Journal_list_fulltitle, collapse = " ")
Journal_list_fulltitle <- str_replace_all(Journal_list_fulltitle, " MedAbbr JournalTitle: ", "|")
Journal_list_fulltitle <- str_replace_all(Journal_list_fulltitle, "\\[.{1,25}\\] ", "")
Journal_list_fulltitle <- str_replace_all(Journal_list_fulltitle, "JournalTitle:", "")
Journal_list_fulltitle <- str_replace_all(Journal_list_fulltitle, " MedAbbr", "")
Journal_list_fulltitle <- str_replace_all(Journal_list_fulltitle, "[[:punct:]]", "")
#str_view_all(test, regex(Journal_list_fulltitle, ignore_case = T))

#journal abbreviation
#str_view_all(Journal_list, regex("MedAbbr:.{0,150}ISSN \\(Print\\)?", ignore_case = T))
Journal_list_abbreviation <- str_extract_all(Journal_list, regex("MedAbbr:.{0,150}ISSN \\(Print\\)?", ignore_case = T))
Journal_list_abbreviation <- paste(unlist(Journal_list_abbreviation), sep = " ")
Journal_list_abbreviation <- paste(Journal_list_abbreviation, collapse = " ")
Journal_list_abbreviation <- str_replace_all(Journal_list_abbreviation, " ISSN \\(Print\\) MedAbbr: ", "|")
Journal_list_abbreviation <- str_replace_all(Journal_list_abbreviation, "\\(.{1,20}\\)", "")
Journal_list_abbreviation <- str_replace_all(Journal_list_abbreviation, " ISSN ", "")
Journal_list_abbreviation <- str_replace_all(Journal_list_abbreviation, "MedAbbr: ", "")
Journal_list_abbreviation <- str_replace_all(Journal_list_abbreviation, "[[:punct:]]", "")
#str_view_all(test, regex(Journal_list_abbreviation, ignore_case = T))

#merge full and abbreviation list
Journal_list_regex <- str_c(Journal_list_fulltitle, Journal_list_abbreviation, sep = "|")

#str_view_all(test, regex(test_regex, ignore_case = T)) #Problem: short and unspecific journal titles are a problem, eg. age, brain, Dent, tic, rn - see line below
#str_extract_all(Journal_list_regex, regex("\\btic\\b", ignore_case = T, match = T))

test <- c("check out Family practice research journal and Proceedings, wowi Am J Econ Sociol, and Simul Synth Med Imagingopa, Annual Management Conference - American Dental Association, it is great")
test_regex <- c("journal of ribi", "Proceedings", "Am J Econ Sociol", "and my mother")
test_regex <- str_c(test_regex, collapse = "|")



###fetch html from webpage using the doi

#doi <- "doi:10.1016/j.jneuroim.2010.08.022"
doi <- "doi:10.1371/journal.pone.0140238"
doi <- substring(doi, 5)
doi <- str_c("https://doi.org/", doi)


paper.html.temp <- read_html(doi)
paper.html.temp <- html_text(paper.html.temp)
paper.html.temp <- strsplit2(paper.html.temp, split="\n")
str_c(paper.html.temp, collapse = ", ")


if(length(grep("Abstract", paper.html.temp, ignore.case = T)) == 0)
{
  Abstract.start <- 1
}else{
  Abstract.start <- grep("Abstract", paper.html.temp, ignore.case = T)  
}

#paper.html.temp <- doi

#paper.html.temp <- 'https://doi.org/10.1371/journal.pone.0140238'


if(length(grep("Introduction", paper.html.temp, ignore.case = T)) == 0)
{
  Introduction.start <- 1
}else{
  Introduction.start <- grep("Introduction", paper.html.temp, ignore.case = T)  
}

if(length(grep("materials and methods", paper.html.temp, ignore.case = T)) == 0)
{
  if(length(grep("experimental procedures", paper.html.temp, ignore.case = T)) == 0)
  {
    if(length(grep("methods", paper.html.temp, ignore.case = T)) == 0)
    {Method.start <- 33}
    else{
      Method.start <- min(grep("methods", paper.html.temp, ignore.case = T))
    }
  }else{
    Method.start <- min(grep("experimental procedures", paper.html.temp, ignore.case = T))
  }    
}else{
  Method.start <- min(grep("materials and methods", paper.html.temp, ignore.case = T))
}

if(length(grep("Discussion", paper.html.temp, ignore.case = T)) == 0)
{
  if(length(grep("Results", paper.html.temp, ignore.case = T)) == 0)
  {
    discussion.start <- 33
  }else{
    discussion.start <-  min(grep("Results", paper.html.temp, ignore.case = T))
  }
}else{
  discussion.start <- min(grep("Discussion", paper.html.temp, ignore.case = T)) 
}


#define range
intro.range <- Introduction.start:Method.start
paper.html.temp.intro <- paper.html.temp[intro.range]

abstract.range <- Abstract.start:Introduction.start
paper.html.temp.abstract <- paper.html.temp[abstract.range]

results.range <- Method.start:discussion.start
paper.html.temp.results <- paper.html.temp[results.range]

temp.x.abstract2 <- paste(temp.x.abstract, collapse = "") %>% strsplit("\\.")
temp.x2 <- paste(temp.x, collapse = "") %>% strsplit("\\. ")
paper.html.temp.results2 <- paste(paper.html.temp.results, collapse = "") %>% strsplit("\\. ")




#extract animal numbers
test5 = c("in our study, we used a 36 female and adult sprague dawley rats. Moreover, 2 groups were considered (n=4) and (n=7).")

if(str_detect(toString(test5), regex(animal_number_regex_total, ignore_case = T)) == T){
  if(str_detect(toString(test5), regex("\\d+", ignore_case = T)) == T){
    animal_number.temp <- str_extract_all(test5, regex(animal_number_regex_total, ignore_case = T))
    animal_number.temp <- Filter(length, animal_number.temp)
    animal_number.temp <- toString(animal_number.temp)
    animal_number.temp <- str_extract_all(animal_number.temp, regex("\\d+", ignore_case = T))
    animal_number <- sapply(animal_number.temp, function(x) sum(as.numeric(x)))
  }else{
    spelled_match.temp <- c()
    animal_number.temp <- c()
    for(numCount.numbers in seq_along(numbers_spelled)) 
    {for(numCount.animals in seq_along(animals))
    {spelled_match.temp <- str_c(numbers_spelled[[numCount.numbers]], animals[[numCount.animals]])
    print(spelled_match.temp)
    }
      animal_number.temp[[numCount.numbers]] <- str_extract_all(test5, regex(spelled_match.temp, ignore_case = T))
    }
    animal_number.temp <- Filter(length, animal_number.temp)
    animal_number.temp <- toString(animal_number.temp)
    animal_number.temp <- str_extract_all(animal_number.temp, regex(animal_number_regex_spelled, ignore_case = T))
    animal_number.temp <- unlist(animal_number.temp, use.names=FALSE)
    
    animal_number.temp2 <- c()
    for(numCount in seq_along(animal_number.temp)) 
    {
      animal_number.temp2[[numCount]] <- animal_number.temp[[numCount]] %>% word2num %>% .[[2]]
    }
    animal_number <- sum(animal_number.temp2)
  }
}else{
  if(str_detect(toString(test5), regex(animal_number_regex_n, ignore_case = T)) == T){
    animal_number.temp <- str_extract_all(test5, regex(animal_number_regex_n, ignore_case = T))
    animal_number.temp <- Filter(length, animal_number.temp)
    animal_number.temp <- toString(animal_number.temp)
    animal_number.temp <- str_extract_all(animal_number.temp, "\\d+")
    animal_number <- sapply(animal_number.temp, function(x) sum(as.numeric(x)))
    
  }else{animal_number <- NA}}
animal_number

#"total" number of animals
#str_view_all(test4, regex(animal_number_regex_total, ignore_case = T))
animal_number_regex_total <- c("\\btotal.{1,50}\\brats\\b", "\\btotal.{1,50}\\bmarmosets\\b",
                               "\\btotal.{1,50}\\bguinea pigs\\b", "\\btotal.{1,50}\\bmacaques\\b",
                               "\\btotal.{1,50}\\brhesus monkeys\\b", "\\btotal.{1,50}[^a] pigs\\b",
                               "\\btotal.{1,50}\\bswines\\b", "\\btotal.{1,50}\\bdogs\\b",
                               "\\btotal.{1,50}\\bcats\\b", "\\btotal.{1,50}\\brabbits\\b",
                               "\\btotal.{1,50}\\bprimates\\b", "\\btotal.{1,50}\\bhamsters\\b", 
                               "\\btotal.{1,50}\\bsheeps\\b", "\\btotal.{1,50}\\bmice\\b", 
                               "\\btotal.{1,50}\\banimals\\b")
animal_number_regex_total <- str_c(animal_number_regex_total, collapse = "|")
animal_number.temp <- str_extract_all(test4, regex(animal_number_regex_total, ignore_case = T))
animal_number.temp <- Filter(length, animal_number.temp)
animal_number.temp <- toString(animal_number.temp)
animal_number.temp <- str_extract_all(animal_number.temp, regex("\\d+", ignore_case = T))
animal_number <- sapply(animal_number.temp, function(x) sum(as.numeric(x)))
animal_number

test4 <- c("we used a total of 34 adult female Sprague Dawley primates for all experiments. A total of 35 male dogs; These wonderful cute hdfkghdfkgdhfkgjdfhgkfdjhgkdjfhkdjgfh and junky rats were distributed in two groups, 26 rats were in group 1 and 8 were in group 2. Also, a total of 24 young adult male macarats were used")



#n=x
str_view_all(test, regex("( |\\()n ?=( \\d|\\d)( |\\d|\\))( |\\d|\\)|,|\\.)", ignore_case = T))

animal_number_regex_n <- c("( |\\()n ?=( \\d|\\d)( |\\d|\\))( |\\d|\\)|,|\\.)")

animal_number.temp <- str_extract_all(test, regex(animal_number_regex_n, ignore_case = T))
animal_number.temp <- Filter(length, animal_number.temp)
animal_number.temp <- toString(animal_number.temp)
animal_number.temp <- str_extract_all(animal_number.temp, "\\d+")
animal_number <- sapply(animal_number.temp, function(x) sum(as.numeric(x)))
animal_number

#spelled numbers
str_view_all(test2, regex(animal_number_regex_spelled, ignore_case = T))

animal_number_regex_spelled <- c("\\bone\\b|\\btwo\\b|\\bthree\\b|\\bfour\\b|\\bfive\\b|
                                 \\bsix\\b|\\bseven\\b|\\beight\\b|\\bnine\\b|\\bten\\b")

animal_number.temp <- str_extract_all(test2, regex(animal_number_regex_spelled, ignore_case = T))
animal_number.temp <- unlist(animal_number.temp, use.names=FALSE)

number.temp <- c()
for(numCount in seq_along(animal_number.temp)) 
{
  
  number.temp[[numCount]] <- animal_number.temp[[numCount]] %>% word2num %>% .[[2]]
}
animal_number <- sum(number.temp)
animal_number

#spelled numbers with nested loop
numbers_spelled <- c("\\btwo\\b", "\\bthree\\b", "\\bfour\\b", "\\bfive\\b", "\\bsix\\b",
                     "\\bseven\\b", "\\beight\\b", "\\bnine\\b", "\\bten\\b")
animals <- c(".{1,15}\\brats\\b", ".{1,15}\\bmarmosets\\b", ".{1,15}\\bguinea pigs\\b",
             ".{1,15}\\bmacaques\\b", ".{1,15}\\brhesus monkeys\\b", ".{1,15}\\bpigs\\b",
             ".{1,15}\\bswines\\b", ".{1,15}\\bdogs\\b", ".{1,15}\\bcats\\b",
             ".{1,15}\\brabbits\\b", ".{1,15}\\bprimatess\\b", ".{1,15}\\bhamsters\\b",
             ".{1,15}\\bsheeps\\b", ".{1,15}\\bmice\\b")

spelled_match.temp <- c()
animal_number.temp <- c()
# for(numCount.numbers in seq_along(numbers_spelled)) 
# {for(numCount.animals in seq_along(animals))
# {spelled_match.temp <- str_c(numbers_spelled[[numCount.numbers]], animals[[numCount.animals]])
# print(spelled_match.temp)
# }
#   animal_number.temp[[numCount.numbers]] <- str_extract_all(test3, regex(spelled_match.temp, ignore_case = T))
# }
# animal_number.temp <- Filter(length > 0, animal_number.temp)
# animal_number.temp <- toString(animal_number.temp)
# animal_number.temp <- str_extract_all(animal_number.temp, regex(animal_number_regex_spelled, ignore_case = T))
# animal_number.temp <- unlist(animal_number.temp, use.names=FALSE)

b<-1
for(i in seq_along(numbers_spelled)) {
  for(j in seq_along(animals)) {
    spelled_match.temp <- str_c(numbers_spelled[[i]], animals[[j]])
    print(spelled_match.temp)
   animal_number.temp[[b]] <- ifelse(length(animal_number.temp) > 0, str_extract_all(test3, regex(spelled_match.temp, ignore_case = T)),NA)
  b<-b+1
  } }
animal_number.temp<-matrix(animal_number.temp,)
animal_number.temp <- na.omit(animal_number.temp)
head(animal_number.temp)


animal_number.temp2 <- c()
for(numCount in seq_along(animal_number.temp)) 
{
  animal_number.temp2[[numCount]] <- animal_number.temp[[numCount]] %>% word2num %>% .[[2]]
}
animal_number <- sum(animal_number.temp2)
animal_number



test <- c("For ex vivo NMR and flow cytometry studies of tissue samples, additional groups were prepared.
          Groups 4 and 5 were EAE (N = 4) and control (N = 5) rats, respectively, and were prepared identically to
          Groups 1 and 2, except were used exclusively for NMR analysis of fixed tissues. Based on published reports,
          we chose 6-week diet with grinded fodder with 0.2% cuprizone for group 1 mice (n=10) [4,7]. Group 2 mice 
          (n=25) received high dose of thewithdrawal, the experimental animals were  neurotoxin (0.6%) with fodder 
          for 4 weeks. After cuprizone transferred to standard fodder. Control group (n=10) received fodder without 
          cuprizone. 8-10 week old female mice were grouped into 2 groups. Control Group contained wild type (WT)
          C57BL/6 mice (B6) n = 103, and the Experimental Group contained 2D2tg mice on C57BL/6 background n = 22.")
test <- str_split(test, pattern = "\n")
test <- paste(unlist(test), sep = " ")
test <- paste(test, collapse = " ")

test2 <- c("Ten mice were used for IHC and eight mice for gene expression analysis in each group, namely placebo, the. In all the experiments, a total of nine affected animals was used")






#extract email
#str_view_all(paper.html.temp, "([_a-z0-9-]+(\\.[_a-z0-9-]+)*@[a-z0-9-]+(\\.[a-z0-9-]+)*(\\.[a-z]{2,4}))", match = T)
if(str_detect(toString(paper.html.temp), regex("([_a-z0-9-]+(\\.[_a-z0-9-]+)*@[a-z0-9-]+(\\.[a-z0-9-]+)*(\\.[a-z]{2,4}))", ignore_case = T)) == F){
  email <- NA
}else{
  email.temp <- str_extract_all(paper.html.temp, regex("([_a-z0-9-]+(\\.[_a-z0-9-]+)*@[a-z0-9-]+(\\.[a-z0-9-]+)*(\\.[a-z]{2,4}))", ignore_case = T))
  email.temp <- Filter(length, email.temp)
  email <- paste(email.temp, collapse = "; ")}
email



#year by using file file title as string
if(str_detect(pdf.test.vector[247], regex("[^\\d| ]\\d\\d\\d\\d[^\\d| ]", ignore_case = T)) == F){
  year <- NA
}else{
  year.temp <- str_extract(pdf.test.vector[247], regex("[^\\d| ]\\d\\d\\d\\d[^\\d| ]", ignore_case = T))
  year.temp <- str_extract(year.temp, regex("\\d\\d\\d\\d", ignore_case = T))
  year <- ifelse(year.temp < 2019 & year.temp > 1900, year.temp, NA)
}
year

#author by using file title as string
if(str_detect(pdf.test.vector[247], regex("/[a-z| ]+-", ignore_case = T)) == F){
  author <- NA
}else{
  author.temp <- str_extract(pdf.test.vector[247], regex("/[a-z| ]+-", ignore_case = T))
  author <- str_extract(author.temp, regex("[a-z]+", ignore_case = T))
}
author

#paper title by using file title as string
if(str_detect(pdf.test.vector[247], regex("-[a-z| ]+\\.", ignore_case = T)) == F){
  paper.title <- NA
}else{
  paper.title.temp <- str_extract(pdf.test.vector[247], regex("-[a-z| ]+\\.", ignore_case = T))
  paper.title <- str_extract(paper.title.temp, regex("[a-z| ]+", ignore_case = T))
}
paper.title

#file title
file.title.temp <- str_extract(pdf.test.vector[247], regex("/.+", ignore_case = T))
file.title <- str_extract(file.title.temp, regex("[^/].+"))


#species
species_list <- c("rat\\b", " rat ", " rat\\b", "rat ", "rats\\b", " rats ", " rats\\b", "rats ", 
                  "marmoset\\b", " marmoset ", " marmoset\\b", "marmoset ", "marmosets\\b", " marmosets ", " marmosets\\b", "marmosets ", 
                  "guinea pig\\b", " guinea pig ", " guinea pig\\b", "guinea pig ", "guinea pigs\\b", " guinea pigs ", " guinea pigs\\b", "guinea pigs ", 
                  "macaque\\b", " macaque ", " macaque\\b", "macaque ", "macaques\\b", " macaques ", " macaques\\b", "macaques ", 
                  "rhesus monkey\\b", " rhesus monkey ", " rhesus monkey\\b", "rhesus monkey ", "rhesus monkeys\\b", " rhesus monkeys ", " rhesus monkeys\\b", "rhesus monkeys ", 
                  "\\bpig\\b", "[^a] pig ", "[^a] pig\\b", "\\bpig ", "\\bpigs\\b", "[^a] pigs ", "[^a] pigs\\b", "\\bpigs ", 
                  "swine\\b", " swine ", " swine\\b", "swine ", "swines\\b", " swines ", " swines\\b", "swines ", 
                  "dog", " dog ", " dog", "dog ", "dogs", " dogs ", " dogs", "dogs ", 
                  "cat\\b", " cat ", " cat\\b", "cat ", "cats\\b", " cats ", " cats\\b", "cats ", 
                  "rabbit\\b", " rabbit ", " rabbit\\b", "rabbit ", "rabbits\\b", " rabbits ", " rabbits\\b", "rabbits ", 
                  "primate", " primate ", " primate", "primate ", "primates", " primates ", " primates", "primates ",
                  "hamster\\b", " hamster ", " hamster\\b", "hamster ", "hamsters\\b", " hamsters ", " hamsters\\b", "hamsters ", 
                  "sheep\\b", " sheep ", " sheep\\b", "sheep ", "sheeps\\b", " sheeps ", " sheeps\\b", "sheeps ", 
                  "mouse\\b", " mouse ", " mouse\\b", "mouse ", "mice\\b", " mice ", " mice\\b", "mice ")
species_match <- str_c(species_list, collapse = "|")
species <- c("rat?", "marmoset?", "guineapig?", "macaque?", "rhesusmonkey?", "pig?", "swine?",
             "dog?", "cat?", "rabbit?", "primate?", "hamster?", "sheep?", "mouse\\b", "mice\\b")
species_match2 <- str_c(species, collapse = "|")

species.temp <- str_subset(test2, species_match)
species.temp <- str_extract_all(test2, regex(species_match, ignore_case = T))
species.temp <- Filter(length, species.temp)
species.temp <- paste(unlist(species.temp), collapse = ";")
species.temp <- gsub(" ", "", species.temp)
species.temp <- str_to_lower(species.temp)
species.temp <- str_extract_all(species.temp, regex(species_match2, ignore_case = T))
species.temp <- unique(unlist(species.temp))
if(sum(str_detect(species.temp, regex("mouse", ignore_case = T))) > 0 &
   sum(str_detect(species.temp, regex("mice", ignore_case = T))) > 0){
  species.temp <- gsub("mice", "mouse", species.temp)
  species.temp <- unique(unlist(species.temp))
}
if(sum(str_detect(species.temp, regex("mice", ignore_case = T))) > 0){
  species.temp <- gsub("mice", "mouse", species.temp)
}
species <- paste(unlist(species.temp), collapse = ";")
species

#alternative approach: order species according to most common appearance, seems better
species_list <- c("\\brat[^ihe]", "\\bmarmoset", "\\bguinea pig", "\\bmacaque", "\\brhesus monkey",
                  "\\bswine", "\\bdog", "\\bcat", "\\brabbit", "\\bprimate",
                  "\\bhamster", "\\bsheep", "\\bmouse", "\\bmice")

species.count.temp <- vector()
species.temp <- vector()
for(numCount in seq_along(species_list)) 
{
  print(str_count(test, regex(numCount)))
  species.count.temp[[numCount]] <- sum(str_count(test, regex(species_list[[numCount]])))
  species.temp[[numCount]] <- str_extract(test, regex(species_list[[numCount]]))
}
species.temp.df <- as.data.frame(t(rbind(species.temp, species.count.temp)), stringsAsFactors = F)
species.temp.df$species.count.temp <- as.numeric(species.count.temp)
species.temp.df$species.temp <- as.character(species.temp)
species.temp.df <- species.temp.df[order(desc(species.temp.df$species.count.temp)),]
species.temp.df <- species.temp.df[!(species.temp.df$species.count.temp == 0),]

species1 <- paste(species.temp.df$species.temp[1])
species_count1 <- paste(species.temp.df$species.count.temp[1])
species2 <- paste(species.temp.df$species.temp[2])
species_count2 <- paste(species.temp.df$species.count.temp[2])

if(species1 == "NA"){species1 <- NA}
if(species2 == "NA"){species2 <- NA}

#str_view_all(test, regex(species_list, ignore_case = T))

##Sex
female_regex <- "\\bfemale\\b"
male_regex <- "\\bmale\\b"

has_female <- str_subset(test3, regex(female_regex, ignore_case = T))
has_female <- str_extract(test3, regex(female_regex, ignore_case = T))
has_female <- paste(unlist(has_female), collapse = ";")
ifelse(str_detect(has_female, regex(female_regex, ignore_case = T)) == F, female <- 0, female <- "female")
female

has_male <- str_subset(test3, regex(male_regex, ignore_case = T))
has_male <- str_extract(test3, regex(male_regex, ignore_case = T))
has_male <- paste(unlist(has_male), collapse = ";")
ifelse(str_detect(has_male, regex(male_regex, ignore_case = T)) == F, male <- 0, male <- "male")
male

if(sum(str_detect(female, regex("female", ignore_case = T))) > 0 &
   sum(str_detect(male, regex("male", ignore_case = T))) == 0){
  sex.test <- "female"
}
if(sum(str_detect(female, regex("female", ignore_case = T))) == 0 &
   sum(str_detect(male, regex("male", ignore_case = T))) > 0){
  sex.test <- "male"
}
if(sum(str_detect(female, regex("female", ignore_case = T))) > 0 &
   sum(str_detect(male, regex("male", ignore_case = T))) > 0){
  sex.test <- "female and male"
}
if(sum(str_detect(female, regex("female", ignore_case = T))) == 0 &
   sum(str_detect(male, regex("male", ignore_case = T))) == 0){
  sex.test <- NA
}
sex.test

test <- c("she is Female but only during nighttime. Male is she during daytime and male she stays female, maleficient group", "female No1", "male")
test2 <- c("she is female")
test3 <- c("male he is male")

str_view_all(test2, regex("\\bmale\\b", ignore_case = T))
str_view_all(test, regex("\\bfemale\\b", ignore_case = T))

###Disease
#define individual disease regex
multiplesclerosis_regex <- c("\\bMS\\b|\\b[M|m]ultiple [S|s]clerosis\\b")
amytrophiclateralsclerosis_regex <- c("\\bALS\\b|\\bamyotrophic lateral sclerosis\\b")
alzheimer_regex <- c("\\bAD\\b|\\b[A|a]lzheimer's [D|d]isease\\b|\\b[A|a]lzheimer's [D|d]ementia\\b|\\b[A|a]lzheimer [D|d]isease\\b")

#create a collection of all disease regex
disease_regex <- c(multiplesclerosis_regex, alzheimer_regex, amytrophiclateralsclerosis_regex)

#mine for the disease and the number of appearances and show the first two most common hits
disease.count.temp <- vector()
disease.temp <- vector()
for(numCount in seq_along(disease_regex)) 
{
  print(str_count(test, regex(numCount)))
  disease.count.temp[[numCount]] <- sum(str_count(test, regex(disease_regex[[numCount]], ignore_case = F)))
  disease.temp[[numCount]] <- str_extract(test, regex(disease_regex[[numCount]], ignore_case = F))
}
disease.temp.df <- as.data.frame(t(rbind(disease.temp, disease.count.temp)), stringsAsFactors = F)
disease.temp.df$disease.count.temp <- as.numeric(disease.count.temp)
disease.temp.df$disease.temp <- as.character(disease.temp)
disease.temp.df <- disease.temp.df[order(desc(disease.temp.df$disease.count.temp)),]

disease1 <- paste(disease.temp.df$disease.temp[1])
disease_count1 <- paste(disease.temp.df$disease.count.temp[1])
disease2 <- paste(disease.temp.df$disease.temp[2])
disease_count2 <- paste(disease.temp.df$disease.count.temp[2])

if(disease1 == "NA"){disease1 <- NA}
if(disease2 == "NA"){disease2 <- NA}

test <- c("The mouse annualized relapse rate was lower with Multiple sclerois (MS) ocrelizumab than with interferon
          beta-1a in trial 1 (0.16 vs. 0.29; 46% lower rate with ocrelizumab; P<0.001) and in trial 2 (0.16 vs. 0.29; 47% lower
          rate; P<0.001). In rats prespecified pooled analyses, the percentage of patients with disability
          progression confirmed at 12 weeks was significantly lower with ocrelizumab than with interferon
          beta-1a (9.1% vs. 13.6%; MS hazard alzheimer's disease ratio, 0.60; 95% confidence interval [CI], 0.45 to 0.81; P<0.001),
          as was the percentage of patients rats with disability progression rat confirmed at 24 weeks multiple
          sclerosis (6.9% vs. 10.5%; ALS hazard ratio, 0.60; 95% CI, 0.43 to 0.84; P=0.003). The mean number of
          gadolinium-enhancing lesions per T1-weighted multiple sclerosis magnetic resonance scan was 0.02 with ocrelizumab
          versus 0.29 with interferon beta-1a in trial 1 (94% AD lower number of lesions with ocrelizumab,
          P<0.001) and 0.02 versus 0.42 in trial 2 (95% lower number of lesions, P<0.001). MS The change in
          the Multiple Sclerosis amyotrophic lateral sclerosis (ALS) rat Functional rat Composite score (a composite measure of walking speed, upper-limb
          movements, and cognition; for this z score, mouse negative values indicate worsening and positive values
          indicate improvement) significantly favored ocrelizumab over interferon beta-1a in trial 2
          (0.28 vs. 0.17, P=0.004) but MS not in trial 1 (0.21 vs. 0.17, P=0.33). Infusion-related reactions
          occurred in 34.3% of the patients treated with ocrelizumab. ALS Serious Multiple sclerosis infection occurred in 1.3% of
          the patients treated with ocrelizumab and in 2.9% of those treated with interferon beta-1a.
          Neoplasms occurred in 0.5% of the patients amyotrophic lateral sclerosis rat treated with ocrelizumab and in 0.2% of those treated
          with interferon beta-1a.", "MS is a disease multiple sclerosis")
test <- str_split(test, pattern = "\n")
test <- paste(unlist(test), sep = " ")
test <- paste(test, collapse = " ")


sum(str_count(test, regex(multiplesclerosis_regex)))
sum(str_count(test, regex(amytrophiclateralsclerosis_regex)))
sum(str_count(test, regex(alzheimer_regex)))


#Disease model
EAE_regex <- c("\\bEAE\\b|\\b[E|e]xperimental [A|a]utoimmune [E|e]ncephalomyelitis\\b")
lysolecithin_regex <- c("\\bLPC\\b|\\b[L|l]ysolecithin\\b")
ethidiumbromide_regex <- c("\\b[E|e]thidium bromide\\b")
cuprizone_regex <- c("\\b[C|c]uprizone\\b")

#create a collection of all model regex
model_regex <- c(EAE_regex, lysolecithin_regex, ethidiumbromide_regex, cuprizone_regex)

#mine for the model and the number of appearances and show the first two most common hits
model.count.temp <- vector()
model.temp <- vector()
for(numCount in seq_along(model_regex)) 
{
  print(str_count(test, regex(numCount)))
  model.count.temp[[numCount]] <- sum(str_count(test, regex(model_regex[[numCount]], ignore_case = F)))#ignore case!
  model.temp[[numCount]] <- str_extract(test, regex(model_regex[[numCount]], ignore_case = F))#ignore case!
}
model.temp.df <- as.data.frame(t(rbind(model.temp, model.count.temp)), stringsAsFactors = F)
model.temp.df$model.count.temp <- as.numeric(model.count.temp)
model.temp.df$model.temp <- as.character(model.temp)
model.temp.df <- model.temp.df[order(desc(model.temp.df$model.count.temp)),]

model1 <- paste(model.temp.df$model.temp[1])
model_count1 <- paste(model.temp.df$model.count.temp[1])
model2 <- paste(model.temp.df$model.temp[2])
model_count2 <- paste(model.temp.df$model.count.temp[2])

if(model1 == "NA"){model1 <- NA}
if(model2 == "NA"){model2 <- NA}



test <- c("EAE and experimental autoimmune encephalomyelitis was our model and Lysolecithin")
str_view_all(test, regex(model_regex, ignore_case = F))

#Risk of bias
#Welfare
if(str_detect(test2, regex("approved|approval|carried out in accordance|conducted in accordance", ignore_case = T)) == T){
  welfare <- "yes"
}else{
  welfare <- "not reported"
}

#Blinding
if(str_detect(test2, regex("blinded|were blind", ignore_case = T)) == T){
  blinding <- "yes"
}else{
  blinding <- "not reported"
}

#Randomization
if(str_detect(test2, regex("random", ignore_case = T)) == T){
  randomization <- "yes"
}else{
  randomization <- "not reported"
}

#sample size calculation
if(str_detect(test2, regex("achieve.{1,10}power|sample size calculation|power analysis", ignore_case = F)) == T){
  powertest <- "yes"
}else{
  powertest <- "not reported"
}

#ARRIVE guidelines
if(str_detect(test2, regex("ARRIVE", ignore_case = F)) == T){#ignore case!
  ARRIVE <- "yes"
}else{
  ARRIVE <- "not reported"
}

#conflict of interest
if(str_detect(test2, regex("conflict[s| ] ?of interest|disclosures|competing interests|nothing to disclose", ignore_case = T)) == T){
  disclosure <- "yes"
}else{
  disclosure <- "not reported"
}

str_view_all(test, regex("conflict[s| ] ?of interest|disclosures|competing interests|nothing to disclose", ignore_case = T))

test <- c("all animal exp were approved by the conflicts of interest, disclose, nothing to disclose, local ethicscommittee and were according to conflict of interest ARRIVE and arrive in accordance")
test2 <- c(" my bonnie is over the ocean approval")
str_view_all(test, regex("approved|carried out in accordance|conducted in accordance", ignore_case = T))
str_view_all(test, regex("ARRIVE", ignore_case = F))









test <- c("my rat and her", "guinea pig", "micery", "These Rats ", "several marmosets were", "the marmosets", "rat", " one mouse", "21 rats were", "rats, mice and pigs", "marmoset", "Marmosets", "ratumbi", "our ratsource", "mihgtice")
test2 <- c(" I had several experiments with rats and guinea pigs but only the experiments with mice worked. Marmoset. Hamster, rats and a monkey was used")
str_view_all(test2, regex(species_match, ignore_case = T))

has_species <- gsub(has_species, "[^ | $]", "a")

str(has_species)
species_match <- str_extract_all(has_species, species_match, simplify = T)

str_view_all(test, regex(species_match, ignore_case = T))
str_view_all(test, regex(" ?rat[s\\s|\\s]\\b", ignore_case = T))
str_view_all(test, regex("rat[s|\\s]\\b", ignore_case = T))
str_view_all(test, regex(" ?marmoset[s\\s|\\s]\\b", ignore_case = T))
str_view_all(test, regex(" ?guinea pig[s\\s|\\s]\\b", ignore_case = T))
str_view_all(test, regex(" ?macaque[s\\s|\\s]\\b", ignore_case = T))
str_view_all(test, regex(" ?rhesus monkey[s\\s|\\s]\\b", ignore_case = T))
str_view_all(test, regex(" ?pig[s\\s|\\s]\\b", ignore_case = T))
str_view_all(test, regex(" ?swine[s\\s|\\s]\\b", ignore_case = T))
str_view_all(test, regex(" ?dog[s\\s|\\s]\\b", ignore_case = T))
str_view_all(test, regex(" ?cat[s\\s|\\s]\\b", ignore_case = T))
str_view_all(test, regex(" ?rabbit[s\\s|\\s]\\b", ignore_case = T))
str_view_all(test, regex(" ?primate[s\\s|\\s]\\b", ignore_case = T))
str_view_all(test, regex(" ?hamster[s\\s|\\s]\\b", ignore_case = T))
str_view_all(test, regex(" ?sheep[s\\s|\\s]\\b", ignore_case = T))
str_view_all(test, regex(" ?m[(ous)|(ic)]e", ignore_case = T))#does not work










str_view_all(test, regex(" rat[s?|\b]", ignore_case = T))



test <- c("Rej@jgh.ch", "hsdkfjsdkg.ggjhg@huhw.com", "@@zut", "riibbii")
test2 <- c("@@@", "hjskjdfhkjds@", "hskfjh", "fjhall", "2")
test3 <- c(" ich war ein kleines ribisohn", " dies ist eine amila sfdfkl@jsfdlkfj.com", "hier keine scheisse")
str(email.temp)

str_view_all(test2, "([_a-z0-9-]+(\\.[_a-z0-9-]+)*@[a-z0-9-]+(\\.[a-z0-9-]+)*(\\.[a-z]{2,4}))")

library(RCurl)
library(XML)
webpage <- getURL("https://doi.org/10.1016/j.jneuroim.2010.08.022")
webpage <- readLines(tc <- textConnection(webpage)); close(tc)
pagetree <- htmlTreeParse(webpage, error=function(...){}, useInternalNodes = TRUE)
# parse the tree by tables
x <- xpathSApply(pagetree, "//*/table", xmlValue)  
# do some clean up with regular expressions
x <- unlist(strsplit(x, "\n"))
x <- gsub("\t","",x)
x <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", x, perl=TRUE)
x <- x[!(x %in% c("", "|"))]

number <- "( |\\() (N =)"
str_view_all(test, "number")

test <- str_sub(doi, 5, max(doi))
test
str_view(paper.html.temp.results, "Fig \\d.", match = T)
str_view(c("grey", "gray"), "gr(e|a)y")
str_view(test, "\\{.+\\}", match = T)

str_view(test, "\\d\\d\\d \\d\\d\\d \\d\\d \\d\\d", match = T)
str_view(test, "\\d{3} \\d{3} \\d{2} \\d{2}", match = T)
str_view(test, "^[^aeiou]{3}", match = T)
str_view(test, "[aeiou]{3,}", match = T)
str_view(fruit, "(.)\\1", match = T)

str_detect(fruit, "an")

str_view(fruit, "^ba", match = T)
sum(str_detect(fruit, "^ba"))


test <- str_c(paper.html.temp, collapse = "")
test_df <- tibble(line = 1:length(test), text = test)
test_df %>%  unnest_tokens(word, test)
test <- c("de2", "at3", "ribi", "076 391 04 01", "{roemer}", "{asabara}", "{asdd")

#str_which = grep
length(str_which(fruit, "berry"))#14 -> number of appearances
length(grep("berry", fruit))#same thing as str_which
#if you want to ignore cases:
length(str_which(fruit, regex("berry", ignore_case = TRUE)))

sum(str_which(fruit, "berry"))#475
str_view(fruit, "berry", match = T)


#str_detect = grepl
str_detect(fruit, regex("berry", ignore_case = T))
grepl("berry", fruit, ignore.case = T)

fruit[str_detect(fruit, "berry$")]
str_subset(fruit, "berry$")#same thing

df <- tibble(word = words, i = seq_along(word))
df %>% filter(str_detect(word, "x$"))

sum(str_detect(words, "x$"))
sum(str_detect(words, "^x"))
sum(str_detect(words, "^x|x$"))
str_view_all(words, "x$", match = T)
str_view(words, "^[eiou].*[aeiou]$", match = T)
max(str_detect(words, "[eaiou]"))
max(str_subset(words, "[eaiou]"))
max(str_count(words, "[eaiou]"))
max(str_count(words, "[^eaiou]"))



str_count(fruit, "a")
str_view(fruit, "a?")
fruit_df <- tibble(word = fruit, i = seq_along(word))

fruit_df %>% mutate(berry = str_detect(word, "berry"), blue = str_detect(word, "blue"))

head(sentences)
colours <- c(" red ", " green ", " purple ", " orange ", " yellow ", "brown")
colour_match <- str_c(colours, collapse = "|")
has_colour <- str_subset(sentences, colour_match)
matches <- str_extract_all(has_colour, colour_match, simplify = T)

str_extract_all(sentences, "^[A-Z][a-z]*\\s", simplify = T)#first sentence of each word
str_extract_all(sentences, ".*ing$", simplify = T)


noun <- "(a|the) ([^ ]+)"
has_noun <- sentences %>% str_subset(noun) %>% head(10)
has_noun %>% str_extract(noun)
has_noun %>% str_match(noun)

sequence1 <- "( one| two| three| four| five| six| seven| eight| nine) ([^ ]+)"
has_sequence1 <- sentences %>% str_subset(sequence1) %>% head(10)
has_sequence1 %>% str_extract(sequence1)
has_sequence1 %>% str_match(sequence1)

has_sequence1 %>% str_replace_all("seven", "7")

sentences %>% str_split(" ", n = 5, simplify = T) %>%  head(3)
sentences %>% str_split(boundary("word"), simplify = T) %>% head(3)#same thing

test <- c("this banana is an apple and is not a banana")
str_split(test, boundary("word"), simplify = T)
str_split(test, "", simplify = T)

sentences %>%  str_locate("The ") %>% head(5)
sentences %>%  head(5)

test <- "Line 1\nLine2 \nLine3"
test %>% str_extract_all("Line", multiline = T)
str_extract_all(test, "Line", multiline = T)


f













#str_which = grep
#str_detect = grepl
#stringi much faster than R base

#grep provides the position of a pattern in the character vector, just as it's equivalent str_which does:
#grep("Lorem", sample_small)
#> [1]  1  9 14 32 45 50 65 93 94
#str_which(sample_small, "Lorem")
#> [1]  1  9 14 32 45 50 65 93 94

#grepl/str_detect on the other hand give you the information for each element of the vector, if it contains the string or not.
#grepl("Lorem", sample_small)
#>   [1]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE
#>  [12] FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [23] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
#>  [34] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [45]  TRUE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE
#>  [56] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
#>  [67] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [78] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [89] FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE
#> [100] FALSE
#str_detect(sample_small, "Lorem")
#>   [1]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE
#>  [12] FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [23] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
#>  [34] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [45]  TRUE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE
#>  [56] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
#>  [67] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [78] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [89] FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE
#> [100] FALSE














f









if(Introduction.start < Method.start){intro.range <- Introduction.start:Method.start
}else{intro.range <- Introduction.start+80}#which value covers 90% of all Intros? 80 is a raw estimate

temp.x.intro <- temp.x[intro.range]
temp.x.intro <- paste(unlist(temp.x.intro), sep = " ")
temp.x.intro <- paste(temp.x.intro, collapse = " ")



test <- read.pdf.EE(pdf.rob.vector[16])
str_count(test, regex("competing interests|open", ignore_case = T))
if(str_detect(test, regex("approved|approval|carried out in accordance|conducted in accordance|in accordance", ignore_case = T)) == T){
  welfare <- "yes"
}else{
  welfare <- "not reported"
}

test <- c("I received ya mother and a fudi")
if(str_detect(test, regex("conflict[s| ] ?of interest|disclosures|competing interests|nothing to disclose|received grant[s| ]|disclosure statement", ignore_case = T)) == T){
  disclosure <- "yes"
}else{
  disclosure <- "not reported"
}




Conflict_list <- c("\\bachieve.{1,10}power", "\\bsample size calculation", "\\bpower analysis", "\\bstudy power|sample size")


a <- data.frame(A = c("Not reported", "Yes", "Not reported", "Not reported"),
                B = c("Yes", "Yes", "Not reported", "Not reported"),
                C = c("Not reported", "Yes", "Not reported", "Yes"))

b <- data.frame(A = c("Yes", "Yes", "Yes", "Not reported"),
                B = c("Yes", "Not reported", "Not reported", "Not reported"),
                C = c("Not reported", "Not reported", "Not reported", "Yes"))
a == b
c <- table(a == b)
c[names(c)==F]
as.numeric(c[names(c)==T])/(nrow(a)*ncol(a))*100



if(str_detect(a, regex("approved|approval|carried out in accordance|conducted in accordance", ignore_case = T)) == T){
  welfare <- "yes"
}else{
  welfare <- "not reported"
}
welfare

a <- c("I love ribi approved and a approval ads")

welfare.count.temp <- vector()
welfare.temp <- vector()
for(numCount in seq_along(welfare_list)) 
{
  welfare.count.temp[[numCount]] <- sum(str_count(a, regex(welfare_list[[numCount]])))
  welfare.temp[[numCount]] <- str_detect(a, regex(welfare_list[[numCount]]))
}
welfare.temp.df <- as.data.frame(t(rbind(welfare.temp, welfare.count.temp)), stringsAsFactors = F)
welfare.temp.df$welfare.count.temp <- as.numeric(welfare.count.temp)
welfare.temp.df$welfare.temp <- as.character(welfare.temp)
welfare.temp.df <- welfare.temp.df[order(desc(welfare.temp.df$welfare.count.temp)),]
welfare.temp.df <- welfare.temp.df[!(welfare.temp.df$welfare.count.temp == 0),]

welfare1 <- paste(welfare.temp.df$welfare.temp[1])
welfare_count1 <- paste(welfare.temp.df$welfare.count.temp[1])

if(welfare1 == "NA"){welfare <- NA}

welfare1
welfare_count1

if(sum(str_count(a, regex("approved|approval|carried out in accordance|conducted in accordance"))) > 1){
  welfare_signal <- "WTF!!"
}else{
  welfare_signal <- "You're good"
}
welfare_signal


welfare_list <- c("\\bapproved", "\\bapproval", "\\bcarried out in accordance", "\\bconducted in accordance")

