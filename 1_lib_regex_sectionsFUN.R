library(pdftools)
library(dplyr)
library(tidyverse)
library(readxl)
library(rvest)
library(stringr)
library(htmltools)
library(htmlwidgets)
library(writexl)
library(formattable)
library(caret)
library(data.table)
library(rvest)
library(reticulate)

######
#REGEX
######

protocol_regex <-  paste0(c("(?i)prospero","(?<!not )(pre(|-)registered)", "prospective register of systematic reviews","prospero", "open science framework", "(\\b|\\()osf(\\b|\\()", "crd(|-| )number", "crd\\d{8,14}\\b"), collapse = "|")
#excludes 'not pre(|-)registered'


####Study goal (only in introduction and abstract)
goal_regex <- paste0(c("(?i)(we|study) aimed at","study aimed to", "study goal", "study aim", "(we|this study) set out to","This review aimed to","We undertook this systematic review in order to","purpose of this study", "study objective(s|)", "study purpose","This systematic review aimed to","This systematic review aims to",
                       "purpose of this systematic review","aim of the study","aim of this systematic review","aim of the systematic review","objective: ","aim: ","aims: ",
                       "we sought to","we aim to","we aimed to","The aim of the present",""), collapse = "|")


####In- and exclusion criteria (method section)
inclusion_regex <- paste0(c("(?i)inclusion[ |-]criteria", "exclusion[ |-]criteria", "eligibility[ |-]criteria",
                            "eligibility[ |-]criteria", "included the following", "excluded the following",
                            "were excluded", "were included", "in(-| )and exclusion",
                            "ex(-| )and inclusion", "inclusion and exclusion"), collapse = "|")


####Database search (method section)
database_regex <- paste0(c("(?i)Medline", "Pubmed", "Embase", "Cinahl", "EBSCO", "Scopus", "CENTRAL",
                           "Web of Science", "Google Scholar", "Science Direct", "ScienceDirect",
                           "OpenGrey", "Cochrane", "Web of Knowledge",
                           "China National Knowledge Infrastructure", "CNKI", 
                           "Oriental Medicine Advanced Searching Integrated System", "OASIS",
                           "China Biology Medicine disc", "CBMdisc", "Wan Fang",
                           "Chinese Clinical Trial Registry", "ClinicalTrials.gov", "BIOSIS",
                           "SciELO", "SportDiscus", "EBSCO", "PsycINFO", "Ovid",
                           "international scientific index", "Central Register of Controlled Clinical Trials",
                           "ICTRP", "IEEE", "Internet Stroke Center", "PEDro", "Biological Abstracts",
                           "Toxicology Data Network", "China Science and Technology Journal", "NLM Gateway",
                           "SpringerLink", "CAB abstracts", "proquest natural science", "Open alex"), collapse = "|")

#CAVE: at least 2 of these should be reported


####Search date (method section)
date_regex <- paste0(c("(\\b\\d{1,2}\\d{0,3})?\\b(?:jan(?:uary)?|feb(?:ruary)?|mar(?:ch)?|apr(?:il)?|may|jun(?:e)?|jul(?:y)?|aug(?:ust)?|sep(?:tember)?|oct(?:ober)?|(nov|dec)(?:ember)?)\\d?(\\d{1,2}\\d?)?\\d?((19[7-9]\\d|20\\d{2})|\\d{2})","\\b(january|february|march|april|may|june|july|august|september|october|november|december)\\s+\\d{1,2},\\s+\\d{4}\\b","search date", "publication date", "search was executed on", "from inception","(?:jan(?:uary)?|feb(?:ruary)?|mar(?:ch)?|apr(?:il)?|may|jun(?:e)?|jul(?:y)?|aug(?:ust)?|sep(?:tember)?|oct(?:ober)?|nov(?:ember)?|dec(?:ember)?)\\s+\\d{4}",
                       "(\\d{1,2}\\s+(?:jan(?:uary)?|feb(?:ruary)?|mar(?:ch)?|apr(?:il)?|may|jun(?:e)?|jul(?:y)?|aug(?:ust)?|sep(?:tember)?|oct(?:ober)?|nov(?:ember)?|dec(?:ember)?)\\s+\\d{4})","\b(?:jan(?:uary)?|feb(?:ruary)?|mar(?:ch)?|apr(?:il)?|may|jun(?:e)?|jul(?:y)?|aug(?:ust)?|sep(?:tember)?|oct(?:ober)?|(?:nov|dec)(?:ember)?)\\s+\\d{4}\b"), collapse = "|")


####Search string (method section)
string_regex <- paste0(c("(?i)search string", "boolean operator", "search was conducted using",
                         "using the following key(| |-)words", "search strategy", 
                         "using the following MeSH", "relevant key(| |-)words", "relevant Mesh(-| )terms",
                         "key(| |-)word", "emtree", "key search term", "search syntax",
                         "pre(|-)defined search", "MeSH(|-)term"), collapse = "|")


####Two reviewer (method section)
induplicate_regex_one <- paste0(c("(?i)(two|2) author", "(two|2) reviewer", "(two|2) screener", "in duplicate", 
                                  "assessed independently", "screened independently", "(two|2) researcher",
                                  "resolved by consensus", "resolved by discussion", "dis(-|)agreement",
                                  "discussed among us", "random check", "(two|2) investigator",
                                  "(two|2) independent author", "(two|2) independent investigator",
                                  "(two|2) independent researcher", "paired reviewer", "second reviewer",
                                  "second authors", "second screener", "under the guidance", "double(|-)checked"), collapse = "|")

induplicate_regex_two <- "([A-Z](\\.?)){2,3}(, | )(and |)([A-Z](\\.?)){2,3}"
#CAVE: Ignore case for induplicate_regex_two = FALSE!


####Risk of bias assessment (entire paper except references)
rob_regex <- paste0(c("(?i)risk of bias", "Coleman methodology", "quality assessment", "appraise",
                      "Newcastle(-| )ottawa", "methodological quality", "bias assessment", "ROBINS",
                      "AMSTAR", "Joanna(-| )briggs", "JBI"), collapse = "|")


####In accordance with guidelines (entire paper except references)
guidelines_regex <- paste0(c("(?i)PRISMA", "Preferred Reporting Items for Systematic Reviews and Meta-Analyses",
                             "CAMARADES", "CONSORT", "QUORUM"), collapse = "|")
#CAVE: Ignore case = FALSE (except for "Preferred Reporting Items for Systematic Reviews and Meta-Analyses")

####Meta-analysis (method section)
metaanalysis_regex <- paste0(c("(?i)random(-| )effect", "fixed(-| )effect", "meta-regression",
                               "DerSimonian", "Laird", "heterogeneity", "Higgins", "Hedges",
                               "Thompson", "subgroup analysis", "publication bias", "SMD",
                               "standardized mean", "p(-| )value"), collapse = "|")

####Flowchart (entire paper except references)
flowchart_regex <- paste0(c("(?i)flow( |-|)chart", "flow( |-|)diagram", "study selection procedure",
                            "PRISMA( |-|)chart", "PRISMA( |-|)diagram", "\\(PRISMA\\) diagram","ow-chart","ow chart" ,"ow-diagram","ow diagram"), collapse = "|")


####Conflict of interest (entire paper including after the reference section)
coi_regex <- paste0(c("(?i)nothing to disclose", "competing interest","declare(d)?(\\s+no(\\s+potential)?)?\\s+conflict(s)?","conflicting interests","no con&#xfb02;ict(s)?",
                      "conflict(-| )of(-| )interest", "financial conflict", "declaration(s|) of interest","no conflicts of interest to declare","no conflicts to declare",
                      "financial support","absence of .* potential con&#xfb02;ict of interest","absence of .* potential conflict of interest"), collapse = "|")


####Protocol deviations (method section)
deviation_regex <- paste0(c("(?i)protocol deviation", "protocol alteration", "deviate from the protocol",
                            "protocol adjustment", "protocol change", "deviations from"), collapse = "|")

####################
#SECTIONING FUNCTION
####################

temp_mapping_function<-function(html_file){
  
  #text to lower case for easier mining, removing "-" as it gives problems in the regex text matching
  html_text <- tolower(readLines(html_file))
  html_text<-gsub("-","",html_text)
  
  #define the pattern to look for in the titles in introduction or background
  pattern_intro <- "<span(.+?)introduction(.+?)span>" 
  pattern_background <- "<span(.+?)background(.+?)span>" 
  
  #extract the pattern
  matched_text_intro <- unique(str_extract_all(html_text, pattern_intro, simplify = TRUE)) 
  
  #extract the pattern in case introduction is not found, look for "background"
  if(dim(matched_text_intro)[1]==0){
    matched_text_intro <- unique(str_extract_all(html_text, pattern_background, simplify = TRUE)) 
    
  } else {
    matched_text_intro<-matched_text_intro
  }
  
  #removing empty matches
  matched_text_intro<-unique(matched_text_intro[matched_text_intro != ""])
  
  #in case multiple matches were found in looking for introduction, selecting the one with the biggest font size
  if (length(matched_text_intro)==1) {
    matched_text_intro<-matched_text_intro
  } else {
    intro_font_size<-paste0(digits(max(data.frame(matched_text_intro) %>%
                                         mutate(font_size=regmatches(matched_text_intro, gregexpr("[0-9]+.[0-9]+pt", matched_text_intro))) %>% 
                                         unnest(font_size) %>% 
                                         mutate(font_size=as.numeric(gsub("pt","",font_size))) %>% 
                                         select(font_size) %>% pull()),1),"pt")
    
    matched_text_intro<-data.frame(matched_text_intro) %>% 
      filter(grepl(intro_font_size,matched_text_intro)) %>% pull()
  }
  
  if (length(matched_text_intro)==1) {
    matched_text_intro<-matched_text_intro
  } else{
    matched_text_intro<-matched_text_intro[1]
  }
  
  #generalize it to all the titles
  general_pattern<-gsub("(>[^<]*)\\b\\d+\\b([^<]*<)", "\\1[0-9]+\\2", matched_text_intro) 
  general_pattern<-gsub("introduction|background", "(.+?)", general_pattern)
  
  #look for all the titles
  matched_texts<-unique(str_extract_all(html_text, general_pattern, simplify = TRUE)) 
  
  #in case of multiple matches, select those that are actual titles by looking at the lines in the html text with the lowest number of characters. This is repeated for all the sections of the paper
  section_intro<-matched_texts[grepl("introduction|background", matched_texts)]
  if (length(section_intro)==1) {
    section_intro<- section_intro
  } else {
    section_intro<-section_intro[which.min(nchar(section_intro))]
  }
  
  section_method<-matched_texts[grepl("method|search strategy", matched_texts)]
  if (length(section_method)==1) {
    section_method<- section_method
  } else {
    section_method<-section_method[which.min(nchar(section_method))]
  }
  
  section_results<-matched_texts[grepl("results", matched_texts)]
  if (length(section_results)==1) {
    section_results<- section_results
  } else {
    section_results<-section_results[which.min(nchar(section_results))]
  }
  
  section_discussion<-matched_texts[grepl("discussion", matched_texts)]
  if (length(section_discussion)==1) {
    section_discussion<- section_discussion
  } else {
    section_discussion<-section_discussion[which.min(nchar(section_discussion))]
  }
  
  sections<-c(section_intro,section_method,section_results,section_discussion)
  sections<-gsub("-","",sections)
  
  #creating a dataframe containg the matches found, and the lines at which they were found 
  sections_df<-
    data.frame(html_text) %>% 
    mutate(match=str_extract_all(html_text, paste0(sections, collapse = "|"), simplify = TRUE)) %>% 
    mutate(start=1:n()) %>% 
    filter(!(match=="")) %>%
    mutate(match=str_extract_all(match, "background|introduction|methodology|methods|results|discussion|references|search strategy", simplify = TRUE)) %>% 
    mutate(match=gsub("background","introduction",match),
           match=gsub("methodology","methods",match)) %>% 
    select(!(html_text)) %>% 
    mutate(end = lead(start)-1) %>% 
    replace_na(list(end=length(html_text))) %>% 
    rbind(data.frame(match="paper",start=0,end=length(html_text)))
  
  #creating the different sections of the paper
  temp.x.intro <- paste(gsub("<.*?>", "", html_text[sections_df$start[sections_df$match=="introduction"]:sections_df$end[sections_df$match=="introduction"]]), collapse = " ")
  temp.x.intro<-gsub("<img(.+?)>","",temp.x.intro)
  
  temp.x.abstract <- paste(gsub("<.*?>", "", html_text[1:sections_df$start[sections_df$match=="introduction"]]), collapse = " ")
  
  temp.x.abstract<-gsub("<img(.+?)>","",temp.x.abstract)
  
  temp.x.method <- paste(gsub("<[^>]+>", "", html_text[sections_df$start[sections_df$match=="methods"]:sections_df$end[sections_df$match=="methods"]]), collapse = " ")
  
  temp.x.method<-gsub("<img(.+?)>","",temp.x.method)  
  
  temp.x.result <- paste(gsub("<.*?>", "", html_text[sections_df$start[sections_df$match=="results"]:sections_df$end[sections_df$match=="results"]]), collapse = " ")
  
  temp.x.result<-gsub("<img(.+?)>","",temp.x.result) 
  
  temp.x.discussion <- paste(gsub("<.*?>", "", html_text[sections_df$start[sections_df$match=="discussion"]:sections_df$end[sections_df$match=="discussion"]]), collapse = " ")
  
  temp.x.discussion<-gsub("<img(.+?)>","",temp.x.discussion) 
  
  temp.x.paper <- paste(gsub("<.*?>", "", html_text[sections_df$start[sections_df$match=="paper"]:sections_df$end[sections_df$match=="paper"]]), collapse = " ")
  
  temp.x.paper<-gsub("<img(.+?)>","",temp.x.paper) 
  
  return(list(abstract=temp.x.abstract,
              introduction=temp.x.intro,
              methods=temp.x.method,
              results=temp.x.result,
              discussion=temp.x.discussion,
              paper=temp.x.paper))
  
}