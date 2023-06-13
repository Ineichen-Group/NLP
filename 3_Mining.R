#Mining dataframe

papers_rob_df<-data.frame()

files<- list.files("html",full.names = TRUE)

for(i in files){
  list[i]<-try({
    paper_mapped<-temp_mapping_function(i)
    
    df<-data.frame(cbind(
      file_name=i,
      
      Rob_protocol=any(str_detect(c(paper_mapped$methods,paper_mapped$introduction,paper_mapped$results,paper_mapped$discussion,paper_mapped$abstract), protocol_regex)),
      
      Rob_question=any(str_detect(c(paper_mapped$introduction,paper_mapped$abstract), goal_regex)),
      
      Rob_inexclusioncriteria=str_detect(c(paper_mapped$methods), goal_regex),
      
      Rob_search=str_count(paper_mapped$methods, database_regex) >= 2,
      
      Rob_date=any(str_detect(c(paper_mapped$methods,paper_mapped$abstract), date_regex)),
      
      Rob_string=str_detect(c(paper_mapped$methods), string_regex),
      
      Rob_rob=any(str_detect(c(paper_mapped$methods,paper_mapped$introduction,paper_mapped$results,paper_mapped$discussion,paper_mapped$abstract), rob_regex)),
      
      Rob_guidelines=any(str_detect(c(paper_mapped$methods,paper_mapped$introduction,paper_mapped$results,paper_mapped$discussion,paper_mapped$abstract), guidelines_regex)),
      
      Rob_metaanalysis=str_detect(c(paper_mapped$methods), metaanalysis_regex),
      
      Rob_flowchart=any(str_detect(c(paper_mapped$methods,paper_mapped$introduction,paper_mapped$results,paper_mapped$discussion,paper_mapped$abstract), flowchart_regex)),
      
      Rob_coi=str_detect(c(paper_mapped$paper), coi_regex),
      
      Rob_deviations=str_detect(c(paper_mapped$results), deviation_regex),
      
      Rob_extraction=str_detect(c(paper_mapped$methods), induplicate_regex_one)
    ))
    
    papers_rob_df<-rbind(papers_rob_df,df)
  }, 
  silent = TRUE)
}


predicted<-papers_rob_df%>% 
  mutate(across(starts_with("Rob_"), ~ifelse(. == TRUE, 1, 0)))  %>% 
  mutate(across(starts_with("Rob_"), ~factor(., levels = c(0,1)))) %>% 
  mutate(file_name=gsub(" ","",file_name)) 

write.csv(predicted,"RoB_predicted.csv")
