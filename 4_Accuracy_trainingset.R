Translation_extraction <- read_excel("Translation-extraction.xlsx",
    sheet = "Analysis") #Manual extraction

actual<-
 Translation_extraction %>%
 filter(!(is.na(Title)|is.na(file))) %>%
 rename(Rob_metaanalysis=Rob_synthesis) %>%
 select(!c(Rob_inclusion)) %>%
 select(file,starts_with("Rob_")) %>%
 mutate(across(starts_with("Rob_"), ~as.numeric(as.character(.)))) %>%
 mutate(across(starts_with("Rob_"), ~gsub("2","0",.))) %>%
 # mutate(across(starts_with("Rob_"), ~ifelse(. == TRUE, 1, 0))) %>%
 mutate(across(starts_with("Rob_"), ~factor(., levels = c(0,1)))) %>%
 mutate(file=gsub(" ","",file)) %>%
 merge(.,predicted[1],by.x="file",by.y="file_name")

#Calculating the accuracy measures for all the items mined

df_confmat<-data.frame()

for(i in names(actual[-1])){
  df_template<-data.frame(cbind(result=c("TN","TP","FN","FP"),na=NA))

df_temp<-
  actual %>% select(i) %>% rename(actual=1) %>%
    cbind(.,predicted %>% select(i))%>% rename(predicted=2) %>%
  drop_na(actual) %>%
  mutate(result=case_when(
    actual==1&predicted==1~"TP",
    actual==1&predicted==0~"FN",
    actual==0&predicted==1~"FP",
    actual==0&predicted==0~"TN",
  )) %>%
  group_by(result) %>% count() %>%
  merge(.,df_template,by="result",all.y=T) %>%
  select(!(na)) %>%
  replace_na(list(n=0)) %>%
  pivot_wider(names_from = result,values_from = n)

df<-data.frame(cbind(
  item=i,
  df_temp,
  sensitivity=df_temp$TP/(df_temp$TP+df_temp$FN),
  specificity=df_temp$TN/(df_temp$TN+df_temp$FP),
  precision=df_temp$TP/(df_temp$TP+df_temp$FP),
  F1=2*df_temp$TP/(2*df_temp$TP+df_temp$FP+df_temp$FN),
    accuracy=(df_temp$TP+df_temp$TN)/(df_temp$TP+df_temp$TN+df_temp$FP+df_temp$FN)
))
  df_confmat<-rbind(df_confmat,df)
}


df_confmat%>%
  mutate_all(~gsub("NaN", NA, .)) %>%
 mutate(across(c(where(is.character), -item), as.numeric)) %>%
  View()