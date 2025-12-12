# This suite of functions extracts healthcare records from the UK Biobank. You will need to have run create_tables.R once in your peojctbefore running this script

packages = c("dplyr", "tidyr", "readr")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
remotes::install_github("lcpilling/ukbrapR")
library('ukbrapR')

dxdownload=function(filepath){
  system(paste0("dx download ", filepath))
}
dld_files=ukbrapR:::ukbrapr_paths[,2]%>%
  lapply(basename)%>%
  lapply(file.exists)%>%
  unlist()
pathlist=ukbrapR:::ukbrapr_paths[!dld_files,2]
lapply(pathlist,dxdownload)


### define helper functions

collapse_codes=function(codes){
  # this function collapses a list of strings into one string with a | for use in a grep
  return(paste(codes, collapse = "|")) 
}

grep_codes=function(codes,filename){
  # this function builds a grep of the file using a code list (should be a list). It outputs the file to temp.tsv and attaches the original header back on
  grep_command=paste("grep -E", shQuote(collapse_codes(codes)), filename, " > temp.tsv")
  header_attach=paste("head -n 1 ", filename, " > temp2.tsv && cat temp.tsv >> temp2.tsv && mv temp2.tsv temp.tsv")
  system(grep_command)
  system(header_attach)
}

df_grep=function(codes,column){
  #find rows of a dataframe column that match on a grep, a bit like filter
  pattern <- paste0("^", collapse_codes(codes))
  return(grepl(pattern, column))
}

baseline_combo <- function(table,date_col){
  join_table=inner_join(table,baseline%>%select(eid, assess_date_initial,date_of_birth))
  join_table$prev=join_table[,date_col]<join_table[,'assess_date_initial']
  join_table$event_age=as.numeric(join_table[,date_col]-join_table[,'date_of_birth'])/365.25
  return(join_table)
}

### Define functions to read from the healthcare records

read_ICD10 <- function(codes){
  # The below block of code produces a grep command  which is used to create temp.tsv. 
  grep_codes(codes,'hesin_diag.tsv')
  # This awk runs a join at system level matching on the record id to create joined.tsv (but only includes the record column not the ICD10). This command was written with AI assistance 
  system("awk -F'\t' 'NR==FNR {keys[$2]; next} FNR==1 || $1 in keys' temp.tsv hesin.tsv > joined.tsv")
  diag_data=read.csv("temp.tsv",sep='\t')
  record_data=read.csv("joined.tsv",sep='\t')
  all_data=inner_join(diag_data,record_data,by = c('dnx_hesin_id', 'eid', 'ins_index'))
  all_data=all_data%>%select('eid','ins_index','arr_index','level','diag_icd9','diag_icd10','epistart','epiend')
  # Just in case, do a final pass to make sure that the ICD10 code came up in the diag_icd10 column
  all_data=all_data[df_grep(codes,all_data$diag_icd10),]
  all_data$epistart=as.Date(all_data$epistart)
  all_data$epiend=as.Date(all_data$epiend)
  system('rm joined.tsv')
  return(baseline_combo(all_data,'epistart'))
}

read_ICD9 <- function(codes){
  # The below block of code produces a grep command  which is used to create temp.tsv. 
  grep_codes(codes,'hesin_diag.tsv')
  # This awk runs a join at system level matching on the record id to create joined.tsv (but only includes the record column not the ICD10). This command was written with AI assistance 
  system("awk -F'\t' 'NR==FNR {keys[$2]; next} FNR==1 || $1 in keys' temp.tsv hesin.tsv > joined.tsv")
  diag_data=read.csv("temp.tsv",sep='\t')
  record_data=read.csv("joined.tsv",sep='\t')
  all_data=inner_join(diag_data,record_data,by = c('dnx_hesin_id', 'eid', 'ins_index'))
  all_data=all_data%>%select('eid','ins_index','arr_index','level','diag_icd9','diag_icd10','epistart','epiend')
  # Just in case, do a final pass to make sure that the ICD9 code came up in the diag_icd10 column
  all_data=all_data[df_grep(codes,all_data$diag_icd9),]
  all_data$epistart=as.Date(all_data$epistart)
  all_data$epiend=as.Date(all_data$epiend)
  system('rm joined.tsv')
  return(baseline_combo(all_data,'epistart'))
}

read_GP <- function(codes,table='clinical') {
  if (!(table%in%c('clinical','scripts'))){
    print('ERROR: table should be either clinical or scripts')
    return()
  }
  file=paste0('gp_',table,'.tsv')
  grep_codes(codes,file)
  data=read.csv("temp.tsv",sep='\t')
  if (table=='clinical'){
    data$event_dt=as.Date(data$event_dt)
    match_read2=df_grep(codes,data$read_2)
    match_read3=df_grep(codes,data$read_3)
    data=data[match_read2|match_read3,]
  }
  if (table=='scripts'){
    data$issue_date  =as.Date(data$issue_date)
    match_read2=df_grep(codes,data$read_2)
    match_bnf=data$bnf_code==codes
    match_dmd=data$dmd_code==codes
    data=data[match_read2|match_bnf|match_dmd,]
    
  }
  return(baseline_combo(data,'event_dt'))
}

read_OPCS <- function(codes) {
  grep_codes(codes,'hesin_oper.tsv')
  data=read.csv("temp.tsv",sep='\t')
  data$opdate=as.Date(data$opdate)
  data=data[df_grep(codes,data$oper4),]
  return(baseline_combo(data,'opdate'))
}

read_cancer <- function(codes,file='cancer_participant.csv') {
  codes <- gsub("\\.", "", codes) # cancer registry data doesn't have the . in the ICD10 code
  grep_codes(codes,'cancer_registry.tsv')
  data=read.csv('temp.tsv',sep='\t')
  long_data <- data %>%
    pivot_longer(
      cols = matches("^p\\d+_i\\d+$"),
      names_to = c("field", "instance"),
      names_pattern = "^p(\\d+)_i(\\d+)$",
      values_to = "value",
      values_transform = list(value = as.character),
      values_ptypes    = list(value = character())
    ) %>%
    mutate(
      field = paste0("p", field),         # restore the “p####” prefix for column names
      instance = as.integer(instance)
    ) %>%
    pivot_wider(
      names_from = field,
      values_from = value
    )
  names(long_data)=c('eid','instance','date','ICD10','age','histology','behaviour')
  long_data$date=as.Date(long_data$date)
  long_data=long_data[df_grep(codes,long_data$ICD10),]
  return(baseline_combo(long_data,'date'))
}

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/baseline_table.R')

print('Thank you for using ukbextractR Version 1.0, by Harry Green and Jiaqi Li, and ukbrapR by Luke Pilling, University of Exeter')

print('For any issues, please contact Harry Green at h.d.green@exeter.ac.uk')






