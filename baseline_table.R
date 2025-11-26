# Jiaqi Li
# Create a baseline table 
baseline <- read.csv("baseline_dates.tsv", sep = "\t")
baseline <- baseline %>%
  rename('assess_date_initial' = p53_i0,
         'assess_date_fist_repeat' = p53_i1,
         'imaging_visit_date' = p53_i2,
         'imaging_visit_first_repeat' = p53_i3,
         'centre_id_initial' = p54_i0,
         'centre_id_first_repeat' = p54_i1,
         'centre_id_imaging_visit' = p54_i2,
         'centre_id_imaging_visit_first_repeat' = p54_i3,
         'month_of_birth' = p52,
         'year_of_birth' = p34,
         'sex' = p31)

## Add date of birth
baseline <- baseline %>%
  arrange(eid) %>%
  mutate(date_of_birth = as.Date(paste(year_of_birth, month_of_birth,
                                       16, sep = '/')))

baseline <- baseline %>%
  mutate(assess_date_initial = as.Date(assess_date_initial),
         assess_date_fist_repeat = as.Date(assess_date_fist_repeat),
         imaging_visit_date = as.Date(imaging_visit_date),
         imaging_visit_first_repeat = as.Date(imaging_visit_first_repeat))

## Calculate the age of assessement
baseline <- baseline %>%
  mutate(assess_age = as.numeric(assess_date_initial - date_of_birth)/365.25,
         sex = ifelse(sex == 1, 'Male', 'Female'))
      
