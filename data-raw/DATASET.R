## code to prepare `tb_clean` dataset goes here


tb <- read.csv("TBdataset-clean.csv", stringsAsFactors = FALSE)
tb_clean <- tb

## Libraries
library(tidyverse)
library(tibble)
library(lubridate)

colnames(tb_clean)<- toupper(colnames(tb))
to_date <- function(x){
  ifelse(grepl(pattern = "/|-| ", x = x), x, NA)
}
tb_clean <- tb %>% mutate(., co = gsub(" COUNTY", "", toupper(county))) %>%
                            mutate_at(vars(contains('DATE')), list(mdy))
table(tb_clean$co)
colnames(tb_clean)


tb_clean <- tb_clean %>%
  select(sex, county,  race, spsmear, hivstatus, homeless,
         PCR.Cluster,
         datecoll) %>%
    filter(PCR.Cluster != "") %>%
    group_by(PCR.Cluster) %>%
  rename(group = PCR.Cluster) %>%
    mutate(rel_time = datecoll -
               min(datecoll, na.rm = TRUE)) %>%
  dplyr::select(-datecoll)


usethis::use_data(tb_clean, overwrite = TRUE)


