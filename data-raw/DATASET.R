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

<<<<<<< HEAD

tb_clean <- tb_clean %>%
  select(sex, county,  race, spsmear, hivstatus, homeless,
         PCR.Cluster,
         datecoll) %>%
    filter(PCR.Cluster != "") %>%
    group_by(PCR.Cluster) %>%
    mutate(rel_time = datecoll -
               min(datecoll, na.rm = TRUE))
=======
library(tidveryse)
tb_clean <- tb_clean %>%
  select(sex, county, ageatrept, race, spsmear, hivstatus, homeless,
         INIT_REGIMEN_START_DATE,
         PCR.Cluster,
         datecoll)
>>>>>>> fd400a1555a7a782265dd95dc3e8cb9f965adce0

usethis::use_data(tb_clean, overwrite = TRUE)


<<<<<<< HEAD
=======
## Clusters
library(dplyr)


tb_clean$hiv <- ifelse(tb_clean$hivstatus == "Negative", -1,
                ifelse(tb_clean$hivstatus == "Positive", 1,
                       0))
clusts_tb <- tb_clean %>% filter(PCR.Cluster != "") %>%
    group_by(PCR.Cluster) %>%
    summarize(size = dplyr::n(),
              n_pos = sum(spsmear == "Positive"),
              n_neg = sum(spsmear == "Negative"),
              n_unk = sum(spsmear == "Unknown"),
              n_hiv_pos = sum(hiv == 1),
              n_hiv_neg = sum(hiv == -1),
              n_hiv_unk = sum(hiv == 0),
              inf_range = diff(range(INIT_REGIMEN_START_DATE, na.rm = TRUE)),
              inf_iqr = IQR(INIT_REGIMEN_START_DATE, na.rm= TRUE)
              )



usethis::use_data(clusts_tb, overwrite = TRUE)
>>>>>>> fd400a1555a7a782265dd95dc3e8cb9f965adce0
