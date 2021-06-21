## SKG
## June 21, 2021
## Revisions for JCGS
## Tables 3 and 4 in manuscript
## Looking at distribution of cluster sizes
## and a cross section of part of the data (cluster 15)

devtools::load_all()
library(tidyverse)
library(knitr)
library(kableExtra)
data(tb_clean)

## Table 3


### grouped individuals
clust_sizes <- tb_clean %>%
    group_by(group) %>%
    summarize(size = n())

tab <- data.frame(table(clust_sizes$size))
colnames(tab) <- c("Cluster Size", "Freq.")

kable(t(tab),
      format = "latex", row.names = TRUE,
      booktabs = TRUE,
      caption = "Distribution of cluster size within the Maryland 2003-2009 TB data.",
      label = "data-clust-sizes")



## Table 4

## Example cluster
tb_clean %>% filter( group == 15) %>%
    select(group, county, rel_time, hivstatus, race, sex, spsmear) %>%
    arrange(datecoll) %>%
    kable(format = "latex",
          col.names = c("Cluster ID", "Co.", "Relative day",
                        "HIV status", "Race", "Sex", "Smear"),
          booktabs = TRUE,
          caption = "Example of individuals and some of their features within a single cluster", label = "ex-clust") %>%
    kable_styling(latex_options="scale_down")
          
