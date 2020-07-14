## SKG
## June 16, 2020
## EDA on data for paper

## size, type

devtools::load_all()
library(tidyverse)
library(knitr)
library(kableExtra)

data(tb_clean)

dim(tb_clean)
range(tb_clean$datecoll, na.rm = TRUE)

### grouped individuals
tb_sub <- tb_clean %>% filter(PCR.Cluster != "")
dim(tb_sub)

clust_sizes <- tb_sub %>%
    group_by(PCR.Cluster) %>%
    summarize(size = n())

tab <- data.frame(table(clust_sizes$size))
colnames(tab) <- c("Cluster Size", "Freq.")

kable(t(tab),
      format = "latex", row.names = TRUE,
      booktabs = TRUE,
      caption = "Distribution of cluster sizes.",
      label = "data-clust-sizes")


### Single cluster
tb_join <- left_join(tb_sub, clust_sizes,
                     by = "PCR.Cluster")

## Example cluster
tb_join %>% filter(size == 5 &
                               PCR.Cluster == 15) %>%
    select(PCR.Cluster, co, datecoll, hivstatus, race, sex, spsmear) %>%
    arrange(datecoll) %>%
    kable(format = "latex",
          col.names = c("Cluster ID", "Co.", "Collection Date",
                        "HIV", "Race", "Sex", "Smear"),
          booktabs = TRUE,
          caption = "Example of individuals and some of their features within a single cluster", label = "ex-clust") %>%
    kable_styling(latex_options="scale_down")
          
