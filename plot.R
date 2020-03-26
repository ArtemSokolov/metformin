library( tidyverse )
library( ggridges )

## Wrangle results
tmap <- list( AB = "A-v-B", AC = "A-v-C", BC = "B-v-C" )
X <- readRDS("metformin.RData") %>%
    select( -Set, -Feats, -Run ) %>%
    mutate( Label = str_c(" p=",pval),
           Label_x = case_when(
               Task == "AB" ~ 0.6,
               Task == "AC" ~ 0.8,
               Task == "BC" ~ 0.75 )) %>%
    mutate_at( "Task", recode, !!!tmap ) %>%
    mutate_at( "Dataset", as_factor )

## Unwrap the background
BK <- X %>% 
    select( Task, Dataset, AUC=BK ) %>%
    unnest(AUC)

gg <- ggplot( BK, aes(AUC, Dataset) ) +
    facet_wrap( ~Task, nrow=1, scales="free_x" ) +
    theme_ridges(center_axis_labels=TRUE) +
    geom_density_ridges2(scale=1.25, size=1, fill="gray", alpha=0.6) +
    geom_segment( aes(x=AUC, xend=AUC, y=as.numeric(Dataset),
                      yend=as.numeric(Dataset)+0.9),
                 data=X, color="tomato", lwd=2 ) +
    geom_text( aes(x=Label_x, y=as.numeric(Dataset)+0.7, label=Label),
              hjust=0, data=X, fontface="bold", size=4 )

ggsave( "metformin-DRIAD.pdf", gg, width=10, height=6 )
