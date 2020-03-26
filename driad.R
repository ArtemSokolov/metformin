## Run the DGL through DRIAD
library( tidyverse )
library( DRIAD )

## wrangleROSMAP("./amp-ad/rosmap")
## wrangleMSBB("./amp-ad/msbb")

## Identify individual datasets
fn <- list(ROSMAP = "amp-ad/rosmap/rosmap.tsv.gz",
           MSBB10 = "amp-ad/msbb/msbb10.tsv.gz",
           MSBB22 = "amp-ad/msbb/msbb22.tsv.gz",
           MSBB36 = "amp-ad/msbb/msbb36.tsv.gz",
           MSBB44 = "amp-ad/msbb/msbb44.tsv.gz")

## Prediction tasks to consider
XX <- enframe( fn, "Dataset", "Filename" ) %>%
    crossing( Task = list("AB", "BC", "AC") ) %>%
    unnest( c(Filename, Task) ) %>%
    mutate( Run = str_c(Dataset, "_", Task),
           Data = map2(Filename, Task, prepareTask) )

## Generate pairs for LPOCV
XP <- XX %>% select( -Filename ) %>%
    mutate( Pairs = map(Data, preparePairs) )

## Load the DGL
dgl <- scan("metformin-dgl.txt", what=character())
f <- function(.x, .y ) evalGeneSets( list(Metformin=dgl), .x, .y, 100 )

## Evaluate the gene set across all datasets / tasks
future::plan( future::multiprocess )
set.seed(100)
fo <- furrr::future_options( seed=TRUE )
RR <- XP %>%
    mutate( Result = furrr::future_map2(Data, Pairs, f, .options=fo) )

## Finalize results and write to file
RS <- RR %>% select( -Data, -Pairs ) %>% unnest( Result )
saveRDS( RS, file="metformin.RData" )
