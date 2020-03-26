library( tidyverse )

DFX <- read_csv( "metformin-dfx.csv" )
c2cp <- DRIAD::read_gmt( "c2.cp.v7.0.symbols.gmt" )
w <- with( DFX, set_names(logFC, Gene) )
gsea <- fgsea::fgsea( c2cp, w, 1e5 ) %>% as_tibble()
