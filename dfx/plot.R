library( tidyverse )

load( "met-dfx.RData" )

## Identify the genes of interest
DFX <- Rmet %>%
    filter(FDR < 0.005,
           sign(logFC40) == sign(logFC10),
           abs(logFC40) > abs(logFC10))

## Isolate the corresponding slice of data
Y <- Ymet %>% rownames_to_column( "Sample" ) %>% as_tibble()
X <- Xmet %>% rownames_to_column( "Gene" ) %>% as_tibble()

RR <- X %>% filter( Gene %in% DFX$Gene ) %>%
    mutate( Med = pmap_dbl(list(S46, S47, S48), lift_vd(median)) ) %>%
    gather( Sample, Counts, -Gene, -Med ) %>%
    inner_join( Y, by="Sample" ) %>% 
    mutate(logFC = log(Counts / Med),
           logCnt = log10(Counts),
           Sample = str_c(Sample, " ", Treatment))

RX <- DFX %>% select( Gene, logFC10, logFC40 ) %>%
    gather( Comparison, logFC, -Gene ) %>%
    mutate_at( "Comparison", recode,
              logFC10 = "Metformin (10uM)",
              logFC40 = "Metformin (40uM)" )

## Plotting elements
pal1 <- rev(RColorBrewer::brewer.pal(n=7, name="YlOrRd"))
pal2 <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))

gg1 <- ggplot( RR, aes(x=Sample, y=Gene, fill=logCnt) ) +
    theme_minimal() + geom_tile() +
    scale_fill_gradientn( colors=pal1, name="log10(Counts)" ) +
    theme( axis.text.x = element_text( angle=90, vjust=0.5, hjust=1 ),
          axis.text.y = element_text(size=9) )

gg2 <- ggplot( RR, aes(x=Sample, y=Gene, fill=logFC) ) +
    theme_minimal() + geom_tile() +
    scale_fill_gradientn( colors=pal2, name="logFC", limits=c(-1.2, 1.2) ) +
    theme( axis.text.x = element_text( angle=90, vjust=0.5, hjust=1 ) )

gg3 <- ggplot( RX, aes(x=Comparison, y=Gene, fill=logFC) ) +
    theme_minimal() + geom_tile() +
    scale_fill_gradientn( colors=pal2, name="logFC", limits=c(-0.85, 0.85) ) +
    theme( axis.text.x = element_text( angle=90, vjust=0.5, hjust=1 ) )
    
gg <- cowplot::plot_grid( gg1, gg2, gg3, ncol=3, rel_widths=c(1.1,1,0.7),
                         labels=str_c("Demo", 1:3), label_x=c(0.6, 0.7, 0.4) )
ggsave( "test.pdf", gg, width=9, height=10 )
