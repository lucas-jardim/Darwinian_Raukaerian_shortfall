
library(PhylogeneticEM)

data(monkeys)
## Run method
# Note: use more alpha values for better results.
res <- PhyloEM(Y_data = monkeys$dat[2,],        ## data
               phylo = monkeys$phy,         ## phylogeny
               process = "BM",            ## scalar OU
               random.root = TRUE,          ## root is stationary
               stationary.root = TRUE,
               K_max = 10,                  ## maximal number of shifts
               nbr_alpha = 4,               ## number of alpha values
               parallel_alpha = TRUE,       ## parallelize on alpha values
               Ncores = 2)
