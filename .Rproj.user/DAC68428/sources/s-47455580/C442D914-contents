# Managing data

tree <- ape::read.tree("Squam_100trees.tree")

traits <- read.csv("squamBase_cleaned.csv", header = TRUE)

traits <-  traits[, c("Scientific.Name", "MaxSVL_mm")]   

traits <- traits[!is.na(traits[, 2]), ]

trait <- traits[, 2]

names(trait) <- gsub(" ", "_", traits[, 1])

output <- list(length(tree))

models_by_tree <- character(length(tree))

# Run models for each tree

for(i in 1:length(tree)){

  print(i)
# Harmonize data
  
dat <- geiger::treedata(tree[[i]], trait)

# Fit models

lambda <- geiger::fitContinuous(phy = dat$phy,
                               dat = dat$data, model = "lambda")

eb <- geiger::fitContinuous(phy = dat$phy,
                               dat = dat$data, model = "EB")

OU <- geiger::fitContinuous(phy = dat$phy,
                            dat = dat$data, model = "EB")

kappa <- geiger::fitContinuous(phy = dat$phy,
                               dat = dat$data, model = "EB")

models_names <- c("lambda", "EB", "OU", "kappa")

models <- list(lambda, eb, OU, kappa)

names(models) <- models_names

lik <- sapply(models, function(m) m$opt$aicc)

names(lik) <- models_names

aicc_table <- geiger::aicw(lik)

best_model <- row.names(aicc_table[1, ])

models_by_tree[i] <- best_model

warp_tree <- geiger::rescale(tree[[i]], model = best_model)

tree_warpped <- warp_tree(models[[best_model]]$opt[[1]])

output[[i]] <- tree_warpped
}

ape::write.tree(output, paste0("warp_tree.tree"))

saveRDS(models_by_tree, "models_by_tree.rds")
