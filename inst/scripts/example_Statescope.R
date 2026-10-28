# This script was used to create the `example_scRNAseq.RData`, `example_signature.RData`,
# `example_selected_genes.RData`, `example_Statescope_Deconvolved.RData`,
# `example_Statescope_Refined.RData` and `example_Statescope_Discovered.RData` files.
# This file contains an example Statescope pipeline of which the intermediate results
# are saved for use in examples and tests.
# For this purpose the SegerstolpePancreasData from the scRNAseq package was used.
# Rare celltypes were excluded, leaving the 5 most common celltypes for analysis.
# After this exclusion, standard preprocessing was done, after which genes were
# selected with AutoGeneS and saved, before running the Statescope framework and saving all steps.
# Key package versions:
#
# scRNAseq   v2.23.0

library(StatescopeR)
library(scRNAseq)
library(scuttle)

## Load SegerstolpePancreas data set
scRNAseq <- SegerstolpePancreasData()

## remove duplicate genes
scRNAseq <- scRNAseq[!duplicated(rownames(scRNAseq)), ]

## Subset to 1 healthy and 2 type 2 diabetes samples
scRNAseq = scRNAseq[,scRNAseq$individual %in% c('H3',
                                                'T2D1','T2D2')]
## remove cells with no cell type label
scRNAseq <- scRNAseq[, !is.na(scRNAseq$`cell type`)]

## remove rare cell types (<100 cells in total data set)
celltypes_to_remove <-
    names(table(scRNAseq$`cell type`)[(table(scRNAseq$`cell type`) < 100)])
scRNAseq <- scRNAseq[, !scRNAseq$`cell type` %in% celltypes_to_remove]

## Gather true fractions and save
true_fractions <- gather_true_fractions(scRNAseq, ids = scRNAseq$individual,
                                        label_col = 'cell type')
save(true_fractions, file = 'inst/extdata/example_true_fractions.RData')

## Normalize (cp10k) and logtransform scRNAseq
cpm(scRNAseq) <- calculateCPM(scRNAseq)
logcounts(scRNAseq) <- log1p(cpm(scRNAseq)/100)

## Create pseudobulk and normalize to cp10k (logging is done within Statescope)
pseudobulk <- aggregateAcrossCells(scRNAseq, ids = scRNAseq$individual)
normcounts(pseudobulk) <- calculateCPM(pseudobulk)/100
pseudobulk = as(pseudobulk, "SummarizedExperiment")
rownames(pseudobulk) = rownames(scRNAseq)

## Create scRNAseq reference/signature with 5 hvg for quick example
signature <- create_signature(scRNAseq, hvg_genes = TRUE, n_hvg_genes =  5L,
                              labels = scRNAseq$`cell type`)

save(signature, file = 'inst/extdata/example_signature.RData')

## select subset of genes for deconvolution (3/5 hvg to make it quick)
selected_genes <- select_genes(scRNAseq, 3L, 5L,
                               labels = scRNAseq$`cell type`)

save(selected_genes, file = 'inst/extdata/example_selected_genes.RData')

## (optional) Create prior expectation using True sc fractions
prior <- gather_true_fractions(scRNAseq,
            ids = scRNAseq$individual, label_col = 'cell type')
prior[rownames(prior) != "ductal cell", ] <- NA #Keep only ductal cells as prior
prior <- t(prior) # Transpose it to nSample x nCelltype

save(prior, file = 'inst/extdata/example_prior.RData')

## Perform Deconvolution with BLADE, refine gene expression estimates
Statescope <- BLADE_deconvolution(
    signature, pseudobulk, selected_genes,
    prior, 1L, Nrep = 1L
)

## Save to RData
save(Statescope, file = 'inst/extdata/example_Statescope_Deconvolved.RData')

Statescope <- Refinement(Statescope, signature, pseudobulk, 2L)

## Save to RData
save(Statescope, file = 'inst/extdata/example_Statescope_Refined.RData')

## Discover states
Statescope <- StateDiscovery(Statescope, k=2L, Ncores = 2L)

## Save to RData
save(Statescope, file = 'inst/extdata/example_Statescope_StateDiscovered.RData')
