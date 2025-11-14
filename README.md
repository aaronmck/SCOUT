### SCOUT: Ornstein–Uhlenbeck modelling of gene expression evolution on single-cell lineage trees
Given a single-cell lineage tree, fit gene expression to evolutionary models to profile selection and drift. 

<img src="https://github.com/hrstuart/SCOUT/blob/main/data/F1_GraphicalAbstract.svg" alt="Graphical abstract of SCOUT" width="450" />

#### Quick Start
##### Installation
```
devtools::install('https://github.com/hrstuart/SCOUT/')
```
#### Run SCOUT
```
# Step 1: Prepare inputs 
inputs <- prepare_data(tree_path = tree, # path to newick file
                       metadata_path = meta, # path to metadata file 
                       outpath = outdir, # output directory 
                       regimes = c("BM1", "OU1", "OUM"), # substitute OUM for column names in metadata with regime annotation
                       normalize = TRUE, # if using counts
                       smoothing = 10) # smoothing param.

# Step 2: Fit models in parallel on prepared data
results <- fitModel(inputs,
                  cores = 4,
                  write = TRUE, 
                  outpath = outdir,
                  testgenes = c('BM1_1', 'OU1_1',  'OUM_1', )) # list of genes to test. Leave empty if want to test all in metadata.

# Step 3: Summarize and calculate fit metrics (AICc weight, delta AICc, etc.)
results <- get_fit_metrics(results, write = FALSE)

# Step 4: Perform model selection. Returns dataframe with final model per gene. 
final_results <- processes_results(results)

```

See examples: 
  - Simulation Vignette
  - Metastasis Vignette 

