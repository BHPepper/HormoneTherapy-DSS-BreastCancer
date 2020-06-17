library(tidyverse)
library(fastDummies)

source("ht_functions.R")
source("ht_nested_cv.R")

# Location of data
root = "..\\Data\\brca_metabric"

# Get Genomic Expression and Clinical Data
express = get_ht_express_dat(root)
express$y = as.factor(express$'HORMONE_THERAPY')
express = express[,!colnames(express) %in% 'HORMONE_THERAPY']

clin = get_ht_clin_dat(root)

# Exclude Cheating Variables
var_remove = c('ER_STATUS', 'ER_IHC', 'PR_STATUS', 'HER2_STATUS', 'HORMONE_THERAPY',
               'INTCLUST', 'VITAL_STATUS', 'HISTOLOGICAL_SUBTYPE', 'OS_STATUS', 
               'OS_MONTHS', 'SAMPLE_ID', 'THREEGENE', 'CANCER_TYPE',
               'CANCER_TYPE_DETAILED', 'SAMPLE_TYPE')
clin_reduced = na.omit(clin[,!colnames(clin) %in% var_remove])

# Dummification
clin_rownames = rownames(clin_reduced)
clin_dummies = dummy_cols(clin_reduced, remove_selected_columns = T, remove_first_dummy=T)
rownames(clin_dummies) = clin_rownames

# Effectively Merge Data sets by PATIENT_ID (rownames)
# Data sets kept separate as elements of a list
common = intersect(rownames(clin_dummies), rownames(express))
clin = clin_dummies[common,]
express = express[common,]
clin$y = express$y

dat = list(express, clin)

# Train & Test Integrated Model
integrated_grid = data.frame(lambda = seq(0.000001, 1, by = 0.05))
res_integrated = nested_cv(cv_k1 = 10, cv_k2 = 10, seed = 1, model = integrated_mod,
                           inner_perf_f = integrated_inner_perf, 
                           outer_perf_f = integrated_outer_perf,
                           score_f = integrated_score, perf_type = 'high',
                           grid = integrated_grid, dat = dat, response = 'y')
res_integrated$Performance
mean(res_integrated$Performance$Perf)

# Build Final Model from Best Parameters
best_integrated_params = mean(integrated_grid[res_integrated$Performance$BestParams,])
mod_integrated = integrated_mod(dat, response = 'y', best_integrated_params)
clin_reduced$GENOMIC_SCORE = lasso_score(mod_integrated$mod1, express, response = 'y')[,1]
clin_reduced <- clin_reduced %>% select(GENOMIC_SCORE, everything())

# Save RDS
saveRDS(mod_integrated, 'integrated_mod.RDS')
saveRDS(res_integrated$Scores, 'integrated_scores.RDS')
saveRDS(res_integrated$Y, 'integrated_labels.RDS')
saveRDS(clin_reduced, 'integrated_dat.RDS')