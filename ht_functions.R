# Benjamin H Pepper
# BHPepper@gmail.com
# linkedin.com/in/BHPepper

library(data.table)
library(bnstruct)
library(glmnet)
library(MLmetrics)

get_ht_express_dat = function(root) {
  dat = fread(paste0(root, '\\data_expression_median.txt'), data.table=F)
  dat = dat[,-2]
  rownames(dat) = dat[,1]
  dat = dat[,-1]
  rownames(dat) = chartr('-', '_', rownames(dat))
  dat = t(dat)
  
  clin = fread(paste0(root, '\\data_clinical_patient.txt'), data.table=F, skip='PATIENT_ID')
  rownames(clin) = clin[,1]
  clin = clin[,-1]
  
  common = intersect(rownames(dat), rownames(clin))
  dat = dat[common,]
  clin = clin[common,]
  
  dat_imp = knn.impute(dat, k = 10)
  dat_imp = data.frame(dat_imp)
  
  dat_imp$'HORMONE_THERAPY' = as.factor(clin$`HORMONE_THERAPY`)
  
  return(dat_imp)
}

get_ht_clin_dat = function(root) {
  dat = fread(paste0(root, '\\data_expression_median.txt'), data.table=F)
  dat = dat[,-2]
  rownames(dat) = dat[,1]
  dat = dat[,-1]
  rownames(dat) = chartr('-', '_', rownames(dat))
  dat = t(dat)
  
  clin1 = fread(paste0(root, '\\data_clinical_patient.txt'), data.table=F, skip='PATIENT_ID')
  rownames(clin1) = clin1[,1]
  clin1 = clin1[,-1]
  
  clin2 = fread(paste0(root, '\\data_clinical_sample.txt'), data.table=F, skip='PATIENT_ID')
  rownames(clin2) = clin2[,1]
  clin2 = clin2[,-1]
  
  common = intersect(rownames(clin1), rownames(clin2))
  clin1 = clin1[common,]
  clin2 = clin2[common,]
  clin = cbind(clin1, clin2)
  
  common = intersect(rownames(dat), rownames(clin))
  dat = dat[common,]
  clin = clin[common,]
  
  return(clin)
}

lasso_mod = function(dat, response, params) {
  X = dat[,-which(colnames(dat) %in% response)]
  y = dat[,response]
  X = as.matrix(X)
  
  mod = glmnet(X, y, family='binomial', lambda = params, alpha=1)
  
  return(mod)
}

lasso_score = function(mod, dat, response) {
  X = dat[,-which(colnames(dat) %in% response)]
  X = as.matrix(X)
  
  preds = predict(mod, X, type = 'response')
  
  return(preds)
}

lasso_outer_perf = function(mod, dat, response) {
  preds = lasso_score(mod, dat, response)
  preds = ifelse(preds > .5, 'YES', 'NO')
  perf = Accuracy(dat[,response], preds)
  
  return(perf)
}

lasso_inner_perf = function(mod, dat, response) {
  perf = lasso_inner_perf(mod, dat, response)
  
  return(perf)
}

lr_mod = function(dat, response, params) {
  f = as.formula(paste0(response, ' ~ .'))
  mod = glm(f, dat, family='binomial')
  
  return(mod)
}

lr_score = function(mod, dat, response) {
  preds = predict(mod, dat, type = 'response')
  
  return(preds)
}

lr_outer_perf = function(mod, dat, response) {
  preds = lr_score(mod, dat, response)
  preds = ifelse(preds > .5, 'YES', 'NO')
  perf = Accuracy(dat[,response], preds)
  
  return(perf)
}

integrated_mod = function(dat, response, params) {
  mod1 = lasso_mod(dat[[1]], response, params)
  dat_integrated = dat[[2]]
  dat_integrated$GENOMIC_SCORE = lasso_score(mod1, dat[[1]], response)[,1]
  mod2 = lr_mod(dat_integrated, response)
  mod = list(mod1 = mod1, mod2 = mod2)

  return(mod)
}

integrated_score = function(mod, dat, response) {
  dat_integrated = dat[[2]]
  dat_integrated$GENOMIC_SCORE = lasso_score(mod$mod1, dat[[1]], response)[,1]
  preds = predict(mod$mod2, dat_integrated, type='response')
  
  return(preds)
}

integrated_outer_perf = function(mod, dat, response) {
  preds = integrated_score(mod, dat, response)
  preds = ifelse(preds > .5, 'YES', 'NO')
  perf = Accuracy(dat[[2]][,response], preds)
  
  return(perf)
}

integrated_inner_perf = function(mod, dat, response) {
  perf = integrated_outer_perf(mod, dat, response)
  
  return(perf)
}
