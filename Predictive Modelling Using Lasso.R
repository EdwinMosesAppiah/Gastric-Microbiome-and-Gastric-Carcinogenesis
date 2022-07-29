phy_std <- readRDS("/media/edwin/Transcend/datasets/dada2/meta_analysis/network/NET/Mf_processed_phyloseq.RDS")
setwd("/media/edwin/Transcend/datasets/dada2/meta_analysis/ml/new")

library(SIAMCAT)
data("feat_crc_zeller", package="SIAMCAT")
data("meta_crc_zeller", package="SIAMCAT")
taxa_names(phy_std) <- tax_table(phy_std)[,"Species"]

condition1 <- "Metaplasia"
condition2 <- "Cancer"
doc <- paste(condition1,condition2, sep = "_")
# c(condition1,condition2)

phy_siam <- subset_samples(phy_std, Conditions %in% c(condition1, condition2))

###convert to relative abundance

physeq_rel_abund = phyloseq::transform_sample_counts(phy_siam, function(x){x / sum(x)})
# phyloseq::otu_table(phy_siam)[1:15, 1:15]
# phyloseq::otu_table(physeq_rel_abund)[1:15, 1:15]

label <- create.label(meta=sample_data(physeq_rel_abund), label="Conditions", case = condition2)

siamcat <- siamcat(phyloseq=physeq_rel_abund, label=label)

sc.obj <- filter.features(siamcat, cutoff=1e-04,
                          filter.method = 'abundance')

## Features successfully filtered
sc.obj <- filter.features(sc.obj, cutoff=0.05,
                          filter.method='prevalence',
                          feature.type = 'filtered')
## Features successfully filtered
sc.obj <- check.associations(sc.obj, detect.lim = 1e-06,
                             alpha=0.001, 
                             max.show = 50,
                             plot.type = 'quantile.rect',
                             mult.corr = 'fdr',
                             sort.by = 'fc',
                             panels = c("fc", "prevalence", "auroc"),
                             verbose = 3,
                             fn.plot = paste("association_",doc,".pdf", sep = ""))

# sc.obj <- check.confounders(
#   sc.obj,
#   fn.plot = 'confounder_HG_plots.pdf',
#   meta.in = NULL,
#   feature.type = 'filtered'
# )

sc.obj <- normalize.features(
  sc.obj,
  norm.method = "log.unit",
  norm.param = list(
    log.n0 = 1e-06,
    n.p = 2,
    norm.margin = 1
  )
)


sc.obj <-  create.data.split(
  sc.obj,
  num.folds = 10,
  num.resample = 20
)

Model_trained <- train.model(
  sc.obj,
  method = "randomForest",
  param.set = 300
)

# get information about the model type
model_type(Model_trained)

## [1] "lasso"

# access the models
models <- models(Model_trained)
models[[1]]

## Model for learner.id=classif.cvglmnet; learner.class=classif.cvglmnet
## Trained on: task.id = data; obs = 112; features = 207
## Hyperparameters: nlambda=100,alpha=1
?make.predictions

sc.obj_try <- make.predictions(Model_trained, sc.obj, normalize.holdout=TRUE)
pred_matrix <- pred_matrix(sc.obj)
head(pred_matrix)
model_list(sc.obj_try) <- model_list(Model_trained) 

sc.obj <-  evaluate.predictions(sc.obj_try)

## Evaluated predictions successfully.

model.evaluation.plot(sc.obj)


model.interpretation.plot(
  sc.obj,
  fn.plot = paste("intepretation_",doc,".pdf", sep = ""),
  consens.thres = 0.5,
  limits = c(-3, 3),
  heatmap.type = 'zscore',
)



