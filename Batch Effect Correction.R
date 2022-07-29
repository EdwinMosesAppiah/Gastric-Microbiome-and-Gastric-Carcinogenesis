library(bladderbatch)
library(Harman)
data(bladderdata)
bladderEset

total = median(sample_sums(phy_std))
## Function to standardize to median sample read count
stand_f = function(x, t=total) round(t * (x / sum(x)))

# Apply to phyloseq object
phy_std = transform_sample_counts(phy, stand_f)
ntaxa(phy_std)

sample_sums(phy_std)

Taxa_filt = filter_taxa(phy_std,function(x) sum(x > 10) > (0.2*length(x)) | sum(x) > 0.001*total, TRUE)

exprs(bladderEset)
edata <- otu_table(Taxa_filt)
sdata <- sample_data(Taxa_filt)
tdata <- tax_table(Taxa_filt)
expt <- sdata$Conditions
batch <- sdata$Project_ID
bladder_harman <- harman(edata, expt=expt, batch=batch)
summary(bladder_harman)
par(mfrow=c(1,1))
plot(bladder_harman)
plot(bladder_harman, 1, 5)

legend("topleft",c("PRJEB21497", "PRJEB26931", "PRJNA313391", "PRJNA375772", "PRJNA389357", "PRJNA413125"),
       cex=0.5,col=c("red","blue","green", "yellow", "brown", "purple"),
       pch=c(1,2,3,4,5,6))

library(limma, quietly=TRUE)
design <- model.matrix(~0 + expt)
colnames(design) <- c("Cancer", "Gastritis", "Healthy_control", "Metaplasia")
CorrectedBladderEset <- reconstructData(bladder_harman) 

CorrectedBladderEset[1000:1010, 1:10]  
CorrectedBladderEset <- as.matrix(CorrectedBladderEset)
edata[1:10,1:10]  
physeeq <- phyloseq(otu_table(CorrectedBladderEset,taxa_are_rows=TRUE),tax_table(tdata),sample_data(sdata))

