setwd("/media/edwin/Transcend/datasets/dada2/meta_analysis/network/NET")
library(NetCoMi)
library(SpiecEasi)
data("amgut1.filt")
data("amgut2.filt.phy")


phy <- readRDS("/media/edwin/Transcend/datasets/dada2/meta_analysis/meta_analysis_phyloseq.RDS")
TREE = read_tree("/media/edwin/Transcend/datasets/dada2/meta_analysis/picrust/meta.tre")
phy_tree(phy) <- TREE

# Check the number of reads per sample
reads <- sample_sums(phy) 
head(reads)
table(tax_table(phy)[, "Kingdom"], exclude = NULL)
ps <- subset_taxa(phy, !is.na(Species) & !Phylum %in% c("", "uncharacterized"))

filterKingdom = c("Archaea", "NA")
phy = subset_taxa(ps, !Kingdom %in% filterKingdom)


# Check the number of reads per sample
reads <- sample_sums(phy) 
reads

# Standardize sample read count to compare between samples
## Set (arbitrary) cutoff for number of acceptable reads/sample
length(which(reads<5000))
## Find median sample read count
total = median(sample_sums(phy))
## Function to standardize to median sample read count
stand_f = function(x, t=total) round(t * (x / sum(x)))

# Apply to phyloseq object
phy_std = transform_sample_counts(phy, stand_f)
ntaxa(phy_std)
ps.taxa <- tax_glom(phy_std, taxrank = 'Species', NArm = TRUE)
M.f = filter_taxa(ps.taxa,function(x) sum(x > 100) > (0.2*length(x)) | sum(x) > 0.01*total, TRUE)
ntaxa(M.f)

taxa_names(M.f) <- tax_table(M.f)[,"Species"]

saveRDS(M.f, "Mf_processed_phyloseq.RDS") 

ps.taxa.Gas <- subset_samples(M.f, Conditions %in% c("Gastritis"))
ps.taxa.Can <- subset_samples(M.f, Conditions %in% c("Cancer"))
ps.taxa.Met <- subset_samples(M.f, Conditions %in% c("Metaplasia"))
ps.taxa.HC <- subset_samples(M.f, Conditions %in% c("Healthy_control"))

write.table(otu_table(ps.taxa.Can), "Cancer.txt", sep = "\t")
write.table(sample_data(ps.taxa.Can), "Cancersample.txt", sep = "\t")
#####SpiecEasi run
#spiec.out=spiec.easi(ps.taxa.Can, method="mb",icov.select.params=list(rep.num=20))

#spiec.graph=adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(M.f)))
#plot_network(spiec.graph, M.f, type='taxa', color="Phylum", label=NULL)

#betaMat=as.matrix(symBeta(getOptBeta(spiec.out)))

#clusters=cluster_fast_greedy(spiec.graph)
#clusterOneIndices=which(clusters$membership==1)
#clusterOneOtus=clusters$names[clusterOneIndices]
#clusterTwoIndices=which(clusters$membership==2)
#clusterTwoOtus=clusters$names[clusterTwoIndices]

net_spring_single <- netConstruct(ps.taxa.Can, verbose = 3,
                                  filtTax = "highestFreq",
                                  filtTaxPar = list(highestFreq = 100),
                                  filtSamp = "totalReads",
                                  filtSampPar = list(totalReads = 1000),
                                  zeroMethod = "none", normMethod = "none",
                                  measure = "spring",
                                  measurePar = list(nlambda = 20, rep.num = 20,
                                                    ncores=12),
                                  sparsMethod = "none",  dissFunc = "signed", 
                                  seed = 20190101)

### Network analysis
netprops_spring_single <- netAnalyze(net_spring_single, 
                                     clustMethod = "cluster_fast_greedy",
                                     hubPar = "eigenvector", hubQuant = 0.95)


#------------------------------------------------------
# Identify strongest positive and negative correlations
assomat <- net_spring_single$assoMat1
assomat[1:5, 1:5]
diag(assomat) <- 0

max(assomat)
which(assomat == max(assomat), arr.ind = TRUE)

min(assomat)
which(assomat == min(assomat), arr.ind = TRUE)
library(dplyr)
library(tidyr)
assomat1 <- data.frame(assomat) %>%
        rownames_to_column() %>%
        gather(key="variable", value="correlation", -rowname) %>%
        filter(abs(correlation) > 0.4)

assomat1[1:5,1:3]

ggplot(data =assomat1, aes(x=rowname, y=variable, fill = correlation)) +
        geom_tile()+ 
        scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                             midpoint = 0, limit = c(-1,1), space = "Lab", 
                             name="Correlation") +theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))


write.csv(assomat, "assomat_Can.csv")

#------------------------------------------------------
# Table S5
capture.output(summary(netprops_spring_single), file = "GAS.txt")

p <- plot(netprops_spring_single, 
          repulsion = 0.86,
          shortenLabels = "simple",
          labelLength = 40L,
          charToRm = "(",
          labelScale = FALSE,
          nodeFilter = "none",
          nodeFilterPar = 20,
          rmSingles = "inboth",
          nodeSize = "eigen",
          nodeSizeSpread = 4,
          nodeColor = "cluster",
          #colorVec = rainbow(6),
          nodeTransp = 65,
          hubTransp = 50,
          hubBorderWidth = 2,
          hubBorderCol = "gray40",
          colorNegAsso = TRUE,
          edgeTranspLow = 70,
          edgeTranspHigh = 30,
          cexNodes = 1,
          cexHubs = 1.3,
          showTitle = FALSE,
          mar = c(1,3,1,3))


net <- netConstruct(
        ps.taxa.pse.sub,
        dataType = "counts",
        group = "Conditions",
        matchDesign = NULL,
        measure = "spieceasi",
        measurePar = NULL,
        jointPrepro = TRUE,
        filtTax = "relFreq",
        filtTaxPar = 0.3,
        filtSamp = "none",
        filtSampPar = NULL,
        zeroMethod = "none",
        zeroPar = NULL,
        normMethod = "clr",
        normPar = NULL,
        sparsMethod = "softThreshold",
        thresh = 0.6,
        alpha = 0.01,
        adjust = "adaptBH",
        trueNullMethod = "lfdr",
        lfdrThresh = 0.2,
        nboot = 1000L,
        cores = 12,
        logFile = "log.txt",
        softThreshType = "signed",
        softThreshPower = NULL,
        softThreshCut = 0.8,
        kNeighbor = 3L,
        knnMutual = FALSE,
        dissFunc = "signed",
        dissFuncPar = NULL,
        simFunc = NULL,
        simFuncPar = NULL,
        scaleDiss = TRUE,
        weighted = TRUE,
        sampleSize = NULL,
        verbose = 2,
        seed = NULL
)

net_single2 <- netConstruct(ps.taxa.Gas,  
                            measure = "spring",
                            normMethod = "clr", 
                            measurePar = "kcld",
                            zeroMethod = "multRepl",
                            sparsMethod = "threshold", 
                            thresh = 0.3,
                            cores = 12,
                            verbose = 3)
summary(net_single2)
props_single2 <- netAnalyze(net_single2, clustMethod = "cluster_fast_greedy")

summary(netprops_spring_single)

plot(netprops_spring_single, 
     layout = "spring",
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     nodeFilter = "highestDegree",
     charToRm = NULL,
     nodeFilterPar = 50,
     hubBorderCol = "black",
     labelScale = FALSE,
     cexNodes = 2, 
     labelLength = 40,
     cexLabels = 2,
     cexHubLabels = 3,
     cexTitle = 5,
     nodeTransp = 80,
     title1 = "Network of Species in Healthy Control Individuals", 
     showTitle = TRUE) 
legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:", 
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)




####Network comparison
ps.taxa.pse.sub <- subset_samples(M.f, Conditions %in% c("Metaplasia", "Cancer"))
amgut_split <- metagMisc::phyloseq_sep_variable(ps.taxa.pse.sub, 
                                                "Conditions")

# Network construction
net_season <- netConstruct(data = amgut_split$Metaplasia, 
                           data2 = amgut_split$Cancer,  
                           filtTax = "highestVar",
                           filtTaxPar = list(highestVar = 200),
                           measure = "spring",
                           measurePar = list(nlambda=10, 
                                             rep.num=10),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 3,
                           seed = 123456)

assomat_compare_MC <- net_season$assoMat1
assomat_compare_CM <- net_season$assoMat2
write.csv(assomat_compare_MC, "assomat_compare_MC.csv")
write.csv(assomat_compare_CM, "assomat_compare_CM.csv")

props_season <- netAnalyze(net_season, 
                           centrLCC = FALSE,
                           avDissIgnoreInf = TRUE,
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree", "between", "closeness"),
                           hubQuant = 0.9,
                           lnormFit = TRUE,
                           normDeg = FALSE,
                           normBetw = FALSE,
                           normClose = FALSE,
                           normEigen = FALSE)


capture.output(summary(props_season), file = "HC_CAN.txt")
plot(props_season, 
     sameLayout = TRUE, 
     layoutGroup = 1,
     rmSingles = "inboth", 
     nodeColor = "cluster",
     nodeSize = "mclr",
     labelScale = FALSE,
     cexNodes = 5, 
     labelLength = 100,
     cexLabels = 2,
     cexHubLabels = 3,
     cexTitle = 5,
     nodeFilter = "highestDegree",
     nodeFilterPar = 20,
     groupNames = c("Gastritis", "Metaplasia"),
     hubBorderCol  = "gray40")

legend("bottom", title = "estimated association:", legend = c("+","-"), 
       col = c("#009900","red"), inset = 0.02, cex = 4, lty = 1, lwd = 4, 
       bty = "n", horiz = TRUE)

diff_season <- diffnet(net_season,
                       diffMethod = "discordant", 
                       adjust = "lfdr",
                       cores = 12)
?diffnet
summary(diff_season)
plot(diff_season, 
     cexNodes = 3, 
     cexLegend = 10,
     cexTitle = 4,
     mar = c(2,2,8,5),
     legendGroupnames = c("Metaplasia", "Cancer"),
     legendPos = c(0.7,1.6))
?plot.diffnet
