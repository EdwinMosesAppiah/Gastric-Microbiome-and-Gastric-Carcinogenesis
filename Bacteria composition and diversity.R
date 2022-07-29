# ##############Alpha diversity##############
library(microbiome)
library(ggpubr)
library(knitr)
library(vegan)
library(tsnemicrobiota)
# load phyloseq data

M.f <- readRDS("/media/edwin/Transcend/datasets/dada2/meta_analysis/filtered_meta_analysis_phyloseq.RDS")

# perform alpha diversity measures
index = alpha(M.f,c("observed", "shannon", "simpson", "chao1", "pielou")) #index=c("observed", "shannon", "simpson", "chao1"))

kable(head(index))

ps1.meta <- meta(M.f)

kable(head(ps1.meta))

ps1.meta$Shannon <- index$diversity_shannon 
ps1.meta$evennesSimpson <- index$evenness_simpson
ps1.meta$dominance <- index$dominance_simpson
ps1.meta$evenpielou <- index$evenness_pielou
ps1.meta$chao1 <- index$chao1

avg <- aggregate(round(ps1.meta$Shannon,2), list(ps1.meta$Conditions), mean, na.rm = TRUE)
avg1 <- aggregate(round(ps1.meta$evennesSimpson,2), list(ps1.meta$Conditions), mean, na.rm = TRUE)
avg2 <- aggregate(round(ps1.meta$dominance,2), list(ps1.meta$Conditions), mean, na.rm = TRUE)
avg3 <- aggregate(round(ps1.meta$dominance,2), list(ps1.meta$Conditions), mean, na.rm = TRUE)

avg[,3] <- avg1$x
avg[,4] <- avg2$x
avg[,5] <- avg3$x

colnames(avg) <- c("Conditions","Shannon","evennesSimpson", "dominance", "dominance" )
avg

t <- kruskal.test(ps1.meta$Shannon~ps1.meta$Conditions)
t1 <- kruskal.test(ps1.meta$evennesSimpson~ps1.meta$Conditions)
t2 <- kruskal.test(ps1.meta$dominance~ps1.meta$Conditions) 
t3 <- kruskal.test(ps1.meta$chao1~ps1.meta$Conditions)


avg[5,2] <- t$p.value
avg[5,3] <- t1$p.value
avg[5,4] <- t2$p.value
avg[5,5] <- t3$p.value


avg[,1]
avg<-data.frame(avg)
kable(avg)
avg$Conditions %>% replace_na("'P-value'")
# create a list of pairwise comaprisons
ps1.meta$Conditions <- as.factor(ps1.meta$Conditions)
gastrointest_disord <- levels(ps1.meta$Conditions) # get the variables

# make a pairwise list that we want to compare.
gastroin <- combn(seq_along(gastrointest_disord), 2, simplify = FALSE, FUN = function(i)gastrointest_disord[i])

print(gastroin)

p1 <- ggboxplot(ps1.meta, x = "Conditions", y = "Shannon",
                color = "Conditions", palette =c("#D14285", "#6DDE88", "#652926", "#FF00FF"),
                add = "jitter", order = c("Healthy_control", "Gastritis", "Metaplasia", "Cancer")) + stat_compare_means(comparisons = gastroin, label = "p.signif")+stat_compare_means(label.y = 10 )  

p2 <- ggboxplot(ps1.meta, x = "Conditions", y = "dominance",
                color = "Conditions", palette =c("#D14285", "#6DDE88", "#652926", "#FF00FF"),
                add = "jitter", order = c("Healthy_control", "Gastritis", "Metaplasia", "Cancer")) + stat_compare_means(comparisons = gastroin, label = "p.signif")+ stat_compare_means(label.y = 2 )  

p3 <- ggboxplot(ps1.meta, x = "Conditions", y = "chao1",
                color = "Conditions", palette =c("#D14285", "#6DDE88", "#652926", "#FF00FF"),
                add = "jitter", order = c("Healthy_control", "Gastritis", "Metaplasia", "Cancer")) + stat_compare_means(comparisons = gastroin, label = "p.signif")+stat_compare_means(label.y = 600 )  

p4 <- ggboxplot(ps1.meta, x = "Conditions", y = "evennesSimpson",
                color = "Conditions", palette =c("#D14285", "#6DDE88", "#652926", "#FF00FF"),
                add = "jitter", order = c("Healthy_control", "Gastritis", "Metaplasia", "Cancer")) + stat_compare_means(comparisons = gastroin, label = "p.signif")+stat_compare_means(label.y = 2 )  

ggarrange(p1, p2, p3, p4 + rremove("x.text"), 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)


# ################################ Bacterial Composition Analysis ##################

#### Phyla prevalence across samples

filterPhyla = c("Aminicenantes", "Aquificae", " Atribacteria", "BRC1", "Caldiserica",
                "Chlorobi", " Cloacimonetes", "Crenarchaeota", " Diapherotrites", 
                "Dictyoglomi", "Elusimicrobia", "Fibrobacteres", " Hydrogenedentes", "Lentisphaerae",
                "Lentisphaerae", "Marinimicrobia", "Microgenomates", " Nanohaloarchaeota",
                " Nitrospinae", "Pacearchaeota", "Poribacteria", "SR1", " Tenericutes", "Thaumarchaeota",
                "Thermotogae", "Woesearchaeota")

ps1 = subset_taxa(M.f, !Phylum %in% filterPhyla)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps1),
               MARGIN = ifelse(taxa_are_rows(ps1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps1),
                    tax_table(ps1))

# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")



filter_abun_Phyla = c("Actinobacteria", "Bacteriodetes", "Firmicutes", "Fusobacteria", "Acidobacteria")
sort(phyloseq::sample_sums(M.f))
M.f1 <- phyloseq::subset_samples(M.f, phyloseq::sample_sums(M.f) > 5000)
M.f1 <- phyloseq::prune_taxa(phyloseq::taxa_sums(M.f1) > 0, M.f1)

#### Phyla changes across various conditions
physeq_rel_abund = phyloseq::transform_sample_counts(M.f1, function(x){x /
    sum(x)})
phyloseq::otu_table(M.f1)[1:15, 1:15]
phyloseq::otu_table(physeq_rel_abund)[1:15, 1:15]

filter_abun_Phyla = c("Actinobacteria","Bacteroidetes", "Firmicutes", "Fusobacteria", "Acidobacteria", "Proteobacteria")
physeq_rel_abund = subset_taxa(physeq_rel_abund, Phylum %in% filter_abun_Phyla)

phyloseq::plot_bar(physeq_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position =
             "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Conditions, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



# ################ Beta diversity ################

#### Permanova test ##
permanova <- adonis(t(otu_table(M.f)) ~ Conditions,
                    data = data.frame(sample_data(ps)), permutations=99, method = "bray")

permanova 

# Anosim test
dune.dist <- vegdist(t(otu_table(M.f)), method = "bray")
dune.ano <- with(sample_data(M.f), anosim(dune.dist, Conditions))
summary(dune.ano)

simpson <- diversity(t(otu_table(M.f)), "simpson") # or assign to var.
simpson[1:5]

# Betadisper
mod <- with(sample_data(M.f), betadisper(dune.dist, Conditions))
anova(mod)
permDisper(mod)
plot(mod)
boxplot(mod)

# TSNE 
tsne_res <- tsne_phyloseq(M.f, distance='wunifrac',
                          perplexity = 8, verbose=0, rng_seed = 3901)
plot_tsne_phyloseq(M.f, tsne_res, axes = 1:2,
                   color = 'Conditions', title='t-SNE (Weighted UniFrac)')


