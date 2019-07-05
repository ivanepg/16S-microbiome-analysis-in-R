#20190502 
Packages <- c("dplyr", "ggplot2", "rstan", "dada2", "DECIPHER", "DESeq2", "ggplot2", "gplots", "gridExtra", "latticeExtra", "permute", "phangorn", "phyloseq", "shiny", "vegan", "structSSI")
lapply(Packages, library, character.only = TRUE)


path<-("/Users/ivanegerasmio/16s_intestine")
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, trimLeft= c(19,20), minLen=100, maxEE=c(2,3), truncQ=2, rm.phix=TRUE, matchIDs=TRUE, compress=TRUE, multithread=TRUE)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]
dadaRs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE, verbose=TRUE)
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "Nonchimera")
rownames(track) <- sample.names
head(track)
write.csv(track, file = "track.csv")

#Assigning taxonomy:
#silva
taxasilva <- assignTaxonomy(seqtab.nochim, "/Users/ivanegerasmio/silva_nr_v132_train_set.fa.gz", multithread=TRUE,tryRC=TRUE)
#gg  
#taxagg <- assignTaxonomy(seqtab.nochim, "/Users/ivanegerasmio/gg_13_8_train_set_97.fa.gz", multithread=TRUE,tryRC=TRUE)
taxa.print <- taxasilva
rownames(taxa.print) <- NULL
head(taxa.print)

#Construct phylogenetic tree:
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign) 
treeNJ <- NJ(dm)
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
#detach("package:phangorn", unload=TRUE)
#try
fitGTR <- optim.pml(fitGTR, rearrangement = "stochastic", ratchet.par = list(iter = 5L, maxit = 5L, prop = 1/3))

#Import mapping file
samdf<-read.csv("/Users/ivanegerasmio/map16s96_intestinev1.csv", header = TRUE)
names(samdf) [1] <- "#SampleID"
samdf
rownames(samdf)
samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out
samdf

#Importing data into phyloseq
ps<-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),sample_data(samdf), tax_table(taxasilva),phy_tree(fitGTR$tree))
ps

pslog <- transform_sample_counts(ps, function(x) log(1 + x))
write.csv(pslog@otu_table, file = "pslog.csv")

out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
#Grp
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "Grp", title="Weighted Unifrac") +
  labs(col = "Grp") +
  coord_fixed(sqrt(evals[2] / evals[1]))
#Tissue
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "Tissue") +
  labs(col = "Tissue") +
  coord_fixed(sqrt(evals[2] / evals[1]))
#Time
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "Time") +
  labs(col = "Time") +
  coord_fixed(sqrt(evals[2] / evals[1]))
#TrmntGrp
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "TrmntGrp") +
  labs(col = "TrmntGrp") +
  coord_fixed(sqrt(evals[2] / evals[1]))

#Bray NMDS
ord.nmds.bray <- ordinate(pslog, method="NMDS", distance="bray")
plot_ordination(pslog, ord.nmds.bray, color="Grp", title="Bray NMDS")

plot_richness(ps, x="Grp", measures=c("Shannon", "Simpson"), color="Grp")
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Grp", fill="Genus") + facet_wrap(~Grp + Time, scales="free_x")

#DESeq2
#use variance stabilizing transformation
ps_dds <- phyloseq_to_deseq2(ps, ~Grp)
varianceStabilizingTransformation(ps_dds, blind = TRUE, fitType = "parametric")
ps_dds <- estimateSizeFactors(ps_dds)
ps_dds <- estimateDispersions(ps_dds)
abund <- getVarianceStabilizedData(ps_dds)

##shorten the names of each microbe
short_names <- substr(rownames(abund), 1, 5)%>%
  make.names(unique = TRUE)
rownames(abund) <- short_names

el <- phy_tree(pslog)$edge
el0 <- el
el0 <- el0[nrow(el):1, ]
el_names <- c(short_names, seq_len(phy_tree(pslog)$Nnode))
el[, 1] <- el_names[el0[, 1]]
el[, 2] <- el_names[as.numeric(el0[, 2])]
unadj_p <- treePValues(el, abund, sample_data(pslog)$Grp)

hfdr_res <- hFDR.adjust(unadj_p, el, .75)
summary(hfdr_res)

abund_sums <- rbind(data.frame(sum = colSums(abund),
                               sample = colnames(abund),
                               type = "DESeq2"),
                    data.frame(sum = rowSums(otu_table(pslog)),
                               sample = rownames(otu_table(pslog)),
                               type = "log(1 + x)"))
ggplot(abund_sums) +
  geom_histogram(aes(x = sum), binwidth = 20) +
  facet_grid(type ~ .) +
  xlab("Total abundance within sample")

plot(hfdr_res, height = 5000) # opens in a browser

options(width=100)
tax <- tax_table(pslog)[, c("Family", "Genus")] %>%
  data.frame()
tax$seq <- short_names
hfdr_res@p.vals$seq <- rownames(hfdr_res@p.vals)
difftax<-tax %>%
  left_join(hfdr_res@p.vals) %>%
  arrange(adjp) %>% head(30)

difftax100<-tax %>%
  left_join(hfdr_res@p.vals) %>%
  arrange(adjp) %>% head(100)


write.csv(ps@otu_table, file = "ps_otutable.csv")
write.csv(abund, file = "abund_shortnames.csv")
write.csv(difftax100, file = "difftax100.csv")
write.csv(difftax, file = "difftax30.csv")

###=====



##Taxonomic filtering
rank_names(ps)
table(tax_table(ps)[, "Phylum"], exclude = NULL)
write.csv(tax_table(ps), file = "taxtable.csv")
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("","uncharacterized"))
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))
#Compute the total and average prevalences of the features in each phylum.
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# Define phyla to filter
#filterPhyla = c("Arthropoda", "FBP","Lentisphaerae","Nitrospirae", "Retaria", "Spirochaetes", "Synergistetes", "Thermotogae", "WS2", "Euglenozoa")
filterPhyla = c("Arthropoda", "Euglenozoa")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps0, !Phylum %in% filterPhyla)
ps1
##Prevalence filtering
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps0)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps0)
ps2
write.csv(ps2@otu_table, file = "ps2int_otutable.csv")

##Agglomerate taxa
# How many genera would be present after filtering?
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
## [1] 187
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
ps3
write.csv(ps3@otu_table, file = "ps3_intestine_otutable.csv")
write.csv(ps3@tax_table, file = "ps3_intestine_taxtable.csv")

h1 = 0.4
ps4 = tip_glom(ps2, h = h1)
ps4
multiPlotTitleTextSize = 8
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
# group plots together
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)

##Abundance value transformation
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Grp",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
# Transform to relative abundance. Save as new object.
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})
plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2, plotBefore, plotAfter)
##Subset by taxonomy
psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)
#Pre-processing
qplot(sample_data(ps)$Grp, geom = "histogram") + xlab("Grp")
qplot(log10(rowSums(otu_table(ps)))) +
  xlab("Logged counts-per-sample")

#=======20190509
#Figures
plot_richness(ps, x="Grp", measures=c("Shannon", "Simpson"), color="Grp")
plot_richness(ps, x="Grp", color="Grp")
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(pslog, ord.nmds.bray, color="Grp", title="Bray NMDS")
plot_ordination(ps.prop, ord.nmds.bray, color="Grp", title="Bray NMDS")
plot_ordination(pslog, ord.nmds.bray, color="Time", title="Bray NMDS")

#subset_samples
summary(sample_data(ps)$Time)
ps_week1 <- subset_samples(ps, Time=="Week 1")
ps_week4 <- subset_samples(ps, Time=="Week 4")
ps_week5 <- subset_samples(ps, Time=="Week 5")

pcoa.wuf.log <- ordinate(ps.prop, method = "PCoA", distance ="wunifrac")
plot_ordination(ps.prop, out.wuf.log, color="Grp", title="PcoA Weighted UniFrac")
plot_ordination(ps.prop, out.wuf.log, color="Time", title="PcoA Weighted UniFrac")
plot_ordination(ps.prop, out.wuf.log, color="TrmntGrp", title="PcoA Weighted UniFrac")

pcoa.uuf.log <- ordinate(ps.prop, method = "PCoA", distance ="uunifrac")
plot_ordination(ps.prop, pcoa.uuf.log, color="Grp", title="PcoA Unweighted UniFrac")


#week1
ps.propw1 <- transform_sample_counts(ps_week1, function(otu) otu/sum(otu))
ord.nmds.brayw1 <- ordinate(ps.propw1, method="NMDS", distance="bray")
plot_ordination(ps.propw1, ord.nmds.brayw1, color="Grp", title="Bray NMDS")
out.wuf.logw1 <- ordinate(ps.propw1, method = "MDS", distance = "wunifrac")
plot_ordination(ps.propw1, out.wuf.logw1, color="Grp", title="Weighted UniFrac")

pcoa.wuf.logw1 <- ordinate(ps.propw1, method = "PCoA", distance ="wunifrac")
plot_ordination(ps.propw1, out.wuf.logw1, color="Grp", title="PcoA Weighted UniFrac")
pcoa.uuf.logw1 <- ordinate(ps.propw1, method = "PCoA", distance ="uunifrac")
plot_ordination(ps.propw1, pcoa.uuf.logw1, color="Grp", title="PcoA Unweighted UniFrac")

#week4
ps.propw4 <- transform_sample_counts(ps_week4, function(otu) otu/sum(otu))
ord.nmds.brayw4 <- ordinate(ps.propw4, method="NMDS", distance="bray")
plot_ordination(ps.propw4, ord.nmds.brayw4, color="Grp", title="Bray NMDS")
out.wuf.logw4 <- ordinate(ps.propw4, method = "MDS", distance = "wunifrac")
plot_ordination(ps.propw4, out.wuf.logw4, color="Grp", title="Weighted UniFrac")

pcoa.wuf.logw4 <- ordinate(ps.propw4, method = "PCoA", distance ="wunifrac")
plot_ordination(ps.propw4, out.wuf.logw4, color="Grp", title="PcoA Weighted UniFrac")
pcoa.uuf.logw4 <- ordinate(ps.propw4, method = "PCoA", distance ="uunifrac")
plot_ordination(ps.propw4, pcoa.uuf.logw4, color="Grp", title="PcoA Unweighted UniFrac")


#week5_stopped1week
ps.propw5 <- transform_sample_counts(ps_week5, function(otu) otu/sum(otu))
ord.nmds.brayw5 <- ordinate(ps.propw5, method="NMDS", distance="bray")
plot_ordination(ps.propw5, ord.nmds.brayw5, color="Grp", title="Bray NMDS")
out.wuf.logw5 <- ordinate(ps.propw5, method = "MDS", distance = "wunifrac")
plot_ordination(ps.propw5, out.wuf.logw5, color="Grp", title="Weighted UniFrac")

pcoa.wuf.logw5 <- ordinate(ps.propw5, method = "PCoA", distance ="wunifrac")
plot_ordination(ps.propw5, out.wuf.logw5, color="Grp", title="PcoA Weighted UniFrac")
pcoa.uuf.logw5 <- ordinate(ps.propw5, method = "PCoA", distance ="uunifrac")
plot_ordination(ps.propw5, pcoa.uuf.logw5, color="Grp", title="PcoA Unweighted UniFrac")


#Richness_Week1,4,5
plot_richness(ps_week1, x="Grp", measures=c("Shannon", "Simpson"), color="Grp")
plot_richness(ps_week4, x="Grp", measures=c("Shannon", "Simpson"), color="Grp")
plot_richness(ps_week5, x="Grp", measures=c("Shannon", "Simpson"), color="Grp")

#Week1
top20w1 <- names(sort(taxa_sums(ps_week1), decreasing=TRUE))[1:20]
ps.top20_w1 <- transform_sample_counts(ps_week1, function(OTU) OTU/sum(OTU))
ps.top20_w1 <- prune_taxa(top20w1, ps.top20_w1)
plot_bar(ps.top20_w1, x="Grp", fill="Genus") + facet_wrap(~Grp + Time, scales="free_x")

top30w1 <- names(sort(taxa_sums(ps_week1), decreasing=TRUE))[1:30]
ps.top30_w1 <- transform_sample_counts(ps_week1, function(OTU) OTU/sum(OTU))
ps.top30_w1 <- prune_taxa(top30w1, ps.top30_w1)
plot_bar(ps.top30_w1, x="Grp", fill="Genus") + facet_wrap(~Grp + Time, scales="free_x")

top15w1 <- names(sort(taxa_sums(ps_week1), decreasing=TRUE))[1:15]
ps.top15_w1 <- transform_sample_counts(ps_week1, function(OTU) OTU/sum(OTU))
ps.top15_w1 <- prune_taxa(top15w1, ps.top15_w1)
plot_bar(ps.top15_w1, x="Grp", fill="Genus") + facet_wrap(~Grp + Time, scales="free_x")



#Week4
top20w4 <- names(sort(taxa_sums(ps_week4), decreasing=TRUE))[1:20]
ps.top20_w4 <- transform_sample_counts(ps_week4, function(OTU) OTU/sum(OTU))
ps.top20_w4 <- prune_taxa(top20w4, ps.top20_w4)
plot_bar(ps.top20_w4, x="Grp", fill="Genus") + facet_wrap(~Grp + Time, scales="free_x")

top15w4 <- names(sort(taxa_sums(ps_week4), decreasing=TRUE))[1:15]
ps.top15_w4 <- transform_sample_counts(ps_week4, function(OTU) OTU/sum(OTU))
ps.top15_w4 <- prune_taxa(top15w4, ps.top15_w4)
plot_bar(ps.top15_w4, x="Grp", fill="Genus") + facet_wrap(~Grp + Time, scales="free_x")


#Week5
top20w5 <- names(sort(taxa_sums(ps_week5), decreasing=TRUE))[1:20]
ps.top20_w5 <- transform_sample_counts(ps_week5, function(OTU) OTU/sum(OTU))
ps.top20_w5 <- prune_taxa(top20w5, ps.top20_w5)
plot_bar(ps.top20_w5, x="Grp", fill="Genus") + facet_wrap(~Grp + Time, scales="free_x")

#Week5
top15w5 <- names(sort(taxa_sums(ps_week5), decreasing=TRUE))[1:15]
ps.top15_w5 <- transform_sample_counts(ps_week5, function(OTU) OTU/sum(OTU))
ps.top15_w5 <- prune_taxa(top15w5, ps.top15_w5)
plot_bar(ps.top15_w5, x="Grp", fill="Genus") + facet_wrap(~Grp + Time, scales="free_x")


plot_richness(ps, x="Grp", measures=c("Shannon", "Simpson"), color="TrmntGrp")
plot_bar(ps.top20, x="Grp", fill="Genus") + facet_wrap(~Grp + Time, scales="free_x")
plot_bar(ps.top20, x="Grp", fill="Genus") + facet_wrap(~TrmntGrp, scales="free_x")
plot_bar(ps.top20, x="TrmntGrp", fill="Genus") + facet_wrap(~TrmntGrp, scales="free_x")

psw1_alpha_div <- estimate_richness(ps_week1, split = TRUE)
psw4_alpha_div <- estimate_richness(ps_week4, split = TRUE)
psw5_alpha_div <- estimate_richness(ps_week5, split = TRUE)

plot_richness(ps_week1, x="Grp", measures=c("Simpson"), color="Grp")
plot_richness(ps_week4, x="Grp", measures=c("Simpson"), color="Grp")
plot_richness(ps_week5, x="Grp", measures=c("Simpson"), color="Grp")


write.csv(track, file = "track.csv")

write.csv(pcoa.wuf.logw1[["vectors"]], file="wuf_vectorsw1.csv")
write.csv(pcoa.wuf.logw4[["vectors"]], file="wuf_vectorsw4.csv")
write.csv(pcoa.wuf.logw5[["vectors"]], file="wuf_vectorsw5.csv")


pslogw1 <- transform_sample_counts(ps_week1, function(x) log(1 + x))
pcoa.wuf.logw1 <- ordinate(pslogw1, method = "PCoA", distance ="wunifrac")
plot_ordination(pslogw1, out.wuf.logw1, color="Grp", title="PcoA Weighted UniFrac")
pcoa.uuf.logw1 <- ordinate(pslogw1, method = "PCoA", distance ="uunifrac")
plot_ordination(pslogw1, pcoa.uuf.logw1, color="Grp", title="PcoA Unweighted UniFrac")

pslogw4 <- transform_sample_counts(ps_week4, function(x) log(1 + x))
pcoa.wuf.logw4 <- ordinate(pslogw4, method = "PCoA", distance ="wunifrac")
plot_ordination(pslogw4, out.wuf.logw4, color="Grp", title="PcoA Weighted UniFrac")
pcoa.uuf.logw4 <- ordinate(pslogw4, method = "PCoA", distance ="uunifrac")
plot_ordination(pslogw4, pcoa.uuf.logw4, color="Grp", title="PcoA Unweighted UniFrac")


pslogw5 <- transform_sample_counts(ps_week5, function(x) log(1 + x))
pcoa.wuf.logw5 <- ordinate(pslogw5, method = "PCoA", distance ="wunifrac")
plot_ordination(pslogw5, out.wuf.logw5, color="Grp", title="PcoA Weighted UniFrac")
pcoa.uuf.logw5 <- ordinate(pslogw5, method = "PCoA", distance ="uunifrac")
plot_ordination(pslogw5, pcoa.uuf.logw5, color="Grp", title="PcoA Unweighted UniFrac")

plot_bar(ps.top20, x="Grp", fill="Genus") + facet_wrap(~Time, scales="free_x")
plot_bar(ps.top20, x="Grp", fill="Phylum") + facet_wrap(~Time, scales="free_x")
#_____


pslog_w1 <- transform_sample_counts(ps_week1, function(x) log(1 + x))
out.wuf.log_w1 <- ordinate(pslog_w1, method = "MDS", distance = "wunifrac")
plot_ordination(pslog_w1, out.wuf.log_w1, color = "Grp", title="MDS Weighted UniFrac")

pcoa.wuf.log_w1 <- ordinate(pslog_w1, method = "PCoA", distance ="wunifrac")
plot_ordination(pslog_w1, out.wuf.log_w1, color="Grp", title="PcoA Weighted UniFrac")

pcoa.uuf.log_w1 <- ordinate(pslog_w1, method = "PCoA", distance ="uunifrac")
plot_ordination(pslog_w1, pcoa.uuf.log_w1, color="Grp", title="PcoA Unweighted UniFrac")



table(tax_table(ps)[, "Phylum"], exclude = NULL)
Phylum<-table(tax_table(ps)[, "Phylum"], exclude = NULL)
write.csv(Phylum, file = "Phylum_intestine.csv")
Class<-table(tax_table(ps)[, "Class"], exclude = NULL)
write.csv(Class, file = "Class_intestine.csv")
Genus<-table(tax_table(ps)[, "Genus"], exclude = NULL)
write.csv(Genus, file = "Genus_intestine.csv")


table(tax_table(ps2)[, "Phylum"], exclude = NULL)
Phylumps2<-table(tax_table(ps2)[, "Phylum"], exclude = NULL)
write.csv(Phylumps2, file = "Phylumps2_intestine.csv")
Classps2<-table(tax_table(ps2)[, "Class"], exclude = NULL)
write.csv(Classps2, file = "Classps2_intestine.csv")
Genusps2<-table(tax_table(ps2)[, "Genus"], exclude = NULL)
write.csv(Genusps2, file = "Genusps2_intestine.csv")


write.csv(ps2@otu_table, file = "ps2_intestineforvenn.csv")


table(tax_table(ps)[, "Phylum"], exclude = NULL)
Phylumps<-table(tax_table(ps)[, "Phylum"], exclude = NULL)
write.csv(Phylumps, file = "Phylumps_intestine.csv")

table(tax_table(ps1)[, "Phylum"], exclude = NULL)
Phylumps1<-table(tax_table(ps1)[, "Phylum"], exclude = NULL)
write.csv(Phylumps1, file = "Phylumps1_intestine.csv")

#----
write.csv(ps@otu_table, file = "psotu.csv")

