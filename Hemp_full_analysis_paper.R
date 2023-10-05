
# ------------- FOLLOWING TO BE RUN IN COMMAND LINE -------------

#####using Cutadapt v2.6
##starting with demultiplexed compressed fastq files
##create file 'samples' with sample names based on file names
##remove primers from forward reads (R1) and reverse reads (R2)
##command line code as follows:
#
# module load anaconda/3
# source activate cutadaptenv
#
#
# ls *R1*.gz | cut -f-7 -d "_" > samples
#
# 
#
# for sample in $(cat $samples)
# 
#   do
# 
#   echo "On sample: $sample"
#   cutadapt -g ^CCTACGGGNGGCWGCAG \
#   -m 220 -M 240 --discard-untrimmed \
#   -o ${sample}_R1_trimmed.fq \
#   ${sample}_R1.fastq.gz \
#   > cutadapt_stats_${sample}_R1.txt 2>&1 || exit 1
# 
#   done
# 
# done
# 
# 
# for sample in $(cat samplename)
# 
#   do
# 
#   echo "On sample: $sample"
#   cutadapt -g ^CCTACGGGNGGCWGCAG \
#   -m 220 -M 240 --discard-untrimmed \
#   -o ${sample}_R2_trimmed.fq \
#   ${sample}_R2.fastq.gz \
#   > cutadapt_stats_${sample}_R2.txt 2>&1 || exit 1
# 
#   done
# 
# done





#------------------------------------------------------------- FOLLOWING TO BE RUN IN R -------------------------------------------------------------


#####################################################################################
#####step 1 - denoise and produce feature table 
library('ggplot2')
library('dada2')
library('phyloseq')

#read sample names (file created in command line above)
samples <- scan("samples", what="character") 

#store names of fastq files
forward_reads <- paste0(samples, "_R1_trimmed.fq.gz")
reverse_reads <- paste0(samples, "_R2_trimmed.fq.gz")

studyname <- paste0("Hemp")

#to store names of filtered fastq files
filtered_forward_reads <- paste0(samples, "_R1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_R2_filtered.fq.gz")

#create function to plot read quality scores
plotQualityProfile_custom <- function(read) {
  plotQualityProfile(read) + scale_x_continuous(breaks = seq(from = 0, to = 300, by = 10)) + theme(axis.text.x = element_text(size = 5))
}
#produce plots with custom function
pdf_name <- paste0(studyname, "_quality_profile_forward_reads.pdf")
pdf(pdf_name)
lapply(forward_reads, plotQualityProfile_custom)
dev.off()
pdf_name <- paste0(studyname, "_quality_profile_reverse_reads.pdf")
pdf(pdf_name)
lapply(reverse_reads, plotQualityProfile_custom)
dev.off

##trim reads down based on a quality threshold based on quality plots
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                              rm.phix=TRUE, truncLen=c(233, 229), compress=TRUE, matchIDs = TRUE)
head(filtered_out) 

#check new plots
pdf_name <- paste0(studyname, "_filtered_quality_profile_forward_reads.pdf")
pdf(pdf_name)
lapply(filtered_forward_reads, plotQualityProfile_custom)
dev.off()
pdf_name <- paste0(studyname, "_filtered_quality_profile_reverse_reads.pdf")
pdf(pdf_name)
lapply(filtered_reverse_reads, plotQualityProfile_custom)
dev.off()

###generating an error model of our data
err_forward_reads <- learnErrors(filtered_forward_reads, multithread = TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread = TRUE)

plotErrors_PDFname <- paste0(studyname, "_plotErrors.pdf")
pdf(plotErrors_PDFname, paper = "a4")
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)
dev.off()
#explanation of plots from https://benjjneb.github.io/dada2/tutorial.html#learn-the-error-rates

###inferring ASVs
#pseudo-pooling used here:
dada_forward <- dada(filtered_forward_reads, err=err_forward_reads, pool="pseudo")
dada_reverse <- dada(filtered_reverse_reads, err=err_reverse_reads, pool="pseudo")

###merging forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, filtered_forward_reads, dada_reverse,
                               filtered_reverse_reads, trimOverhang=TRUE, minOverlap=20)
head(merged_amplicons)

###generating a count table
seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

###chimera identification
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", verbose=TRUE)
dim(seqtab.nochim) 
sum(seqtab.nochim)/sum(seqtab)

#check quantity of reads filtered at each step
getN <- function(x) sum(getUniques(x))
track <- cbind(filtered_out, sapply(dada_forward, getN), sapply(dada_reverse, getN), sapply(merged_amplicons, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- samples
head(track)

#save feature table
saveRDS(seqtab.nochim, paste0(studyname, "_feature_table.rds"))


#####################################################################################
#####step 2 - assign taxa and produce phylogenetic tree

library('ggplot2')
library('dada2')
library('phyloseq')
library('dplyr')
library('DECIPHER')
library('phangorn')

###load objects
ft<-readRDS("Hemp_feature_table.rds")
samples <- scan("samples", what="character")

###assigning taxonomy
#using the DECIPHER package
library(DECIPHER)
taxa <- assignTaxonomy(ft, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
#add species
taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")
#save object
saveRDS(taxa,"Hemp_taxa_table.rds")


#check output
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)

save.image("Hemp_step2_taxaAssigned.rdata")

#generate phylogenetic tree
ft<-readRDS("Hemp_feature_table.rds")
seqs<-getSequences(ft)
names(seqs) <- seqs 
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) 
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
saveRDS(fitGTR, "Hemp_phylo_tree.rds")

#####################################################################################
#####step 3 - pass to phyloseq

library('ggplot2')
library('dada2')
library('phyloseq')
library('dplyr')
library('stringr')

studyname<-paste0("Hemp")

ft<-readRDS("Hemp_feature_table.rds")

samples.out <- rownames(ft)
samp <- sapply(strsplit(samples.out, "_"), `[`, 2)

earl.ref <- sapply(strsplit(samples.out, "_"), `[`, 1); earl.ref

genotype <- sapply(strsplit(samp, "S"), `[`, 1)

sample_number <- sapply(strsplit(samp, "S"), `[`, 2)

genz.df <- data.frame(Samp = samp, Genotype = genotype, Sample_number = sample_number, Earlham_ref = earl.ref)

#Label blanks
genz.df$Sample_or_control <- "True Sample"
genz.df$Sample_or_control[genz.df$Samp=="BK1"] <- "Control"

genz.df$Sample_or_control <- "True Sample"
genz.df$Sample_or_control[grep("^BK", genz.df$Samp)] <- "Control" 

#import codes to assign genotype names
codes.file <- read.csv("Sample_codes.csv")
codes.file[complete.cases(codes.file), ]
class(codes.file) #is dataframe

#combine codes.file and genz.df by common column
genz.df$Genotype_name<-codes.file[match(paste(genz.df$Genotype),paste(codes.file$Code)),"name"]
genz.df$Genotype_name[which(genz.df$Genotype_name=="Markant")]<-"Markant 3"
genz.df

genz.df$Sample_number[is.na(genz.df$Sample_number)] <- "Blank"
genz.df$Genotype_name[is.na(genz.df$Genotype_name)] <- "Blank"

genz.df$Sample_number[is.na(genz.df$Sample_number)] <- "Blank"
genz.df$Genotype_name[is.na(genz.df$Genotype_name)] <- "Blank"

#add additional genotype info
pheno.data<-read.csv("Hemp_metadata.csv", na.strings=c("","NA"))
colnames(pheno.data)[1]<-"name"
colnames(pheno.data)
#Origin
genz.df$Origin<-pheno.data[match(paste(genz.df$Genotype_name),paste(pheno.data$name)),"Origin"]
genz.df$Genotype_name[is.na(genz.df$Origin)]
#Derived from
genz.df$Derived.from<-pheno.data[match(paste(genz.df$Genotype_name),paste(pheno.data$name)),"Derived.from"]
genz.df$Genotype_name[is.na(genz.df$Derived.from)]
#Supplier
genz.df$Supplier<-pheno.data[match(paste(genz.df$Genotype_name),paste(pheno.data$name)),"Supplier"]
genz.df$Genotype_name[is.na(genz.df$Supplier)]
#Flowering.type
genz.df$Flowering.type<-pheno.data[match(paste(genz.df$Genotype_name),paste(pheno.data$name)),"Flower.type"]
genz.df$Genotype_name[is.na(genz.df$Flowering.type)]
#Maturity
genz.df$Maturity<-pheno.data[match(paste(genz.df$Genotype_name),paste(pheno.data$name)),"Maturity.short.form"]
genz.df$Genotype_name[is.na(genz.df$Maturity)]
#Country.of.origin
genz.df$Country.of.origin<-pheno.data[match(paste(genz.df$Genotype_name),paste(pheno.data$name)),"Country.of.origin.lit.search"]
genz.df$Genotype_name[is.na(genz.df$Country.of.origin)]
#Supplier.country
genz.df$Supplier.country<-pheno.data[match(paste(genz.df$Genotype_name),paste(pheno.data$name)),"Supplier.country"]
genz.df$Genotype_name[is.na(genz.df$Supplier.country)]

rownames(genz.df) <- samples.out #this doesn't work with just one sample
colnames(genz.df)

#creating phyloseq object from DADA2 outputs
taxa<-readRDS("Hemp_taxa_table.rds")
fitGTR<-readRDS("Hemp_phylo_tree.rds")
ps <- phyloseq(otu_table(ft, taxa_are_rows=FALSE), 
               sample_data(genz.df), 
               tax_table(taxa), 
               phy_tree(fitGTR$tree))
#store DNA sequences in the refseq slot of the phyloseq object
  #and then rename our taxa to ASV1, ASV2 etc
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


#check changes in read counts after chloroplast/mitochondria removals
#remove samples that are removed at a later stage - Futura 75 and Santhica 70
used_samples<-rownames(sample_data(ps)[-which(sample_data(ps)$Genotype_name == "Futura 75" | sample_data(ps)$Genotype_name == "Santhica 70" | sample_data(ps)$Sample_or_control == "Control"),])
ps3<-prune_samples(used_samples, ps)
#remove from ps3
ps.no.c<-subset_taxa(ps3, (Order!="Chloroplast") | is.na(Order))
ps.no.c.m<-subset_taxa(ps.no.c, (Family!="Mitochondria") | is.na(Family))
ps.rc<-as.data.frame(sample_sums(ps3)); colnames(ps.rc)<-"read.count"
ps.rc$read.count
#plot
options(scipen=10000)
g1<-ggplot(data = ps.rc, aes(y=read.count)) +
  geom_boxplot() +
  scale_y_log10(limits = c(50, 2250000)) +
  ggtitle("All amplicons"); g1
sample_sums(ps.no.c); ps.no.c.rc<-as.data.frame(sample_sums(ps.no.c)); colnames(ps.no.c.rc)<-"read.count"
g2<-ggplot(data = ps.no.c.rc, aes(y=read.count)) +
  geom_boxplot() +
  scale_y_log10(limits = c(50, 2250000)) +
  ggtitle("Chloroplasts removed"); g2
ps.no.c.m.rc<-as.data.frame(sample_sums(ps.no.c.m)); colnames(ps.no.c.m.rc)<-"read.count"
g3<-ggplot(data = ps.no.c.m.rc, aes(y=read.count)) +
  geom_boxplot() +
  scale_y_log10(limits = c(50, 2250000)) +
  ggtitle("Chloroplasts and mitochondria removed"); g3
library(cowplot)
plot_grid(g1, g2, g3, ncol = 3)
#ggsave("Read_count_after_organelle_removals.png", plot = last_plot(), device = "png", units = "cm", width = 32, height = 12)
#stats
mean(ps.rc$read.count)
mean(ps.no.c.rc$read.count)
mean(ps.no.c.m.rc$read.count)
mean(ps.no.c.m.rc$read.count)
mean((ps.rc$read.count - ps.no.c.rc$read.count)/ps.rc$read.count) #mean reads that are chloroplasts
mean((ps.no.c.rc$read.count - ps.no.c.m.rc$read.count)/ps.rc$read.count) #mean reads that are mitochondria

#remove chloroplasts before saving the phyloseq object
ps.no.c<-subset_taxa(ps, (Order!="Chloroplast") | is.na(Order))
#and remove mitochondria
ps.no.c.m<-subset_taxa(ps.no.c, (Family!="Mitochondria") | is.na(Family))

#continue with organelles removed
ps<-ps.no.c.m

#save phyloseq object
saveRDS(ps, paste0("phyloseq_object_", studyname, ".rds"))

#####################################################################################
#####step 4 - remove contaminants

library(decontam)
library(phyloseq)
library(ggplot2)
library(microViz)

ps<-readRDS("phyloseq_object_Hemp.rds")

###identify contaminants using decontam package's prevalence method
#visualise library sizes - group controls seperately
df <- as.data.frame(sample_data(ps)) 
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_control)) + geom_point()
#identify contaminant based on prevalence in samples and prevalence in negative controls
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_control == "Control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant) #number of contaminants
head(which(contamdf.prev$contaminant))
#plot prevalence in samples vs prevalence in controls
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_control == "True Sample", ps.pa)
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#generate phyloseq object of contaminants only
contams<-rownames(contamdf.prev[which(contamdf.prev$contaminant),])
cps <- prune_taxa(contams, ps)

#generate phyloseq object excluding contaminants for further analysis
notcontams<-rownames(contamdf.prev[-which(contamdf.prev$contaminant),])
ps <- prune_taxa(notcontams, ps)
#remove controls 
ps <- ps %>% ps_filter(Sample_or_control == "True Sample")
#check (and remove if needed) taxa with NA at phylum level
ps <- subset_taxa(ps, Phylum != "NA")
saveRDS(ps, "Hemp_ps_processed_step4_decontam_prevalence.rds")

#obtain list of asvs in >0 controls but not removed as contaminants
df.pa
npa<-df.pa[-which(df.pa$contaminant),] #remove contaminants
npa<-npa[which(npa$pa.neg>0 & npa$pa.pos>0),] #remove ASVs present only in controls
npa.asvs<-rownames(npa)
nps<-prune_taxa(npa.asvs, ps) #keep only non-contaminant that are present in >0 control
Genos_toRemove<-c("Futura 75", "Santhica 70")
Samples_toRemove<-rownames(sample_data(nps)[which(sample_data(nps)$Genotype_name %in% Genos_toRemove)],)
nps<-prune_samples(!sample_names(nps) %in% Samples_toRemove, nps)
nps #number of taxa is number of ASVs
nps.glom<-tax_glom(nps, taxrank = "Genus", NArm = FALSE) 
df.tt<-as.data.frame(tax_table(nps.glom))
sort(unique(df.tt$Genus)) #genera
df.tt[which(is.na(df.tt$Genus)),] #which NA at genus level
npstt<-as.data.frame(tax_table(nps))
npsrs<-as.data.frame(refseq(nps))
npsdf<-merge(npstt, npsrs, by = "row.names")
colnames(npsdf)<-c("ASV","Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",  "Sequence" )
library(openxlsx)
write.xlsx(npsdf, "NonContaminants_butInControls_table.xlsx")


#look at contaminants
cps #number of taxa is number of ASVs
cps.glom<-tax_glom(cps, taxrank = "Genus", NArm = FALSE) 
df.ctt<-as.data.frame(tax_table(cps.glom))
sort(unique(df.ctt$Genus)) #genera
df.ctt[which(is.na(df.ctt$Genus)),] #which NA at genus level
cpstt<-as.data.frame(tax_table(cps))
cpsrs<-as.data.frame(refseq(cps))
cpsdf<-merge(cpstt, cpsrs, by = "row.names")
colnames(cpsdf)<-c("ASV","Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",  "Sequence" )
library(openxlsx)
write.xlsx(cpsdf, "Contaminants_table.xlsx")



#####################################################################################
#####step 5 - filter samples by read count and filter ASVs by prevelance

library(phyloseq)
library(ggplot2)
library(dplyr)
library(microViz)
library(openxlsx)

ps<-readRDS("Hemp_ps_processed_step4_decontam_prevalence.rds")
studyname <- paste0("Hemp")

###analyse read counts
ps.rc<-as.data.frame(sample_sums(ps)); colnames(ps.rc)<-"read.count"
ps.rc$Samp <- rownames(ps.rc)
rownames(ps.rc)<-NULL
ps.rc$Sample <- sapply(strsplit(ps.rc$Samp, "_"), `[`, 2)
ps.rc$Genotype <- sapply(strsplit(ps.rc$Sample, "S"), `[`, 1)
ps.rc$Sample_number <- sapply(strsplit(ps.rc$Sample, "S"), `[`, 2)

ggplot(data = ps.rc, aes(x=reorder(Sample, -read.count), y=read.count)) + 
  geom_bar(stat="identity") +
  ggtitle("Hemp - Read counts per sample") +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "")
#by genotype
ggplot(data = ps.rc, aes(x=Sample_number, y=read.count)) + 
  geom_bar(stat="identity") +
  facet_grid(~Genotype) +
  ggtitle("Hemp - Read counts per sample") +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "")
#ggsave(filename = "Hemp_readcounts_by_genotype.png", plot = last_plot(), device = "png", units = "cm", width = 100, height = 20)

###H06 and H13 have much lower read counts across all replicates than other genotypes
#remove these to improve data quality
ps<-ps %>% ps_filter(!Genotype == "H06") %>% ps_filter(!Genotype == "H13")

###filter by prevalence
##remove taxa appearing in <3 samples 
ps
ps <- filter_taxa(ps, function (x) {sum(x > 0) > 2}, prune=TRUE)

###obtain stats on total, median and mean read counts
ps.rc2<-as.data.frame(sample_sums(ps)); colnames(ps.rc2)<-"read.count"
ps.rc2$Samp <- rownames(ps.rc2)
rownames(ps.rc2)<-NULL
ps.rc2$Sample <- sapply(strsplit(ps.rc2$Samp, "_"), `[`, 2)
ps.rc2$Genotype <- sapply(strsplit(ps.rc2$Sample, "S"), `[`, 1)
ps.rc2$Sample_number <- sapply(strsplit(ps.rc2$Sample, "S"), `[`, 2)
nrow(ps.rc2)
sum(ps.rc2$read.count)
median(ps.rc2$read.count)
mean(ps.rc2$read.count) 

###save phyloseq object that will be used for all subsequent analysis
saveRDS(ps, "Hemp_ps_processed_step5_readyForAnalysis.rds")

#produce list of ASVs
final_pstt<-as.data.frame(tax_table(ps))
nrow(final_pstt)
final_psrs<-as.data.frame(refseq(ps))
nrow(final_psrs)
rownames(final_pstt)
rownames(final_psrs)
final_psdf<-merge(final_pstt, final_psrs, by = "row.names")
nrow(final_psdf)
head(final_psdf)
colnames(final_psdf)<-c("ASV","Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",  "Sequence" )
nrow(final_psdf)
##add negative control info
ncic<-read.xlsx("NonContaminants_butInControls_table.xlsx")
ncic$ASV
#find ASVs that were in negative control
final_psdf[which(final_psdf$ASV %in% ncic$ASV),]$ASV #fewer than reported because some probably since filtered out by prevalence etc
#check it is the same with pstt
#final_pstt2<-final_pstt; final_pstt2$ASV<-rownames(final_pstt2); rownames(final_pstt2)<-NULL; sort(final_pstt2[which(final_pstt2$ASV %in% ncic$ASV),]$ASV)
#add column about presence in negative control
final_psdf$In.negative.control<-"No"
final_psdf[which(final_psdf$ASV %in% ncic$ASV),]$In.negative.control<-"Yes"
#check number of ASVs
nrow(final_psdf)
#reorder to sensible ASV order
final_psdf$ASV_number <- as.numeric(gsub(".*?([0-9]+).*", "\\1", final_psdf$ASV)) ; final_psdf <- final_psdf[order(final_psdf$ASV_number),]; final_psdf$ASV_number <- NULL 
#remove sequences
final_psdf$Sequence<-NULL
#add contaminants
final_psdf$Contaminant<-"No"
con<-read.xlsx("Contaminants_table.xlsx")
con$Sequence<-NULL
con$In.negative.control<-""
con$Contaminant<-"Yes"
colnames(final_psdf)
colnames(con)
cf<-rbind(final_psdf, con)
cf$In.negative.control
colnames(cf)
cf<-cf[,c(1:8,10,9)]
table(cf$Contaminant)
cf[which(cf$Contaminant=="Yes"),]$In.negative.control
names(cf)[names(cf) == 'In.negative.control'] <- 'Present in a negative control'
write.xlsx(cf, "ASV_final_dataset_table.xlsx")


#####################################################################################
#####step 6 - main data analysis
############################################################
####step 6.1 - summary stats on phyla and genera
library(microbiome)  
library(phyloseq)
library(ggplot2)
library(dplyr)
options(scipen=999)

##phyla
ps<-readRDS("Hemp_ps_processed_step5_readyForAnalysis.rds")
ps.glom.phyla<-tax_glom(ps, taxrank = "Phylum", NArm = FALSE) #no removal of ASVs not classifed to phylum; follows Wallen 2021
#mean abundances
ps.t.p <- transform_sample_counts(ps.glom.phyla, function(OTU) OTU/sum(OTU))
mps.t.p<-psmelt(ps.t.p)
mps.t.p<-mps.t.p[c(19, 3)]
#means
mps.t.p.mean<-aggregate(mps.t.p['Abundance'], by=mps.t.p['Phylum'], mean)
mps.t.p.mean<-mps.t.p.mean[order(-mps.t.p.mean$Abundance),]
mps.t.p.mean
#samples where each phylum is most abundant phylum
ps.glom.phyla<-tax_glom(ps, taxrank = "Phylum", NArm = FALSE) #no removal of ASVs not classifed to Phylum; follows Wallen 2021
mps.t.p<-psmelt(ps.glom.phyla)
tps<-as_tibble(mps.t.p) %>% 
  dplyr::rename(ASV=OTU) %>% 
  dplyr::select(Samp, Abundance, Genotype_name, Phylum) %>% 
  dplyr::group_by(Samp) %>% 
  dplyr::slice_max(Abundance, n=1) %>% 
  dplyr::ungroup() %>%
  dplyr::arrange(Phylum, Samp, Abundance)
dftps<-as.data.frame(tps)
print(tps[order(tps$Samp),], n=48)
unique(tps$Phylum) #phyla that dominate at least one sample
table(tps$Phylum)/48 #% samples were each phylum is most abundant

##genera
ps.rel <- microbiome::transform(ps, "compositional")
ps.glom<-tax_glom(ps.rel, taxrank = "Genus", NArm = FALSE) #no removal of ASVs not classifed to genus; follows Wallen 2021
mps<-psmelt(ps.glom)
#mean proportion of reads per genus
mps.means<- mps %>%
  group_by(Genus) %>%
  dplyr::summarize(Mean = mean(Abundance, na.rm=FALSE))
sorted.means<-mps.means[order(mps.means$Mean, decreasing = TRUE),]
print(sorted.means, n=15)
#mean proportion of reads collectively represented by top 3, top 5, top 10
sum(sorted.means[1:3,]$Mean)
sum(sorted.means[1:5,]$Mean)
sum(sorted.means[1:10,]$Mean)
#plot mean proportion of reads per sample cumulatively represented by subsets of genera of increasing size. 
  #Genera are added to the subset in descending order of mean relative abundance.
colnames(mps)
p<-mps[,c(4,22,23,1,3)]
p.mm<-p %>%
  group_by(Family, Genus, OTU)%>%
  summarize(Abundance.median = median(Abundance), Abundance.mean = mean(Abundance))
p.mm<-arrange(p.mm, desc(Abundance.mean))
out<-c()
for (i in 1:88){
  print(i)
  top<-p.mm$OTU[1:i]
  p2<-(filter(p, OTU %in% top))
  ab<-sum(p2$Abundance)/48
  print(ab)
  out<-append(out, ab)
}
out
df<-data.frame(1:88, out)
colnames(df)<-c("No_taxa", "Mean_proportion_of_community_represented")
ggplot(df, aes(x=No_taxa, y=Mean_proportion_of_community_represented)) + 
  theme_bw() +
  geom_line() +
  geom_point() +
  xlab("Number of genera") +
  ylab("Mean proportion of community\nrepresented by subset of genera (%)") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0.2,1), breaks = seq(0.2, 1, by = 0.2)) 
#ggsave("Hemp_analysis_representation_by_number_of_top_taxa_by_mean_abundance.png", plot = last_plot(), device = "png", units = "cm", width = 15, height = 10)
##stats on top 3 genera
top3<-p.mm$OTU[1:3]
g.ps.core<-prune_taxa(top3, ps.glom)
g.core.prop<-as.data.frame(sample_sums(g.ps.core))
colnames(g.core.prop) <- "Proportion.represented"
g.core.prop$Sample <- rownames(g.core.prop)
g.core.prop$Samp <- sapply(strsplit(g.core.prop$Sample, "_"), `[`, 2)
g.core.prop$Genotype <- sapply(strsplit(g.core.prop$Samp, "S"), `[`, 1)
g.core.prop$Sample<-NULL; rownames(g.core.prop) <- NULL
g.core.prop[order(g.core.prop$Proportion.represented),]
mean(g.core.prop$Proportion.represented) #mean
median(g.core.prop$Proportion.represented) #median
min(g.core.prop$Proportion.represented) #min
max(g.core.prop$Proportion.represented) #max
#by cultivar
gcpc<-g.core.prop %>%
  group_by(Genotype)%>%
  summarize(Prop.mean = mean(Proportion.represented))
arrange(gcpc, desc(Prop.mean))
#samples where top 3 make up >50% of community
gcpBIN<-g.core.prop%>%mutate(Prop.bins = cut(Proportion.represented, breaks = c(-Inf,0.501,Inf)))
table(gcpBIN$Prop.bins) 
##look at how reads assigned to each common genus are distributed across ASVs
flavoASVs<-sorted.ASVmeans[which(sorted.ASVmeans$Genus=="Flavobacterium"),]; flavoASVs
flavoASVs$Mean[1]/sum(flavoASVs$Mean) #proporton of all Flavobacterium reads represented by most abundant ASV
pseudoASVs<-sorted.ASVmeans[which(sorted.ASVmeans$Genus=="Pseudomonas"),]; pseudoASVs
pseudoASVs$Mean[1]/sum(pseudoASVs$Mean) 
pantASVs<-sorted.ASVmeans[which(sorted.ASVmeans$Genus=="Pantoea"),]; pantASVs
pantASVs$Mean[1]/sum(pantASVs$Mean) 
brevASVs<-sorted.ASVmeans[which(sorted.ASVmeans$Genus=="Brevibacillus"),]; brevASVs
brevASVs$Mean[1]/sum(brevASVs$Mean) 
acidASVs<-sorted.ASVmeans[which(sorted.ASVmeans$Genus=="Acidovorax"),]; acidASVs
acidASVs$Mean[1]/sum(acidASVs$Mean) 
hermASVs<-sorted.ASVmeans[which(sorted.ASVmeans$Genus=="Herminiimonas"),]; hermASVs
hermASVs$Mean[1]/sum(hermASVs$Mean) 
massASVs<-sorted.ASVmeans[which(sorted.ASVmeans$Genus=="Massilia"),]; massASVs
massASVs$Mean[1]/sum(massASVs$Mean) 
kasASVs<-sorted.ASVmeans[which(sorted.ASVmeans$Genus=="Kosakonia"),]; kasASVs
kasASVs$Mean[1]/sum(kasASVs$Mean) 
chrysASVs<-sorted.ASVmeans[which(sorted.ASVmeans$Genus=="Chryseobacterium"),]; chrysASVs
chrysASVs$Mean[1]/sum(chrysASVs$Mean) 
bacASVs<-sorted.ASVmeans[which(sorted.ASVmeans$Genus=="Bacillus"),]; bacASVs
bacASVs$Mean[1]/sum(bacASVs$Mean) 
##look at genera that are specific to one, two or three cultivars (and present in all replicates of a cultivar)
ps.rel <- microbiome::transform(ps, "compositional")
ps.glom<-tax_glom(ps.rel, taxrank = "Genus", NArm = FALSE) 
mmps<-psmelt(ps.glom)
mps<-mmps[,c(4, 9, 23, 1, 3)]
glist<-unique(mps[-which(is.na(mps$Genus)),]$Genus) #look only at classified taxa at genus level, no NAs
#produce output showing, for each genus, the samples where it appears at >0.1% abundance 
gens<-c()
samps.out<-c()
for (i in glist){
  print(i)
  gens<-append(gens, i)
  gs<-mps[which(mps$Genus == i & mps$Abundance > 0.001),]
  print(gs)
  samps<-nrow(mps[which(mps$Genus == i & mps$Abundance > 0.001),])
  print(samps)
  samps.out<-append(samps.out, samps)
}
#for genera found to be specific to one, two or three cultivars, how abundant are they where they are present?
#list of genera specific to 3 or fewer cultivars
cslist<-c("Bacillus","Atlantibacter","Lachnoclostridium","Siccibacter","Sediminihabitans","SN8","Stenotrophomonas","Neorhizobium","Variovorax","Cellulomonas","Raoultella","Methylobacterium-Methylorubrum","Simplicispira","Serratia","Acinetobacter","Brevibacillus","Enterococcus","Kluyvera","Paenibacillus","Sanguibacter-Flavimobilis","Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","Erwinia")
length(cslist)
gens<-c()
csout<-c()
for (i in cslist){
  gs<-mps[which(mps$Genus == i & mps$Abundance > 0.001),]
  samps<-nrow(mps[which(mps$Genus == i & mps$Abundance > 0.001),])
  print(paste0("Where present, mean abundance of ",i," is: ",mean(gs$Abundance)))
  cs_out<-mean(gs$Abundance)
  csout<-append(csout, cs_out)
  gens<-append(gens, i)
}
df<-data.frame(gens, csout); colnames(df) <- c("Genus", "Abudance_where_present")
df<-df[order(df$Abudance_where_present, decreasing = TRUE),]
median(df$Abudance_where_present)


############################################################
####step 6.2 - taxa barplots
library(phyloseq)
library(ggplot2)
library(cowplot)
library(microViz)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(ggh4x)

ps<-readRDS("Hemp_ps_processed_step5_readyForAnalysis.rds")

###plot phyla

###plot phyla
#transform
ps.t.p <- transform_sample_counts(ps.glom.phyla, function(OTU) OTU/sum(OTU))
#set colours
PhPalette<-c()
PhPalette["Actinobacteriota"] <- "#D55E00"
PhPalette["Bacteroidota"] <- "#CC79A7"
PhPalette["Cyanobacteria"] <- "#F0E442"
PhPalette["Firmicutes"] <- "#56B4E9"
PhPalette["Proteobacteria"] <- "#009E73"
legend.order<-sort(names(PhPalette))
#prepare data frame
mps.t.p<-psmelt(ps.t.p)
mps.t.p[which(mps.t.p$Genotype_name == "CS Carmagnola"),]$Genotype_name <- "CS"
mps.t.p[which(mps.t.p$Genotype_name == "Eletta Campana"),]$Genotype_name <- "Eletta C."
mps.t.p[which(grepl("Bial", mps.t.p$Genotype_name, fixed = TRUE)),]$Genotype_name <- "Bialobrzeskie"
#plot
hp2<-ggplot(mps.t.p, aes(fill=Phylum, x=Sample_number, y=Abundance)) +
  theme_minimal() +
  theme(axis.ticks.x = element_line(linewidth = 0.1), 
        axis.ticks.y = element_line(linewidth = 0,1),
        axis.ticks.length = unit(2, "pt")) +
  theme(strip.background = element_rect(fill = "lightgrey", color = "lightgrey"),
        strip.text.x = element_text(size = 12)) +
  facet_nested_wrap(~Supplier + Genotype_name, nrow = 2) +
  geom_bar(position = "fill", stat="identity") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme(legend.text = element_text(size = 14)) + 
  theme(legend.title = element_text(size = 14, face = "bold")) + 
  theme(legend.position = "bottom", legend.box.spacing = margin(-10)) +
  theme(legend.key.size=unit(14,"point")) +
  scale_fill_manual(values= PhPalette, breaks=legend.order) +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") +
  ylab("Relative abundance (%)") +
  theme(axis.title = element_text(size=14, face = "bold")) +
  theme(axis.text = element_text(size=12)) +  
  guides(fill=guide_legend(ncol=2)); hp2
#ggsave("Hemp_Phylum_barplot.png", plot = hp2, device = png, height = 24, width = 24, units = "cm", limitsize = TRUE, bg = "white")

###plot top 10 most abundant genera
ps.t <- transform_sample_counts(ps.glom, function(OTU) OTU/sum(OTU))
#find top 10
top10 <- names(sort(taxa_sums(ps.t), decreasing=TRUE))[1:10]
#one is NA at genus level - discount for this plot as it may represent several genera within the family
dfpst<-as.data.frame(tax_table(ps.t)); dfpst2<-(filter(dfpst, row.names(dfpst) %in% top10))
top10 <- names(sort(taxa_sums(ps.t), decreasing=TRUE))[1:12]
top10 <- top10[! top10 %in% c('ASV28', 'ASV91')]
dfpst<-as.data.frame(tax_table(ps.t)); dfpst2<-(filter(dfpst, row.names(dfpst) %in% top10)); dfpst2 #no more NAs
top10_ps <- prune_taxa(top10, ps.t)
top10 <- names(sort(taxa_sums(top10_ps), decreasing=TRUE))
mps.t<-psmelt(ps.t)
mps.t[!mps.t$OTU %in% top10, ]$Genus<-"Other"
unique(mps.t[mps.t$OTU %in% top10, ]$Genus)
mps.t[which(is.na(mps.t$Genus)),]$Genus
getPalette = colorRampPalette(brewer.pal(11, "Set1"))
#to avoid most abundant genera being too similar in colour, randomise order
g.df<-as.data.frame(unique(mps.t$Genus))
set.seed(100)
g.df<-g.df[sample(nrow(g.df)),]
g.df
#manually change some colours for aesthetics
GenusList<-c("Massilia","Herminiimonas", "Kosakonia", "Acidovorax", "Bacillus","Pseudomonas", "Brevibacillus", "Pantoea", "Chryseobacterium", "Flavobacterium", "Other")
GenusPalette = getPalette(length(GenusList))
names(GenusPalette) = GenusList
GenusPalette["NA"] <- "grey55"
GenusPalette["Other"] <- "grey20"
GenusPalette["Bacillus"] <- "olivedrab3"
GenusPalette["Pantoea"] <- "tomato"
GenusPalette["Flavobacterium"] = "skyblue3" 
GenusPalette["Brevibacillus"] = "gold3" 
GenusPalette["Chryseobacterium"] = "darkred" 
GenusPalette["Massilia"] = "pink4" 
GenusPalette["Pseudomonas"] = "orange" 
GenusPalette["Kosakonia"] = "plum" 
GenusPalette
legend.order<-c(sort(GenusList[1:9]), GenusList[10:11])
#plot with nested facets
mps.t[which(mps.t$Genotype_name == "CS Carmagnola"),]$Genotype_name <- "CS"
mps.t[which(mps.t$Genotype_name == "Eletta Campana"),]$Genotype_name <- "Eletta C."
mps.t[which(grepl("Bial", mps.t$Genotype_name, fixed = TRUE)),]$Genotype_name <- "Bialobrzeskie"
hg2<-ggplot(mps.t, aes(fill=Genus, x=Sample_number, y=Abundance)) +
  theme_minimal() +
  theme(axis.ticks.x = element_line(linewidth = 0.1), 
        axis.ticks.y = element_line(linewidth = 0,1),
        axis.ticks.length = unit(2, "pt")) +
  facet_nested_wrap(~Supplier + Genotype_name, nrow = 2) +
  theme(strip.background = element_rect(fill = "lightgrey", color = "lightgrey"),
        strip.text.x = element_text(size = 12)) +
  geom_bar(position = "fill", stat="identity") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme(legend.text = element_text(size = 14)) + 
  theme(legend.title = element_text(size = 14, face = "bold")) + 
  theme(legend.position = "bottom", legend.box.spacing = margin(-10)) +
  theme(legend.key.size=unit(14,"point")) +
  scale_fill_manual(values= GenusPalette, breaks=legend.order) +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") +
  ylab("Relative abundance (%)") +
  theme(axis.title = element_text(size=14, face = "bold")) +
  theme(axis.text = element_text(size=12)) +
  guides(fill=guide_legend(ncol=3)); hg2
#ggsave("Hemp_Genus_barplot_top10_nested_facets.png", plot = hg2, device = png, height = 24, width = 24, units = "cm", limitsize = TRUE, bg = "white")


############################################################
####step 6.3 - alpha diversity analyses

library(phyloseq)
library(ggplot2)
library(cowplot)
library(microViz)
library(dplyr)
library(RColorBrewer)
library(multcomp)
library(Rmisc)  
library(FSA)
library(rcompanion)
library(Rmisc)  
library(btools)


ps<-readRDS("Hemp_ps_processed_step5_readyForAnalysis.rds")

#use rarefaction curve to inform rarefaction
library(phyloseq.extended)
p <- ggrare(ps, step = 1000, color = "Genotype_name", se = FALSE)  +
  ylab("ASV richness"); p #ASV richness
#ggsave("Hemp_diversity_alpha_rarefaction.png", plot = p, device = png, height = 15, width = 20, units = "cm", limitsize = TRUE)

#obtain min read count to inform rarefication
ps.rc<-as.data.frame(sample_sums(ps)); colnames(ps.rc)<-"read.count"
min(ps.rc$read.count) #21667
hist(ps.rc$read.count)

#rarefy to min read count
rps<-rarefy_even_depth(ps, replace = FALSE, sample.size = 21667, rngseed = 100)

#calculate ASV richness and Shannon index
table_rps <- estimate_richness(rps, split = TRUE, measures = c("Observed", "Shannon"))
data_rps <- cbind(sample_data(rps), table_rps)
data_rps[which(data_rps$Genotype_name=="CS Carmagnola"),]$Genotype_name <- "CS"
rps_obs_hist <- ggplot(data = data_rps, aes(x = Observed)) +
  geom_histogram(binwidth = 5) +
  theme_classic() +
  xlab("Observed ASVs") +
  ylab("# samples") +
  ggtitle("ASV richness per sample"); rps_obs_hist
#ggsave("Hemp_diversity_alpha__observedASVs_histogram.png", plot = rps_obs_hist, device = "png", units = "cm", height = 10, width = 10)
#summary stats on ASV richness and Shannon index
median(data_rps$Observed); mean(data_rps$Observed); min(data_rps$Observed); max(data_rps$Observed)
median(data_rps$Shannon); mean(data_rps$Shannon); min(data_rps$Shannon); max(data_rps$Shannon)

###ASV richness
##statistical tests - anova
data_rps$Genotype_name<-as.factor(data_rps$Genotype_name)
aov.go<-aov(data=data_rps, Observed ~ Genotype_name)
summary(aov.go)
#check assmuptions
  #normality 
qqnorm(aov.go$residuals)
qqline(aov.go$residuals)
  #variance
plot(aov.go, 1) 
##pairwise comparisons - multcomp package
  #reorder by mean so that CLDs will be in right order in plot
    #get order of richness
order_obs<- data_rps %>%
  group_by(Genotype_name) %>%
  dplyr::summarize(Mean.richness = mean(Observed, na.rm=FALSE))
list_order_obs<-order_obs[order(order_obs$Mean.richness, decreasing = TRUE),]$Genotype_name; list_order_obs
    #change factor level order
data_rps2 <- data_rps
data_rps2$Genotype_name<-as.factor(data_rps2$Genotype_name)
data_rps2$Genotype_name<-factor(data_rps2$Genotype_name, levels = list_order_obs)
#pairwise comparisons
model<-lm(Observed ~ Genotype_name, data=data_rps2)
anova(model) 
summary(model) 
mc = glht(model,
          mcp(Genotype_name = "Tukey"))
mcs = summary(mc, test=adjusted("BH"))
mcs
#obtain CLDs
model.sigs<-cld(mcs,
                level=0.05,
                decreasing=FALSE) 
df.cld<-as.data.frame(model.sigs$mcletters$Letters); colnames(df.cld)<-"Letters"
df.cld$Genotype_name <- rownames(df.cld); rownames(df.cld) <- NULL; df.cld 
richness_letters<-df.cld$Letters
#obtain standard errors
Data1 = summarySE(data=data_rps2,
                  "Observed",
                  groupvars="Genotype_name",
                  conf.interval = 0.95)
Tabla = as.table(Data1$Observed)         
rownames(Tabla) = Data1$Genotype_name
Tabla
offset.v = -2.0     # offsets for CLDs
offset.h = 0.5
#obtain info on supplier
Data1$Supplier<-data_rps2[match(paste(Data1$Genotype_name),paste(data_rps2$Genotype_name)),]$Supplier
#plot
ggobs<-ggplot(Data1,
              aes(x=reorder(Genotype_name,-Observed,FUN = mean), y = Observed, fill = Supplier,
                  ymax=90, ymin=0))+
  geom_bar(stat="identity",
           width = 0.7)  +
  geom_errorbar(aes(ymax=Observed+se, ymin=Observed-se),
                width=0.0, linewidth=0.5, color="black")  +
  geom_text(aes(label=richness_letters,
                hjust=offset.h, vjust=offset.v), size = 4) +             
  labs(x = "Hemp genotype",
       y = "ASV richness")  +
  theme_classic()  +
  theme(panel.grid.major.x = element_blank(),
        axis.title.y = element_text(vjust= 1.8, size = 14, face = "bold"),
        axis.title.x = element_text(vjust= -0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)) +
  scale_fill_discrete(labels=c("CREA", "HEMPit", "Hempoint", "PSTS", "Vandinter Semo")) +
  scale_y_continuous(breaks = seq(0, 200, 20), expand = expansion(mult = c(0,0.15))); ggobs
#ggsave("Hemp_diversity_alpha_richness_ASV.png", plot = last_plot(), device = png, height = 10, width = 15, units = "cm", limitsize = TRUE)

###shannons index
##stats tests
data_rps$Genotype_name<-as.factor(data_rps$Genotype_name)
aov.gs<-aov(data=data_rps, Shannon ~ Genotype_name)
summary(aov.gs)
#check assmuptions
  #normality
qqnorm(aov.gs$residuals)
qqline(aov.gs$residuals) 
shapiro.test(data_rps$Shannon) 
  #variance
plot(aov.gs, 1) #no clear relationship
##need to do non-parametric test
  #kruskal wallis test
kruskal.test(Shannon ~ Genotype_name, data = data_rps)
#cldList below doesn't like names with hyphen so changing one name - restored later
data_rps3<-data_rps
data_rps3$Genotype_name<-as.character(data_rps3$Genotype_name)
data_rps3[which(data_rps3$Genotype_name == "Uso-31"),]$Genotype_name <- "Uso31"
data_rps3$Genotype_name<-as.factor(data_rps3$Genotype_name)
#dunn test
dt<-dunnTest(Shannon ~ Genotype_name, data = data_rps3,
             method="bh", list = TRUE)
print(dt, dunn.test.results = TRUE)
#restore name
data_rps3$Genotype_name<-as.character(data_rps3$Genotype_name)
data_rps3[which(data_rps3$Genotype_name == "Uso31"),]$Genotype_name <- "Uso-31"
data_rps3$Genotype_name<-as.factor(data_rps3$Genotype_name)
#obtain CLDs
shannon_letters<-cldList(comparison = dt$res$Comparison,
                         p.value    = dt$res$P.adj,
                         threshold  = 0.05)
shannon_letters<-shannon_letters$Letter
#get standard errors
Data3 = summarySE(data=data_rps,
                  "Shannon",
                  groupvars="Genotype_name",
                  conf.interval = 0.95)
Tabla3 = as.table(Data3$Shannon)         
rownames(Tabla3) = Data3$Genotype_name
Tabla3
offset.v = -2.0     #offsets for CLDs
offset.h = 0.5
#obtain info on supplier
Data3$Supplier<-data_rps[match(paste(Data3$Genotype_name),paste(data_rps$Genotype_name)),]$Supplier
#plot
ggshan<-ggplot(Data3,
               aes(x=reorder(Genotype_name,-Shannon,FUN = mean), y = Shannon, fill = Supplier,
                   ymax=3, ymin=0))+
  geom_bar(stat="identity",
           width = 0.7)  +
  geom_errorbar(aes(ymax=Shannon+se, ymin=Shannon-se),
                width=0.0, linewidth=0.5, color="black")  +
  geom_text(aes(label=shannon_letters, 
                hjust=offset.h, vjust=offset.v), 
            size = 4, angle = 0) +           
  labs(x = "Hemp genotype",
       y = "Shannon diversity index")  +
  theme_classic()  +
  theme(panel.grid.major.x = element_blank(),
        axis.title.y = element_text(vjust= 1.8, size = 14, face = "bold"),
        axis.title.x = element_text(vjust= -0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)) +
  scale_fill_discrete(labels=c("CREA", "HEMPit", "Hempoint", "PSTS", "Vandinter Semo")) +
  scale_y_continuous(breaks = seq(0, 3.5, 0.5), expand = expansion(mult = c(0,0.15))); ggshan
#ggsave("Hemp_diversity_alpha_Shannon_ASV.png", plot = last_plot(), device = png, height = 15, width = 12, units = "cm", limitsize = TRUE)


###Faith's phylogenetic diversity
###stats tests - anova
#calculate Faith's PD
table_PD<-estimate_pd(rps)
pdsd <- cbind(sample_data(rps), table_PD)
pdsd[which(pdsd$Genotype_name=="CS Carmagnola"),]$Genotype_name <- "CS"
pdsd$Genotype_name<-as.factor(pdsd$Genotype_name)
aov.pd<-aov(data=pdsd, PD ~ Genotype_name)
summary(aov.pd)
#check assmuptions
  #normality
qqnorm(aov.pd$residuals)
qqline(aov.pd$residuals) 
#variance
plot(aov.pd, 1) 
###pairwise comparisons - multcomp package
#reorder by mean so that letters will be in right order
#get order of richness
order_pd<- pdsd %>%
  group_by(Genotype_name) %>%
  dplyr::summarize(Mean.PD = mean(PD, na.rm=FALSE))
list_order_pd<-order_pd[order(order_pd$Mean.PD, decreasing = TRUE),]$Genotype_name; list_order_pd
#change factor level order - to ensure CLDs are alphabetical
pdsd2 <- pdsd
pdsd2$Genotype_name<-as.factor(pdsd2$Genotype_name)
pdsd2$Genotype_name<-factor(pdsd2$Genotype_name, levels = list_order_pd)
model<-lm(PD ~ Genotype_name, data=pdsd2)
anova(model)
summary(model)
#pairwise comparisons
mc = glht(model,
          mcp(Genotype_name = "Tukey"))
mcs = summary(mc, test=adjusted("BH"))
mcs
#obtain CLDs
model.sigs<-cld(mcs,
                level=0.05,
                decreasing=FALSE) 
df.cld<-as.data.frame(model.sigs$mcletters$Letters); colnames(df.cld)<-"Letters"
df.cld$Genotype_name <- rownames(df.cld); rownames(df.cld) <- NULL; df.cld 
pd_letters<-df.cld$Letters
#obtain standard errors
Data4 = summarySE(data=pdsd2,
                  "PD",
                  groupvars="Genotype_name",
                  conf.interval = 0.95)
#obtain info on supplier
Data4$Supplier<-pdsd2[match(paste(Data4$Genotype_name),paste(pdsd2$Genotype_name)),]$Supplier
#plot
ggpd<-ggplot(Data4, aes(x = reorder(Genotype_name,-PD, FUN = mean), y = PD, fill = Supplier,
                        ymax = 12, ymin = 0)) + 
  geom_bar(stat="identity",
           width = 0.7) +
  geom_errorbar(aes(ymax=PD+se, ymin=PD-se),
                width=0.0, linewidth=0.5, color="black")  +
  geom_text(aes(label=pd_letters,
                hjust=offset.h, vjust=offset.v), size = 4) +
  labs(x = "Hemp genotype",
       y = "Faith's PD index")  +
  theme_classic()  +
  theme(panel.grid.major.x = element_blank(),
        axis.title.y = element_text(vjust= 1.8, size = 14, face = "bold"),
        axis.title.x = element_text(vjust= -0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)) +
  scale_y_continuous(breaks = seq(0, 20, 2), expand = expansion(mult = c(0,0.15))) +
  scale_fill_discrete(labels=c("CREA", "HEMPit", "Hempoint", "PSTS", "Vandinter Semo")); ggpd
#ggsave("Hemp_diversity_alpha_PD.png", plot = last_plot(), device = png, height = 15, width = 12, units = "cm", limitsize = TRUE)

###plot all diversity metrics together
plot_grid(ggobs, ggshan, ggpd, ncol = 1, scale = 0.90, labels = "auto", label_size = 14, label_fontface = "bold")
#ggsave("Hemp_diversity_alpha_richnessShannonPD.png", plot = last_plot(), device = png, height = 32, width = 24, units = "cm", limitsize = TRUE, bg = "white")


############################################################
####step 6.4 - beta diversity analyses
library(phyloseq)
library(ggplot2)
library(cowplot)
library(microViz)
library(dplyr)
library(vegan)

ps<-readRDS("Hemp_ps_processed_step5_readyForAnalysis.rds")

ps

#using proportions rather than rarefication 
rps <- microbiome::transform(ps, "compositional")
#McKnight et al. reports that proportions work best for these analyses
  #DOI: 10.1111/2041-210X.13115

##PERMANOVA - supplier and genotype
all_ps_bray <- phyloseq::distance(rps, method = "bray") 
all_ps_df <- data.frame(sample_data(rps)); colnames(all_ps_df)
all_ps_df$Supplier[is.na(all_ps_df$Supplier)]<-"Unknown"
ad_all_ps<-adonis2(all_ps_bray ~ Supplier + Genotype_name, data = all_ps_df, permutations = 10000) 
print(ad_all_ps)