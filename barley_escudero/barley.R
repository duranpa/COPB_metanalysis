##gwas switchgrass edwards et al. 2023 current biology
setwd("C:/Users/plduran/Desktop/proposal/Core microbiota/figure1_withinspecies/barley_escudero")
design<- read.table("design.txt", header=TRUE, sep="\t")
escudero<- readRDS("JH15_JH07_rare_ASV_10K.rds")
OTU1 = as(otu_table(escudero), "matrix")
if(taxa_are_rows(escudero)){OTU1 <- t(OTU1)}
OTUdf = as.data.frame(OTU1)
tax<- tax_table(escudero)
taxdf = as.data.frame(tax)
write.table(taxdf, "taxonomy.txt", sep="\t")
write.table(OTUdf, "otu_table.txt", sep="\t")

asv_table<- read.table("otu_table.txt", header=TRUE, sep="\t")
asv_table<- t(asv_table)
write.table(asv_table, "t_asv_table.txt", sep="\t")
asv_table<- read.table("t_asv_table.txt", header=TRUE, sep="\t")

taxonomy<- read.table("taxonomy.txt", header=TRUE, sep="\t")
design<- read.table("Map_JH07_phyloseq_3.txt", header=TRUE, se="\t")


asv_table_norm <- apply(asv_table, 2, function(x) x/sum(x))

threshold <- .05
idx <- rowSums(asv_table_norm * 100 > threshold) >= 1
asv_table <- asv_table[idx, ]
asv_table_norm <- asv_table_norm[idx, ]

idx <- taxonomy$ASVID %in% rownames(asv_table)
taxonomy<- taxonomy[idx,]

level <- which(colnames(taxonomy)=="Order")
asv_table_order<- aggregate(asv_table_norm, by=list(taxonomy[, level]), FUN=sum)
write.table(asv_table_order, "asv_table_order.txt", sep="\t")

asv_table_order_core<- read.table("asv_table_order_core.txt", header=TRUE, sep="\t")
otu_long <- melt(asv_table_order_core, id.vars = "Group.1", variable.name = "Sample")
points <- cbind(otu_long, design[match(otu_long$Sample, design$SampleID), ])
points<- na.omit(points)

sd<- ddply(points, c("Group.1", "Genotype"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))

p<- ggplot(sd, aes(x = Genotype, y = mean, fill = Group.1)) +     geom_bar(stat = "identity", colour="black")+
  theme(axis.text.x = element_text(size=10, angle=45))	 
