##gwas maize he et al. 2023 biorxiv
setwd("C:/Users/plduran/Desktop/proposal/Core microbiota/figure1_withinspecies/he_maize")
design<- read.table("design_root.txt", header=TRUE, sep="\t")
asv_table_bac<- read.table("asv_table_bacteria_root.txt", header=TRUE, sep="\t")
asv_table_fun<- read.table("asv_table_fungi_root.txt", header=TRUE, sep="\t")

tax_bac<- read.table("tax_bacteria.txt",header=TRUE, sep="\t")
tax_fun<- read.table("tax_fungi.txt",header=TRUE, sep="\t")

##calculated first for rs but too many orders were common
#asv_table<- t(asv_table)
#write.table(asv_table, "t_asv_table.txt", sep="\t")

asv_table_norm_bac <- apply(asv_table_bac, 2, function(x) x/sum(x))
asv_table_norm_fun <- apply(asv_table_fun, 2, function(x) x/sum(x))

threshold <- .05
idx <- rowSums(asv_table_norm_bac * 100 > threshold) >= 1
asv_table_bac <- asv_table_bac[idx, ]
asv_table_norm_bac <- asv_table_norm_bac[idx, ]

idx <- tax_bac$ASV %in% rownames(asv_table_bac)
tax_bac<- tax_bac[idx,]


idx <- rowSums(asv_table_norm_fun * 100 > threshold) >= 1
asv_table_fun <- asv_table_fun[idx, ]
asv_table_norm_fun <- asv_table_norm_fun[idx, ]

idx <- tax_fun$ASV %in% rownames(asv_table_fun)
tax_fun<- tax_fun[idx,]


level <- which(colnames(tax_bac)=="Order")
asv_table_order_bac<- aggregate(asv_table_norm_bac, by=list(tax_bac[, level]), FUN=sum)
write.table(asv_table_order_bac, "asv_table_order_bac.txt", sep="\t")

asv_table_order_core_bac<- read.table("asv_table_order_bac_core_root.txt", header=TRUE, sep="\t")
otu_long <- melt(asv_table_order_core_bac, id.vars = "Group.1", variable.name = "Sample")
points <- cbind(otu_long, design[match(otu_long$Sample, design$SampleID), ])
#points<- na.omit(points)
sd<- ddply(points, c("Group.1", "Germplasm", "Genotype"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))

p<- ggplot(sd, aes(x = Genotype, y = mean, fill = Group.1)) +     geom_bar(stat = "identity", colour="black")+
  theme(axis.text.x = element_text(size=10, angle=45))	 
p+facet_grid(~Germplasm, scales="free_x", space="free")

##fungi
threshold <- .05
idx <- rowSums(asv_table_norm_fun * 100 > threshold) >= 1
asv_table_fun<- asv_table_fun[idx, ]
asv_table_norm_fun <- asv_table_norm_fun[idx, ]

idx <- tax_fun$ASV %in% rownames(asv_table_fun)
tax_fun<- tax_fun[idx,]


idx <- rowSums(asv_table_norm_fun * 100 > threshold) >= 1
asv_table_fun <- asv_table_fun[idx, ]
asv_table_norm_fun <- asv_table_norm_fun[idx, ]

idx <- tax_fun$ASV %in% rownames(asv_table_fun)
tax_fun<- tax_fun[idx,]


level <- which(colnames(tax_fun)=="Order")
asv_table_order_fun<- aggregate(asv_table_norm_fun, by=list(tax_fun[, level]), FUN=sum)
write.table(asv_table_order_fun, "asv_table_order_fun.txt", sep="\t")

asv_table_order_core_fun<- read.table("asv_table_order_fun_core_root.txt", header=TRUE, sep="\t")
otu_long <- melt(asv_table_order_core_fun, id.vars = "Group.1", variable.name = "Sample")
points <- cbind(otu_long, design[match(otu_long$Sample, design$SampleID), ])
#points<- na.omit(points)
sd<- ddply(points, c("Group.1", "Germplasm", "Genotype"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))

p<- ggplot(sd, aes(x = Genotype, y = mean, fill = Group.1)) +     geom_bar(stat = "identity", colour="black")+
  theme(axis.text.x = element_text(size=10, angle=45))	 
p+facet_grid(~Germplasm, scales="free_x", space="free")



