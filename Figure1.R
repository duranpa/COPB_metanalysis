#core microbiota combined
#barley
core_escudero<-read.table("asv_table_order_core_reduced.txt", header=TRUE, sep="\t")
long_escudero <- melt(core_escudero, id.vars = "Group.1", variable.name = "Sample")
design_escudero<- read.table("design.txt", header=TRUE, sep="\t")
points_escudero <- cbind(long_escudero, design_escudero[match(long_escudero$Sample, design_escudero$SampleID), ])
points_escudero<- na.omit(points_escudero)
write.table(points_escudero, "points_escudero.txt", sep="\t")
#maize_bac
core_he<-read.table("asv_table_order_bac_core_root_reduced.txt", header=TRUE, sep="\t")
long_he <- melt(core_he, id.vars = "Group.1", variable.name = "Sample")
design_he<- read.table("design_root.txt", header=TRUE, sep="\t")
points_he <- cbind(long_he, design_he[match(long_he$Sample, design_he$SampleID), ])
points_he<- na.omit(points_he)
write.table(points_he, "points_he.txt", sep="\t")

core_he<-read.table("asv_table_order_fun.txt", header=TRUE, sep="\t")
long_he <- melt(core_he, id.vars = "Group.1", variable.name = "Sample")
design_he<- read.table("design_root.txt", header=TRUE, sep="\t")
points_he <- cbind(long_he, design_he[match(long_he$Sample, design_he$SampleID), ])
points_he<- na.omit(points_he)
write.table(points_he, "points_he_fun.txt", sep="\t")

#thaliana_karasov
core_karasov<-read.table("asv_table_order_core.txt", header=TRUE, sep="\t")
long_karasov <- melt(core_karasov, id.vars = "Group.1", variable.name = "Sample")
design_karasov<- read.table("design.txt", header=TRUE, sep="\t")
points_karasov <- cbind(long_karasov, design_karasov[match(long_karasov$Sample, design_karasov$SampleID), ])
points_karasov<- na.omit(points_karasov)
write.table(points_karasov, "points_karasov.txt", sep="\t")

#switchgrass
core_edwards<-read.table("asv_table_norm_order_core_reduced.txt", header=TRUE, sep="\t")
long_edwards <- melt(core_edwards, id.vars = "Group.1", variable.name = "Sample")
design_edwards<- read.table("design.txt", header=TRUE, sep="\t")
points_edwards <- cbind(long_edwards, design_edwards[match(long_edwards$Sample, design_edwards$Plate_site), ])
points_edwards<- na.omit(points_edwards)
write.table(points_edwards, "points_edwards.txt", sep="\t")

#arabidopsis_thiergart_bac
core_thiergart<-read.table("asv_table_order_bac_core_reduced.txt", header=TRUE, sep="\t")
long_thiergart <- melt(core_thiergart, id.vars = "Group.1", variable.name = "Sample")
design_thiergart<- read.table("design_bac.txt", header=TRUE, sep="\t")
points_thiergart <- cbind(long_thiergart, design_thiergart[match(long_thiergart$Sample, design_thiergart$SampleID), ])
points_thiergart<- na.omit(points_thiergart)
write.table(points_thiergart, "points_thiergart.txt", sep="\t")

#fungi
core_thiergart<-read.table("asv_table_order_fun_core_reduced.txt", header=TRUE, sep="\t")
long_thiergart <- melt(core_thiergart, id.vars = "Group.1", variable.name = "Sample")
design_thiergart<- read.table("design_fun.txt", header=TRUE, sep="\t")
points_thiergart <- cbind(long_thiergart, design_thiergart[match(long_thiergart$Sample, design_thiergart$SampleID), ])
points_thiergart<- na.omit(points_thiergart)
write.table(points_thiergart, "points_thiergart_fun.txt", sep="\t")

#merged files manually (because different factors each)
#fungi
combined_core<- read.table("combined_core_fun.txt", header=TRUE, sep="\t")

shapes <- data.frame(group=c("Root"),
                     shape=c(16))		

sd<- ddply(combined_core, c("Article", "Microhabitat","Study_type", "Group.1", "Phylum"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))

p<- ggplot(sd, aes(x = Article, y = Group.1, color = Phylum, shape=Microhabitat, size=mean))+
  geom_point(alpha=0.8) +
  scale_color_viridis(discrete=TRUE)+
  scale_shape_manual(values=shapes$shape) +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1), axis.text.y = element_text(size=10),
        legend.text = element_text(size=10), 
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(size=12), strip.text.y = element_text(size=12, angle=0))+
    #geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+
  labs(y=paste("Average relative abundance (%)"))
p+facet_grid(Phylum~Study_type, scale="free", space="free")

#Bacteria
combined_core<- read.table("combined_core.txt", header=TRUE, sep="\t")

shapes <- data.frame(group=c("Leaf","Rhizosphere", "Root"),
                     shape=c(18,15,16))		

sd<- ddply(combined_core, c("Article", "Microhabitat","Study_type", "Group.1", "Phylum"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))

p<- ggplot(sd, aes(x = Article, y = Group.1, color = Phylum, shape=Microhabitat, size=mean))+
  geom_point(alpha=0.8) +
  scale_color_viridis(discrete=TRUE, option="magma")+
  scale_shape_manual(values=shapes$shape) +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1), axis.text.y = element_text(size=10),
        legend.text = element_text(size=10), 
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(size=12), strip.text.y = element_text(size=12, angle=0))+
  #geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+
  labs(y=paste("Average relative abundance (%)"))
p+facet_grid(Phylum~Study_type, scale="free", space="free")

##calculated the number of core groups in each study/species
number<- read.table("species_number_bac.txt", header=TRUE, sep="\t")

colors <- data.frame(group=c("Multiple genotypes_multiple habitats","Multiple genotypes_one habitat", "Multiple habitats"),
                     color=c("#d8d8d8","#b2b2b2","#8c8c8c"))		

p<- ggplot(number, aes(x = Article, y =number_species , fill = Study_type))+
  geom_bar(stat = "identity", colour="black")+
  scale_fill_manual(values=as.character(colors$color)) +
  theme(axis.text.x = element_text(size=7, angle=45, hjust=1))+
  labs(y=paste("Number of core groups"))
p+facet_grid(~Study_type, scales="free_x", space="free")

number<- read.table("species_number_fun.txt", header=TRUE, sep="\t")

colors <- data.frame(group=c("Multiple genotypes_one habitat", "Multiple habitats"),
                     color=c("#b2b2b2","#8c8c8c"))		

p<- ggplot(number, aes(x = Article, y =number_species , fill = Study_type))+
  geom_bar(stat = "identity", colour="black")+
  scale_fill_manual(values=as.character(colors$color)) +
  theme(axis.text.x = element_text(size=7, angle=45, hjust=1))+
  labs(y=paste("Number of core groups"))
p+facet_grid(~Study_type, scales="free_x", space="free")


