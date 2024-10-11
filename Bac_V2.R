
cran_packages <- c("knitr", "qtl", "bookdown", "ggplot2", "grid", "gridExtra", "tidyverse", "xtable", "devtools", "kableExtra", "remotes")
bioc_packages <- c("phyloseq", "dada2", "DECIPHER", "phangorn", "ggpubr", "DESeq2", "genefilter", "philr")
git_source <- c("twbattaglia/btools", "gmteunisse/Fantaxtic", "MadsAlbertsen/ampvis2", "opisthokonta/tsnemicrobiota")
git_packages <- c("btools", "fantaxtic")

#ls()
#rm(list = ls())
set.seed(135885)
library(readxl)

samples.out <- rownames(seqtab.nochim)

subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
#gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
#day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject)
#samdf$When <- "Early"
#samdf$When[samdf$Day>100] <- "Late"

#Ajuste de directorio de trabajo
path <- "C:/Users/kuris/OneDrive/Escritorio/16S_bac"
setwd(path)#ajusta el directio con el objeto path
#Listando archivos del directorio de trabajo
#list.files(path)#vemos la lista

#Phyloseq
##################################

#Invertir OTU.table
path <- "C:/Users/kuris/OneDrive/Escritorio/16S_bac"
setwd (path)

rownames(samdf) <- samples.out

sample_aa <- read_excel("C:/Users/kuris/OneDrive/Escritorio/16S_bac/Metadata_AN.xlsx")
#sample_aa <- read_excel("C:/Users/kuris/OneDrive/Escritorio/16S_bac/met_2.xlsx")
#load("C:/Users/kuris/OneDrive/Escritorio/16S_bac/pyloseqBact AN.RData")
library(phyloseq)

samdfaa<-cbind(samdf,sample_aa)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE, tax_table(taxonomy)), 
               sample_data(samdfaa), 
               tax_table(taxonomy))
ps
table(tax_table(ps)[, "Kingdom"], exclude=NULL)

ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remover potenciales muestras sintéticas
ps = subset_taxa(ps, Kingdom =="Bacteria")
ps = subset_taxa(ps, Kingdom != "Cyanobacteria" & Kingdom != "Tenericutes")

ps
View(tax_table(ps))

#View(otu_table(ps))
samdfaa

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps2 <- merge_phyloseq(ps, dna)
taxa_names(ps2) <- paste0("ASV", seq(ntaxa(ps2)))
ps2
View(tax_table(ps2))
     
library(ape)
random_tree = rtree(ntaxa(ps2), rooted=TRUE, tip.label=taxa_names(ps2))
#plot(random_tree)
ps2
ps2 = merge_phyloseq(ps2, sample_data, random_tree)
ps2


ps3 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(samdfaa), 
                tax_table(taxonomy),
                phy_tree( random_tree))


ps3
phy_high

#conrol de calidad
# Computamos prevalencia para cada feature y la guardamos en un data frame
prevdf = apply(X = otu_table(ps2),
               MARGIN = ifelse(taxa_are_rows(ps2), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Le agregamos la taxonomía
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps2),
                    tax_table(ps2))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) -> dfprev
library(kableExtra)
kable(dfprev)
#Al examinar la tabla, es evidente que algunos Phylum aunque presentes, están muy poco representados. 
#La columna 1 representa la media de read counts para ese Phylum, mientras que la columna 2 representa la suma.
###Medidas de riqueza, uniformidad, dominancia, diversidad filogenética (diversidad alfa)

# Definimos taxa a filtrar
filterPhyla = c("10bav-F6", "Abditibacteriota", "Calditrichota", "Dadabacteria", "Edwardsbacteria", "Deinococcota", "Fusobacteriota", "GAL15", "Spirochaetes", "NKB15", "PAUC34f", "Schekmanbacteria", "TA06", "NA")

# Procedemos a filtrar
(ps3 = subset_taxa(ps2, !Phylum %in% filterPhyla))



# Además aprovechamos a remover taxa que no corresponde a microorganismos como cloroplastos, mitocondrias y otros

filterPhyla2 <- c("Chloroplast", "Mitochondria", "Eukaryota")
ps3 <- subset_taxa(ps3, !Kingdom %in% filterPhyla2)
ps3 <- subset_taxa(ps3, !Phylum %in% filterPhyla2)
ps3 <- subset_taxa(ps3, !Class %in% filterPhyla2)
ps3 <- subset_taxa(ps3, !Order %in% filterPhyla2)
ps3 <- subset_taxa(ps3, !Family %in% filterPhyla2)
ps3 <- subset_taxa(ps3, !Genus %in% filterPhyla2)


#Además del filtrado que acabamos de realizar, existen otros tipos de filtrado que tienen que ver con la media 
#de read counts por taxa, con la distribución de éstas, y con filtrar muestras bajo un número mínimo de reads.

# Filtramos taxa de acuerdo a un umbral de número medio de _read counts_, en este caso 1e-5
psd2 <- filter_taxa(ps3, function(x) mean(x) > 1e-5, TRUE)
psd2
# También podemos remover taxa que no se observe más de X veces en al menos 10% de las muestras
psd3 <- filter_taxa(psd2, function(x) sum(x > 2) > (0.1*length(x)), TRUE)
psd3
# Y finalmente filtrar muestras con menos de 1000 reads
psd4 = prune_samples(sample_sums(psd3) > 1000, psd3)

psd4

# Seleccionamos las taxa de interés
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(psd4, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps2),color=Phylum)) +
  # Agregamos una línea para nuestro umbral
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# Definimos el umbral de prevalencia a un 5%
(prevalenceThreshold = 0.05 * nsamples(psd4))

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
(psd5 = prune_taxa(keepTaxa, psd4))

# Reemplazamos las secuencias por un nombre genérico
taxa_names(psd5) <- paste0("ASV", seq(ntaxa(psd5)))

psd5
View(tax_table(psd5))
# gráficar la distribución de read counts por número de muestra de forma de tener una idea sobre la distribución de éstas.

sample_sum_df <- data.frame(sum = sample_sums(psd5))
sample_sum_df

library(ggplot2)
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) 

##################################################################################
#Finalmente, calculamos curvas de rarefacción para cada muestra, de manera tal que podamos determinar si la profundidad 
#de secuenciación fue sufuciente o si tal vez necesitemos secuenciar más. En otras palabras, este análisis nos permitiría 
#averiguar si al secuenciar más observaríamos más OTUs o ASVs.

# Primero cargamos algunos scripts de manera remota
scripts <- c("graphical_methods.R",
             "tree_methods.R",
             "plot_merged_trees.R",
             "specificity_methods.R",
             "ternary_plot.R",
             "richness.R",
             "edgePCA.R",
             "copy_number_correction.R",
             "import_frogs.R",
             "prevalence.R",
             "compute_niche.R")
urls <- paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/", scripts)

for (url in urls) {
  source(url)
}


# Y graficamos
p <- ggrare(psd5, step = 100, color = "Treatment", label = "name", se = TRUE)
(p <- p + facet_wrap(~species))


library(kableExtra)

# Esta es la tabla de cuentas o read counts
head(otu_table(psd5)) %>%
  kable(format = "html", col.names = colnames(otu_table(psd5))) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "350px")

(tax_table(psd5)) %>%
  kable(format = "html", col.names = colnames(tax_table(psd5))) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "320px")

plot_tree(psd5, method = "treeonly", ladderize = "left")

View(tax_table(psd5))
############################################################################################################
meta16S.epi = data.frame(sample_data(psd5))
min(sample_sums(psd5))
epi_rare = phyloseq::rarefy_even_depth(psd5, rngseed = 123, replace = F)


div_sp2 = phyloseq::estimate_richness(epi_rare, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))


head(div_sp2)

write.table(div_sp2, file="richness.txt", sep=",")
library(openxlsx)
write.xlsx(div_sp2, file="richness_1.xlsx", sep=",")

mrg.epi = merge(meta16S.epi, div_sp2, by.x="row.names", by.y="row.names")
head(mrg.epi)
write.xlsx(mrg.epi, file="richness.xlsx", sep=",")


epi_adiv2 = mrg.epi %>%
  tidyr::gather(key = "measure", value = "value", Observed:InvSimpson)

# select "Observed" for observed richness, or "Shannon" for shannon index, etc. 
epi.df = epi_adiv2 %>%
  filter(measure == "Observed")

aov.epi = epi.df %>%
  aov(value ~ Straw*Time*Treatment, data=.) 

# see anova result
summary(aov.epi)



# Calcular índices de diversidad alfa
alpha_div <- estimate_richness(psd5)
# Promedio de la diversidad alfa
mean(alpha_div$observed)

# Estimación de la diversidad alfa por especie
alpha_div_species <- estimate_richness(psd5, measures = c("Observed", "Chao1", "Shannon", "ACE", "Simpson"))
alpha_div_species
alpha_div_species <- estimate_richness(ps3, measures = c("Observed", "Chao1", "Shannon", "ACE", "Simpson"))
alpha_div_species
#dopblox
plot_richness(psd5, x= "Treatment", color="Time", measures=c("Shannon", "Simpson", "InvSimpson", "Chao1"))

richness <- estimate_richness(psd5)

richness


#Alpha diversity
#dopblox
library(ggplot2)
library("ggpubr")
plot_richness(psd5, x="Treatment", measures=c("Chao1","Shannon"))+ geom_boxplot(aes(fill = Treatment), alpha=.6)+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))
plot_richness(psd5, x="Treatment", measures=c("ACE","Shannon"))+ geom_boxplot(aes(fill = Treatment), alpha=.6)+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

plot_richness(psd5, x="Treatment", color = "Time", measures=c("observed", "ACE", "Shannon", "Simpson", "InvSimpson"))+ 
  geom_point(size=2)+theme(axis.title.x=element_text(size=3))


newSTorder = c("C2T0", "C2TF", "BIO", "ME", "MAG")
a_my_comparisons <- list(c("C1T0", "C1TF"), c("C2T0", "C2TF"), c("ME", "MAG"), c("MAG", "BIO"), c("BIO", "ME"), c("BIO", "C1TF"), c("ME", "C1TF"), c("MAG", "C1TF"), c("C2TF", "MAG"))
#
palette <- c("#20DE8B","#CCDE8B", "#FFDE8B", "#FFA88B", "#FF6A8B", "#AD8CFF", "#90CFFF")
a_my_comparisons <- list(c("C2T0", "C2TF"), c("ME", "MAG"), c("C2TF", "MAG", c("MAG", "BIO"), c("BIO", "ME")))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
p2 <- plot_richness(psd5, x="Treatment", measures=c("ACE", "Shannon", "Simpson","Evenness"))+ geom_boxplot(aes(fill = Treatment), alpha=.6)+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))
  scale_fill_manual(values=palette)+
  
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)


p2$data$Treatment<-as.character(p2$data$Treatment)#cambiar orde variable x
p2$data$Treatment <- factor(p2$data$Treatment, levels=newSTorder)
p2

newSTorder = c("C1T0", "C1TF","C2T0", "C2TF", "BIO", "ME", "MAG")
a_my_comparisons <- list(c("C1T0", "C1TF"), c("C2T0", "C2TF"), c("ME", "MAG"), c("MAG", "BIO"), c("BIO", "ME"), c("BIO", "C1TF"), c("ME", "C1TF"), c("MAG", "C1TF"), c("C2TF", "MAG"))
#
palette <- c("#20DE8B","#CCDE8B", "#FFDE8B", "#FFA88B", "#FF6A8B", "#AD8CFF", "#90CFFF")
a_my_comparisons <- list(c("C1T0", "C1TF"), c("C2T0", "C2TF"), c("ME", "MAG"), c("MAG", "BIO"), c("BIO", "ME"), c("BIO", "C1TF"), c("ME", "C1TF"), c("MAG", "C1TF"), c("C2TF", "MAG"), c("C2TF", "BIO"), c("C2TF", "ME"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
p2 <- plot_richness(psd5, x="Treatment", measures=c("ACE", "Shannon", "Simpson"))+ geom_boxplot(aes(fill = Treatment), alpha=.6)+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
  scale_fill_manual(values=palette)
  
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)


p2$data$Treatment<-as.character(p2$data$Treatment)#cambiar orde variable x
p2$data$Treatment <- factor(p2$data$Treatment, levels=newSTorder)
p2

a_my_comparisons <- list( c("BIO", "C2T0"), c("C2T0", "MAG"), c("C2T0", "ME"), c("C1TF", "C2TF"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

plot_richness(psd5, x="Treatment", measures="Shannon", color = "Treatment")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
  stat_compare_means(method = "kruskal.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)


library(datasets)
library(ggplot2)
library(multcompView)
library(dplyr)

#TUKEY SHAnNON

anova.sh = aov(richness$Shannon ~ sample_data(mrg.epi)$Treatment)
summary(anova.sh)
tukey = TukeyHSD(anova.sh)
tukey
write.table(TukeyHSD, file="TukeyHSD.txt", sep=",")

library(agricolae)
tukey.test2 <- HSD.test(anova.sh, "sample_data(mrg.epi)$Treatment")
with(mrg.epi, tukey.test2)
tukey.test2

#TUKEY ace

anova.ace = aov(richness$ACE ~ sample_data(mrg.epi)$Treatment)
summary(anova.ace)
tukey = TukeyHSD(anova.ace)
tukey
write.table(TukeyHSD, file="TukeyHSD.txt", sep=",")

library(agricolae)
tukey.test2 <- HSD.test(anova.ace, "sample_data(mrg.epi)$Treatment")
with(mrg.epi, tukey.test2)
tukey.test2


kruskal.test(richness$Shannon ~ sample_data(mrg.epi)$Treatment)
kruskal.test(richness$Observed ~ sample_data(mrg.epi)$Treatment)
kruskal.test(richness$ACE ~ sample_data(mrg.epi)$Treatment)

#TUKEY SIMSOMP

anova.sim = aov(richness$Simpson ~ sample_data(mrg.epi)$Treatment)
summary(anova.sim)
tukey = TukeyHSD(anova.sim)
tukey
write.table(TukeyHSD, file="TukeyHSD.txt", sep=",")

library(agricolae)
tukey.test2 <- HSD.test(anova.sim, "sample_data(mrg.epi)$Treatment")
with(mrg.epi, tukey.test2)
tukey.test2


hist(richness$Shannon, main="Shannon index", xlab="")



#Wilcox test Shannon, Observed and ACE


pairwise.wilcox.test(richness$Shannon, sample_data(mrg.epi)$Treatment, p.adj = "bonf")
pairwise.wilcox.test(richness$Observed, sample_data(mrg.epi)$Treatment, p.adj = "bonf")
pairwise.wilcox.test(richness$ACE, sample_data(mrg.epi)$Treatment, p.adj = "bonf")




###Diversidad beta y escalamiento multidimensional (Bray-Curtis, UniFrac)###


#phyloseq con tree


sample_aa <- read_excel("C:/Users/kuris/OneDrive/Escritorio/16S_bac/Metadata_AN.xlsx")
#sample_aa <- read_excel("C:/Users/kuris/OneDrive/Escritorio/16S_bac/met_2.xlsx")
samdfaa<-cbind(samdf,sample_aa)

random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
#plot(random_tree)
ps2
ps1 = merge_phyloseq(ps, sample_data, random_tree)
ps1


ps3 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(samdfaa), 
                tax_table(taxonomy),
                phy_tree( random_tree))


ps3

relab_genera = transform_sample_counts(epi_rare, function(x) x / sum(x) * 100) 
pslog <- transform_sample_counts(psd5, function(x) log(1 + x))
#head(otu_table(relab_genera)[,1:3])


#nMDS#


wUF.ordu = ordinate(relab_genera, method="NMDS", distance="unifrac", weighted=TRUE)
wUF.ordu<- as.matrix(wUF.ordu)

### Double Principle Coordinate Analysis, Non-metric MultiDimenstional Scaling, y MDS/PCoA.###

library(microbiome)
library("plyr")
library(phyloseq)
wunifrac_dist <- wunifrac(relab_genera)
ordination = ordinate(relab_genera, method="NMDS", distance= wunifrac_dist)
plot_ordination(relab_genera, ordination, color="Treatment") + theme(aspect.ratio=1)





psd5.mds.unifrac <- ordinate(relab_genera, method = "MDS", distance = "unifrac")
evals <- psd5.mds.unifrac$values$Eigenvalues
pord1 <- plot_ordination(relab_genera, psd5.mds.unifrac, color = "Treatment") +
  labs(col = "geo_loc_name") +
  coord_fixed(sqrt(evals[2] / evals[1]))

psd5.mds.bray <- ordinate(relab_genera, method = "MDS", distance = "bray")
evals <- psd5.mds.bray$values$Eigenvalues
pord2 <- plot_ordination(relab_genera, psd5.mds.bray, color = "Treatment") +
  labs(col = "Treatment") +
  coord_fixed(sqrt(evals[2] / evals[1]))
library(gridExtra)
grid.arrange(pord1, pord2)



#Unifrac unweighted#


wunifrac_dist = phyloseq::distance(relab_genera, method="unifrac", weighted=F)
wunifrac_dist <- as.matrix(wunifrac_dist)
head(wunifrac_dist)[,1:6]

write.table(wunifrac_dist, file="wunifrac_dist.txt", sep=",")



wunifrac_dist = phyloseq::distance(relab_genera, method="unifrac", weighted=F)
ordination = ordinate(relab_genera, method="PCoA", distance=wunifrac_dist)
plot_ordination(relab_genera, ordination, color="Treatment") + theme(aspect.ratio=1)




ord = ordinate(relab_genera, method="PCoA", distance = "unifrac")


plot_ordination(relab_genera, ord, color = "Treatment", shape="Time") + 
  geom_point(size=3) + 
  stat_ellipse(aes(group=Treatment))
  
                



#Permanova#

options(repos = c(CRAN = "https://cran.irsn.fr/"))

install.packages("vegan")
library(vegan)
samples <- data.frame(sample_data(relab_genera))
adonis2(wunifrac_dist ~ Treatment, data = samples)





#Unifrac weighted#


pslog <- microbiome::transform(relab_genera, "compositional")
wunifrac_dist = phyloseq::distance(pslog, method="unifrac", weighted=T)
ord = ordinate(pslog, method="PCoA", distance = "unifrac", weighted=T)



plot_ordination(pslog, ord, color = "Treatment", shape="Time") + 
  geom_point(size=3) + 
  stat_ellipse(aes(group=Treatment))



#Permanova#


samples_b <- data.frame(sample_data(relab_genera))
adonis2(wunifrac_dist ~ Treatment, data = samples_b)



#Bray-curtis#


abrel_bray <- phyloseq::distance(relab_genera, method = "bray")
abrel_bray <- as.matrix(abrel_bray)
head(abrel_bray)[,1:6]

write.table(abrel_bray, file="abrel_bray.txt", sep=",")


ord = ordinate(relab_genera, method="PCoA", distance = "bray")


plot_ordination(relab_genera, ord, color = "Treatment", shape="Time") + 
  geom_point(size=3) + 
  stat_ellipse(aes(group=Treatment))
beta=plot_ordination(relab_genera, ord, color = "Treatment", shape="Time", title="Permanova F-value:5.0532 - R-square:0.60254 - p-value:0.001") + 
  geom_point(size=5) + 
  stat_ellipse(aes(group=Treatment))+
  theme_bw()+
  theme (axis.text.x = element_text(size=rel(1.5)))+
  theme (axis.text.y = element_text(size=rel(1.5)))+
  theme (legend.title = element_text(size = 15))+
  theme (legend.text = element_text(size = 15))+
  theme (axis.title = element_text(size=rel(1.5)))+
  theme(plot.title = element_text(size= 15))
beta
grid.arrange(p2,beta)

#Permanova test using the adonis function from the vegan package.##
#Bray test

library(tidyverse)
samples_b <- data.frame(sample_data(relab_genera))
adonis2(abrel_bray ~ Treatment, data = samples_b)
samples_b 


devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
#por grupos

group=pairwise.adonis(abrel_bray, factors = samples_b$Treatment)
group

write.table(group, file="Permanova_bacteria.txt", sep=",")

#Homogeneidad#


ps.disper <- betadisper(as.dist(abrel_bray), samples$Treatment)
permutest(ps.disper, pairwise = TRUE)


##Abundancias relativas##


library(tidyverse)
library(magrittr)


ps_genus = subset_taxa(psd5, Kingdom == "Bacteria")
ps_genus = prune_taxa(taxa_sums(ps_genus) > 1000, psd5)
# Merge sample by their sampling site
merged_ps = merge_samples(psd5, "name") # summed
# Merge taxa by 'Phylum' rank
normalized_ps_phylum = tax_glom(merged_ps, "Phylum", NArm = TRUE)
# Calculate relative abundance
normalized_ps_phylum_relabun = transform_sample_counts(normalized_ps_phylum, function(OTU) OTU/sum(OTU))

# Show taxa_sums
taxaSums = data.frame(tax_table(normalized_ps_phylum_relabun)[,"Phylum"],
                      taxa_sums = taxa_sums(normalized_ps_phylum_relabun)) %>%
  arrange(desc(taxa_sums)) # reverse sort
taxaSums


# Select top10 phyla
top10 = head(rownames(taxaSums), 10)
top10
# Remove unselected phyla
y = prune_taxa(top10, normalized_ps_phylum_relabun)
y
# Build data.frame
df1 = data.frame(ID = c(taxa_names(y), "Other"), Phylum = c(tax_table(y)[,"Phylum"], "Other"))
df1

df2 = t(cbind(otu_table(y), data.frame(Other = 1 - sample_sums(y))))
df2

df = cbind(df1, df2) %>%
  pivot_longer(-c(ID, Phylum), names_to = "SampleType", values_to = "Abundance") %>%
  as.data.frame
head(df)

# Re-level Phylum by abundance, "Other" placed at the last
df$Phylum = as.factor(df$Phylum)
levels(df$Phylum) # Before ordering

df$Phylum = factor(df$Phylum, levels = c(taxaSums[top10,"Phylum"], "Other"))
levels(df$Phylum) # After ordering

ggplot(df, aes(SampleType, Abundance, fill = Phylum)) +
  geom_col(color = "black") + scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Sample type", y = "Average relative abundance")



ps_genus = subset_taxa(psd5, Kingdom == "Bacteria")
ps_genus = prune_taxa(taxa_sums(psd5) > 1000, ps_genus)
# Merge sample by their sampling site
merged_ps = merge_samples(ps_genus, "name") # summed
# Merge taxa by 'Phylum' rank
normalized_ps_phylum = tax_glom(merged_ps, "Genus", NArm = TRUE)
# Calculate relative abundance
normalized_ps_phylum_relabun = transform_sample_counts(normalized_ps_phylum, function(OTU) OTU/sum(OTU))

# Show taxa_sums
taxaSums = data.frame(tax_table(normalized_ps_phylum_relabun)[,"Genus"],
                      taxa_sums = taxa_sums(normalized_ps_phylum_relabun)) %>%
  arrange(desc(taxa_sums)) # reverse sort
taxaSums


# Select top10 phyla
top10 = head(rownames(taxaSums), 10)
top10
# Remove unselected phyla
y = prune_taxa(top10, normalized_ps_phylum_relabun)
y
# Build data.frame
df1 = data.frame(ID = c(taxa_names(y), "Other"), Phylum = c(tax_table(y)[,"Genus"], "Other"))
df1

df2 = t(cbind(otu_table(y), data.frame(Other = 1 - sample_sums(y))))
df2

df = cbind(df1, df2) %>%
  pivot_longer(-c(ID, Phylum), names_to = "SampleType", values_to = "Abundance") %>%
  as.data.frame
head(df)

# Re-level Phylum by abundance, "Other" placed at the last
df$Phylum = as.factor(df$Phylum)
levels(df$Phylum) # Before ordering

df$Phylum = factor(df$Phylum, levels = c(taxaSums[top10,"Genus"], "Other"))
levels(df$Phylum) # After ordering

ggplot(df, aes(SampleType, Abundance, fill = Phylum)) +
  geom_col(color = "black") + scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Sample type", y = "Average relative abundance")




#Familia#

ps_fa = subset_taxa(ps2, Kingdom == "Bacteria")
ps_fa = prune_taxa(taxa_sums(psd5) > 1000, ps_fa)
# Merge sample by their sampling site
merged_ps = merge_samples(psd5, "name") # summed
# Merge taxa by 'Phylum' rank
normalized_ps_fa = tax_glom(merged_ps, "Family", NArm = TRUE)
# Calculate relative abundance
normalized_ps_fa_relabun = transform_sample_counts(normalized_ps_fa, function(OTU) OTU/sum(OTU))

# Show taxa_sums
taxaSums = data.frame(tax_table(normalized_ps_phylum_relabun)[,"Family"],
                      taxa_sums = taxa_sums(normalized_ps_phylum_relabun)) %>%
  arrange(desc(taxa_sums)) # reverse sort
taxaSums


# Select top10 phyla
top10 = head(rownames(taxaSums), 10)
top10
# Remove unselected phyla
y = prune_taxa(top10, normalized_ps_fa_relabun)
y
# Build data.frame
df1 = data.frame(ID = c(taxa_names(y), "Other"), Family = c(tax_table(y)[,"Family"], "Other"))
df1

df2 = t(cbind(otu_table(y), data.frame(Other = 1 - sample_sums(y))))
df2

df = cbind(df1, df2) %>%
  pivot_longer(-c(ID, Family), names_to = "SampleType", values_to = "Abundance") %>%
  as.data.frame
head(df)

# Re-level Phylum by abundance, "Other" placed at the last
df$Family = as.factor(df$Family)
levels(df$Family) # Before ordering


ggplot(df, aes(SampleType, Abundance, fill = Family)) +
  geom_col(color = "black") + scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Sample type", y = "Average relative abundance")


#phylum top 15·

ps_genus = subset_taxa(psd5, Kingdom == "Bacteria")
ps_genus = prune_taxa(taxa_sums(ps_genus) > 1000, ps_genus)
# Merge sample by their sampling site
merged_ps = merge_samples(psd5, "Treatment") # summed
# Merge taxa by 'Phylum' rank
normalized_ps_phylum = tax_glom(merged_ps, "Phylum", NArm = TRUE)
# Calculate relative abundance
normalized_ps_phylum_relabun = transform_sample_counts(normalized_ps_phylum, function(OTU) OTU/sum(OTU))

# Show taxa_sums
taxaSums = data.frame(tax_table(normalized_ps_phylum_relabun)[,"Phylum"],
                      taxa_sums = taxa_sums(normalized_ps_phylum_relabun)) %>%
  arrange(desc(taxa_sums)) # reverse sort
taxaSums

# Select top15 phyla
top15 = head(rownames(taxaSums), 15)
top15
# Remove unselected phyla
y = prune_taxa(top15, normalized_ps_phylum_relabun)
y
# Build data.frame
df1 = data.frame(ID = c(taxa_names(y), "Other"), Phylum = c(tax_table(y)[,"Phylum"], "Other"))
df1

df2 = t(cbind(otu_table(y), data.frame(Other = 1 - sample_sums(y))))
df2

df = cbind(df1, df2) %>%
  pivot_longer(-c(ID, Phylum), names_to = "SampleType", values_to = "Abundance") %>%
  as.data.frame
head(df)

# Re-level Phylum by abundance, "Other" placed at the last
df$Phylum = as.factor(df$Phylum)
levels(df$Phylum) # Before ordering

df$Phylum = factor(df$Phylum, levels = c(taxaSums[top15,"Phylum"], "Other"))
levels(df$Phylum) # After ordering

ggplot(df, aes(SampleType, Abundance, fill = Phylum)) +
  geom_col(color = "black") + scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Sample", y = "Average relative abundance")




###########################################################################################################
# Get top 20 genera
top20 <- names(sort(taxa_sums(normalized_ps_fa_relabun), decreasing=TRUE))[1:20]
physeq_top20 <- transform_sample_counts(psd5, function(OTU) OTU/sum(OTU)*100)
physeq_top20 <- prune_taxa(top20, physeq_top20)

# Convert into dataframe
taxa_abundance_table_genus <- psmelt(physeq_top20)

StackedBarPlot_genus <- taxa_abundance_table_genus %>% 
  ggplot(aes(x =name, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity",position="fill") +
  
  scale_fill_manual(values=c("#8AD6E4", "#E5A8A0","#B6EEE2","#CAB8E5","#C6EDC5","#A1BDE6","#EBD4AB","#A0CDE2","#CFB192","#87C4B8","#E6BBC7","#B2C69A","#E2C5E0","#A3C0A6",
                             "#B5B4C8","#DBE7C5","#A4BEB8", "#F0DCD0","#CDE4E9","#D0BCAD")) +
  labs(x = "Samples",
       y = "Relative Abundance",
       title = "Genus Relative Abundance") +
  facet_grid(~ Treatment, scales = "free") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10, face = 'italic'),
    strip.text = element_text(size = 12)
  )

StackedBarPlot_genus

plot_heatmap(physeq_top20, sample.label="Treatment", sample.order="Treatment", taxa.label="Genus") # base plot
plot_heatmap(physeq_top20, sample.label="Treatment", sample.order="name", taxa.label="Genus")+
  scale_fill_distiller("Abundance", palette = "RdYlBu") + theme_bw()+ # Change color
  theme(axis.text.y = element_text(colour = 'black', 
                                   size = 10, 
                                   face = 'italic')) + # Make bacterial names italics
  facet_grid(~Treatment, scales = "free") + rremove("x.text")  # Make seperate samples based on main varaible
################################################################################################################

# Get top 20 genera
top20 <- names(sort(taxa_sums(normalized_ps_fa_relabun), decreasing=TRUE))[1:20]
physeq_top20 <- transform_sample_counts(ps_genus, function(OTU) OTU/sum(OTU)*100)
physeq_top20 <- prune_taxa(top20, physeq_top20)

# Convert into dataframe
taxa_abundance_table_genus <- psmelt(physeq_top20)

StackedBarPlot_genus <- taxa_abundance_table_genus %>% 
  ggplot(aes(x =name, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity",position="fill") +
  
  scale_fill_manual(values=c("#8AD6E4", "#E5A8A0","#B6EEE2","#CAB8E5","#C6EDC5","#A1BDE6","#EBD4AB","#A0CDE2","#CFB192","#87C4B8","#E6BBC7","#B2C69A","#E2C5E0","#A3C0A6",
                             "#B5B4C8","#DBE7C5","#A4BEB8", "#F0DCD0","#CDE4E9","#D0BCAD")) +
  labs(x = "Samples",
       y = "Relative Abundance",
       title = "Genus Relative Abundance") +
  facet_grid(~ Treatment, scales = "free") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10, face = 'italic'),
    strip.text = element_text(size = 12)
  )

StackedBarPlot_genus

plot_heatmap(physeq_top20, sample.label="Treatment", sample.order="Treatment", taxa.label="Phylum") # base plot
plot_heatmap(physeq_top20, sample.label="Treatment", sample.order="name", taxa.label="Phylum")+
  scale_fill_distiller("Abundance", palette = "RdYlBu") + theme_bw()+ # Change color
  theme(axis.text.y = element_text(colour = 'black', 
                                   size = 10, 
                                   face = 'italic')) + # Make bacterial names italics
  facet_grid(~Treatment, scales = "free") + rremove("x.text")  # Make seperate samples based on main varaible

###################################################################################################################

##Tree por generos abundantes##


Thiobacillus = subset_taxa(psd5, Genus == "Thiobacillus")
Thiobacillus

plot_tree(Thiobacillus, plot.margin = 0.5, ladderize = TRUE)

plot_tree(Thiobacillus, color = "Time", label.tips = "Treatment", size = "abundance", plot.margin = 0.5, ladderize = TRUE)

plot_tree(Thiobacillus, color = "Treatment", label.tips = "Time", size = "abundance", plot.margin = 0.5, ladderize = TRUE)





Desulfoba = subset_taxa(ps3, Genus == "Desulfobacca")
Desulfoba
plot_tree(Desulfoba, color = "Time", label.tips = "Treatment", size = "abundance", plot.margin = 0.5, ladderize = TRUE)



clorstridium = subset_taxa(ps3, Genus == "Clostridium sensu stricto 1")
clorstridium
plot_tree(clorstridium, color = "Time", label.tips = "Treatment", size = "abundance", plot.margin = 0.5, ladderize = TRUE) 

plot_tree(ps3, color = "Treatment", shape = "Family", label.tips = "Genus", 
          size = "abundance", plot.margin = 0.5, ladderize = TRUE)

##Mapas de calor##

library(microbiomeutilities)
library(viridis)
library(RColorBrewer)
library(microViz)

#ASV prevalencia

heat.sample <- plot_taxa_heatmap(psd5, taxonomic.level = "Genus", subset.top = 10,
                                 VariableA = "Treatment",
                                 heatcolors = brewer.pal(10, "Blues"),
                                 transformation = "log10")




heat.sample <- plot_taxa_heatmap(psd5, subset.top = 10,
                                 VariableA = "Time",
                                 heatcolors = brewer.pal(10, "Blues"),
                                 transformation = "log10")



ps1.rel <- microbiome::transform(psd5, "compositional")
ps.rel.f <- format_to_besthit(ps1.rel)

#Set different detection levels and prevalence
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- 10^seq(log10(1e-3), log10(.2), length = 9)
#(1e-3) = 0.001% abundance; change "-3" to -2 to increase to 0.01%
detections <- seq(from = 20, to = round(max(abundances(psd5))/10, -1), by = 100)


p <- plot_core(ps.rel.f, plot.type = "heatmap", 
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = 0.6, 
               prevalences = prevalences, 
               detections = detections) +
  xlab("Detection Threshold (Relative Abundance (%))")
print(p)

ps1.f <- format_to_besthit(psd5)

top_otu <- top_taxa(ps1.f, 10)

print(top_otu)

p <- plot_select_taxa(ps1.f, select.taxa= top_otu, "Treatment", "Paired", plot.type = "stripchart")
p

library(dplyr)
psd5 %>%
  tax_transform("compositional", rank = "Genus") %>%
  comp_heatmap()



install.packages("gplots")
library(gplots)
heatmap.2(as.matrix(abrel_bray))
heatmap(as.matrix(abrel_bray))


heatmap_object <- heatmap.2(as.matrix(abrel_bray), col=heat.colors(20),
                            col.clust="black", Rowv=NA, Colv=NA,
                            key=TRUE, breaks=seq(0, 1, 0.05),
                            main="Bray-Curtis dissimilarity heatmap")

# Create a heatmap object
heatmap_object <- heatmap(as.matrix(abrel_bray))

# Add the color scale
heatmap_object <- key(heatmap_object, dendrogram="none",
                      breaks=seq(0, 1, 0.05), col=heat.colors(20))

# Add the histogram
heatmap_object <- hist(as.matrix(abrel_bray), main="Histogram of Bray-Curtis dissimilarity",
                       ylab="Frequency", xlab="Bray-Curtis dissimilarity")

# Add the barplot
heatmap_object <- barplot(colMeans(as.matrix(abrel_bray)), main="Mean Bray-Curtis dissimilarity by sample",
                          ylab="Mean Bray-Curtis dissimilarity", xlab="Sample")

# Plot the heatmap object
plot(heatmap_object)

# Plot the heatmap object
plot(heatmap_object)

###Redes###
#Diversidad beta


plot_net(psd5, color="Treatment", distance="bray")

#####################################################################################################
library(ampvis2)

psd6= transform_sample_counts(psd5, function(x) x / sum(x) * 100) 

psrare_r2_relabT <- psmelt(psd6)

psrare_r2_relabT <- psmelt(psd5)
library(writexl)

write_xlsx(psrare_r2_relabT, "psrare_r2_relabT.xlsx") 



psd6
# Necesitamos extraer la tabla de read counts y la tabla de taxonomía del objeto ps3
# Generamos una copia para no sobreescribir ps3
obj <- psd5
obj
# Cambiamos la orientación de la otu_table
t(otu_table(obj)) -> otu_table(obj)
# Extraemos las tablas
otutable <- data.frame(OTU = rownames(phyloseq::otu_table(obj)@.Data),
                       phyloseq::otu_table(obj)@.Data,
                       phyloseq::tax_table(obj)@.Data,
                       check.names = FALSE
)

# Extraemos la metadada
metadata <- data.frame(phyloseq::sample_data(obj), 
                       check.names = FALSE
)


# ampvis2 requiere que 1) los rangos taxonómicos sean siete y vayan de Kingdom a Species y 2) la primera columna de la metadata sea el identificador de cada muestra
# Entonces duplicamos la columna Género y le cambiamos el nombre a Especie
otutable$Species = otutable$Genus
# Reordenamos la metadata
metadata <- metadata[,c("name","Treatment","Time","Straw")]

my_tree <- phy_tree(psd6)

ps1.av2 <- amp_load(otutable = otutable, metadata = metadata, tree = my_tree)

ps1.av2

# finalmente generamos el objeto ampvis


av2 <- amp_load(otutable, metadata,)




amp_rankabundance(av2, log10_x = T, group_by = "Treatment")


amp_heatmap(ps1.av2, 
            group_by = "Treatment",
            facet_by = "Time", 
            plot_values = TRUE,
            tax_show = 20,
            tax_aggregate = "Genus",
            tax_add = "Phylum",
            plot_colorscale = "sqrt",
            plot_legendbreaks = c(1, 5, 10))

mp_heatmap(ps1.av2, 
           group_by = "Treatment",
           facet_by = "Time", 
           plot_values = TRUE,
           tax_show = 20,
           tax_aggregate = "Class",
           tax_add = "Phylum",
           plot_colorscale = "sqrt",
           plot_legendbreaks = c(1, 5, 10))

amp_heatmap(ps1.av2, 
            group_by = "Treatment",
            facet_by = "Time", 
            plot_values = TRUE,
            tax_show = 10,
            tax_aggregate = "Phylum",
            
            plot_colorscale = "sqrt",
            plot_legendbreaks = c(1, 5, 10))+
  theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
        axis.text.y = element_text(size=12),
        legend.position="right")

amp_heatmap(ps1.av2, 
            group_by = "Treatment",
            facet_by = "Time", 
            plot_values = TRUE,
            tax_show = 10,
            tax_aggregate = "Phylum",
            tax_add = "Genus",
            plot_colorscale = "sqrt",
            plot_legendbreaks = c(1, 5, 10))+
  theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
        axis.text.y = element_text(size=12),
        legend.position="right")

a= amp_heatmap(av2, 
            group_by = "name",
            facet_by = "Treatment", 
            plot_values = TRUE,
            tax_show = 10,
            tax_aggregate = "Class",
            tax_add = "Phylum",
            plot_colorscale = "sqrt",
            plot_legendbreaks = c(1, 5, 10))

grid.arrange(a, b)


TGroup <- tax_glom(psd5, taxrank = "Phylum")
PGroup <- transform_sample_counts(TGroup, function(x)100* x / sum(x))
OTUg <- otu_table(PGroup)
TAXg <- tax_table(PGroup)[,"Phylum"]
AverageD <- as.data.frame(rowMeans(OTUg))
names(AverageD) <- c("Mean")
GTable <- merge(TAXg, AverageD, by=0, all=TRUE)
GTable$Row.names = NULL
GTable <- GTable[order(desc(GTable$Mean)),]

mean_PGroup = sapply(levels(SampleType),function(i){
  rowMeans(otu_table(PGroup)[,SampleType==i])
})

phy = tax_table(PGroup)[rownames(mean_PGroup ),"Phylum"]
rownames(mean_PGroup) = phy
head(sort(rowMeans(mean_PGroup),decreasing=TRUE))
head(GTable)
View(GTable)
####################################################################################
library("DESeq2")

# Creamos un objeto DESeq2 con la función `phyloseq_to_deseq2
diagdds = phyloseq_to_deseq2(psd5, ~Treatment)
# Calculamos los factores de tamaño como parte de la normalización de las muestras
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds, fitType = "parametric")
diagdds <- nbinomWaldTest(diagdds)

# Normalizamos y realizamos el test paramétrico de Wald para determinar taxa diferencialmente abundante.
diagdds = DESeq(diagdds, test="Wald", fitType="local")


res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps3)[rownames(sigtab), ], "matrix"))
#head(sigtab)
write.table(sigtab, file="Deseq2.txt", sep=",")
view(sigtab)


# Guardamos los resultados en el objeto res
res = results(diagdds, cooksCutoff = FALSE)
# hacemos un poco de aseo y ordenamos la tabla de resultados según p-value, y dejamos los valores NA al final
res = res[order(res$padj, na.last=NA), ]

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
view(posigtab)
write.table(posigtab, file="log2Fold.txt", sep=",")

# Manipulaciones varias para finalmente graficar los resultados
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size = 10), axis.text.y = element_text(size = 13), legend.text = element_text(size = 13) )

####################################################################################

library(devtools)
devtools::install_github("zdk123/SpiecEasi")
#Installing package: SpiecEasi
remotes::install_github("zdk123/SpiecEasi")
remotes::install_cran("SpiecEasi")
library(remotes)
# install_github("zdk123/SpiecEasi")
library(SpiecEasi)
se.mb.psd4 <- spiec.easi(psd4, method='mb', lambda.min.ratio=1e-2,
                         nlambda=20, icov.select.params=list(rep.num=50))
ig2.mb <- adj2igraph(getRefit(se.mb.psd4),  vertex.attr=list(name=taxa_names(psd4)))
plot_network(ig2.mb, psd4, type='taxa', color="Phylum")

########################################################################################
#efse_analysis_cleaned_prevalence
#Loading packages and the second phyloseq object
set.seed(134591)
library(tidyverse)

library(phyloseq)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiomeMarker")


if (!requireNamespace("remotes", quietly=TRUE))
install.packages("remotes")
remotes::install_github("yiluheihei/microbiomeMarker")
devtools::install_github("yiluheihei/microbiomeMarker")

library(microbiomeMarker)
ps3 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(samdfaa), 
                tax_table(taxonomy),
                phy_tree( random_tree))

library(phyloseq)
library("scales")
suppressPackageStartupMessages(library(microViz))
set.seed(6692)

# Filtramos taxa de acuerdo a un umbral de número medio de _read counts_, en este caso 1e-5
lfps <- filter_taxa(psd5, function(x) mean(x) > 1e-5, TRUE)

# También podemos remover taxa que no se observe más de X veces en al menos 10% de las muestras
psd3 <- filter_taxa(lfps, function(x) sum(x > 2) > (0.1*length(x)), TRUE)

# Y finalmente filtrar muestras con menos de 1000 reads
psd4 = prune_samples(sample_sums(psd3) > 1000, psd3)

psd4

# Seleccionamos las taxa de interés
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(psd4, "Phylum"))


# Definimos el umbral de prevalencia a un 5%
(prevalenceThreshold = 0.05 * nsamples(psd4))

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
(psd5lf = prune_taxa(keepTaxa, psd4))
psd5lf


#Remove taxonomic levels filled with NAs

psd7 <- tax_fix(
  psd6,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
)


set.seed(1245)

psd5lf = transform_sample_counts(psd5, function(x) x / sum(x) * 100) 
psd5lf
phy_high<-subset_samples(psd5lf, Treatment %in% c("BIO","C2TF"))

phy_high


table(sample_data(phy_high)$Treatment)
group <- sample_data(phy_high)$Treatment
#Using CPM normalization is a common practice in microbiome analysis, and it is recommended to use 
#it when comparing the abundance of different taxa across multiple samples. 
#This ensures that the results are more reliable and comparable.
set.seed(33699)

lef_out<-run_lefse(phy_high, group = "Treatment", norm = "CPM",
                   taxa_rank = "Genus",
                   kw_cutoff = 0.05, lda_cutoff = 3)

lef_out
library(tidyverse)
library(knitr)

dat<-marker_table(lef_out) %>% data.frame() %>% select(1:5)

head(dat)
write.table(dat, file="lef_BIO_C2TF.txt", sep=",")

dat %>% filter(enrich_group=="BIO") %>% head()
dat %>% kable(align = "c")
view(dat)
lefBIO=plot_ef_bar(lef_out)
lefBIO

plot_ef_bar+ scale_fill_manual(values = c("C2TF" = "grey", "BIO" = "red"))
p_abd <- plot_abundance(lef_out, group = "Treatment")
p_abd
p_abd + scale_fill_manual(values = c("C2TF" = "grey", "BIO" = "red"))

mm_edger <- run_edger(
  phy_high,
  group = "Class",
  pvalue_cutoff = 0.1,
  p_adjust = "fdr"
)
mm_edger

plot_cladogram(lef_out, color = c(C2TF = "darkgreen", BIO = "red")) +
  theme(plot.margin = margin(0, 0, 0, 0))
plot_heatmap(lef_out, transform = "log10p")

pht <- run_posthoc_test(psd5, group = "Class")
pht

p_pht <- plot_postHocTest(pht, feature = "p__Bacteroidetes|g__Bacteroides")
p_pht

table(tax_table(psd7)[, "Kingdom"], exclude=NULL)

plot_cladogram(lef_out, color = c("red","blue", "#FFB90F", "#76EEC6", "darkslategray1", "deeppink1", "#EEEEE0"))
plot_cladogram(lef_out, color = c("red","blue"), clade_label_level = 4)
############################################################################

phy_high_1<-subset_samples(psd5lf, Treatment %in% c("ME","C2TF"))

phy_high_1
table(sample_data(phy_high_1)$Treatment)

#Using CPM normalization is a common practice in microbiome analysis, and it is recommended to use 
#it when comparing the abundance of different taxa across multiple samples. 
#This ensures that the results are more reliable and comparable.

lef_out_1<-run_lefse(phy_high_1, group = "Treatment", norm = "CPM",
                   taxa_rank = "Genus",
                 
                   kw_cutoff = 0.05, lda_cutoff = 3)

lef_out_1
library(tidyverse)
library(knitr)

dat<-marker_table(lef_out_1) %>% data.frame() %>% select(1:5)

head(dat)
write.table(dat, file="lef_ME_C2TF.txt", sep=",")

dat %>% filter(enrich_group=="ME") %>% head()
dat %>% kable(align = "c")

view(dat)
lefMe=plot_ef_bar(lef_out_1)
lefMe
similarity_matrix <- dist.taxa(otu_table, method = "unifrac")


plot_cladogram(lef_out_1, color = c("red","blue"), clade_label_level = 4,   trim.names = TRUE)

phy_high_2<-subset_samples(psd5lf, Treatment %in% c("C2TF", "MAG"))

phy_high_2
table(sample_data(phy_high_2)$Treatment)



lef_out_2<-run_lefse(phy_high_2, group = "Treatment", norm = "CPM",
                     taxa_rank = "Genus",
                     
                     kw_cutoff = 0.05, lda_cutoff = 3)

lef_out_2

library(gridExtra)
grid.arrange(lef_out,lef_out_1,lef_out_2)




dat<-marker_table(lef_out_2) %>% data.frame() %>% select(1:5)

head(dat)
write.table(dat, file="lef_MAG_C2TF.txt", sep=",")

dat %>% filter(enrich_group=="ME") %>% head()
dat %>% kable(align = "c")
view(dat)
lefMag=plot_ef_bar(lef_out_2)
lefMag


pht <- run_posthoc_test(phy_high_1, group = "ME")
pht
##########
phy_high_3<-subset_samples(psd5lf, Treatment %in% c("MAG","BIO"))

phy_high_3
table(sample_data(phy_high_3)$Treatment)



lef_out_3<-run_lefse(phy_high_3, group = "Treatment", norm = "CPM",
                     taxa_rank = "Genus",
                     
                     kw_cutoff = 0.05, lda_cutoff = 3)

lef_out_3


dat<-marker_table(lef_out_3) %>% data.frame() %>% select(1:5)

head(dat)
write.table(dat, file="lef_BIO_MAG.txt", sep=",")

dat %>% filter(enrich_group=="ME") %>% head()
dat %>% kable(align = "c")
view(dat)
lef_ma_bio=plot_ef_bar(lef_out_3)
lef_ma_bio

phy_high_4<-subset_samples(psd5lf, Treatment %in% c("MAG","ME"))

phy_high_4
table(sample_data(phy_high_4)$Treatment)
lef_out_4<-run_lefse(phy_high_4, group = "Treatment", norm = "CPM",
                    taxa_rank = "Genus",
                    
                    kw_cutoff = 0.05, lda_cutoff = 3)

lef_out_4
library(tidyverse)
library(knitr)

dat<-marker_table(lef_out_4) %>% data.frame() %>% select(1:5)

head(dat)
write.table(dat, file="lef_ME_MAG.txt", sep=",")

dat %>% filter(enrich_group=="ME") %>% head()
dat %>% kable(align = "c")
view(dat)
lef_me_ma=plot_ef_bar(lef_out_4)

lef_me_ma



phy_high_6<-subset_samples(psd5lf, Treatment %in% c("ME","BIO", "MAG","C2TF"))

phy_high_6
table(sample_data(phy_high_6)$Treatment)



lef_out_6<-run_lefse(phy_high_6, group = "Treatment", norm = "CPM",
                     taxa_rank = "Family",
                     kw_cutoff = 0.05,
                     lda_cutoff = 3)

lef_out_6
library(tidyverse)
library(knitr)

dat<-marker_table(lef_out_6) %>% data.frame() %>% select(1:5)
write.table(dat, file="lef_Treat_C2TF.txt", sep=",")
head(dat)

dat %>% filter(enrich_group=="MAG") %>% head()
dat %>% kable(align = "c")
View(dat)
plot_ef_bar(lef_out_6)
p_abd <- plot_abundance(lef_out_6, group = "Treatment")
p_abd
plot_ef_dot(lef_out_6)
library(microeco)
# clade_label_level 5 represent phylum level in this analysis
# require ggtree package
plot_diff_cladogram(lef_out_6, use_taxa_num = 200, use_feature_num = 50, clade_label_level = 3, group_order = c("ME","BIO", "MAG","C2TF"))
#################
plot_heatmap(lef_out_6, transform = "log10p", group = "Treatment", color = AcidPlots::YelloBlue)

p_pht <- plot_postHocTest(lef_out_6, feature = "p__Bacteroidetes|g__Bacteroides")
p_pht

plot_cladogram(lef_out_6, color = c(BIO = "darkgreen", ME = "red", MAG = "pink", C2TF = "blue")) +
  theme(plot.margin = margin(0, 0, 0, 0))

plot_cladogram(lef_out_6, color = c("red","blue","green", "grey"), clade_label_level = 4)

phy_high_7<-subset_samples(psd5lf, Treatment %in% c("ME","BIO", "MAG"))

phy_high_7
table(sample_data(phy_high_7)$Treatment)



lef_out_7<-run_lefse(phy_high_7, group = "Treatment", norm = "CPM",
                     taxa_rank = "Genus",
                     kw_cutoff = 0.05, lda_cutoff = 3)

lef_out_7
library(tidyverse)
library(knitr)

dat<-marker_table(lef_out_7) %>% data.frame() %>% select(1:5)
write.table(dat, file="lef_treat.txt", sep=",")
head(dat)

dat %>% filter(enrich_group=="MAG") %>% head()
dat %>% kable(align = "c")
View(dat)
plot_ef_bar(lef_out_7)


phy_high_5<-subset_samples(psd5lf, Treatment %in% c("ME","BIO"))

phy_high_5
table(sample_data(phy_high_5)$Treatment)



lef_out_5<-run_lefse(phy_high_5, group = "Treatment", norm = "CPM",
                     taxa_rank = "Genus",
                     
                     kw_cutoff = 0.05, lda_cutoff = 3)

lef_out_5
library(tidyverse)
library(knitr)

dat<-marker_table(lef_out_5) %>% data.frame() %>% select(1:5)

head(dat)
write.table(dat, file="lef_BIO_ME.txt", sep=",")

dat %>% filter(enrich_group=="ME") %>% head()
dat %>% kable(align = "c")
view(dat)
lef_me_bio=plot_ef_bar(lef_out_5)

p_abd <- plot_abundance(lef_out_5, group = "Treatment")
p_abd

phy_high_8<-subset_samples(psd5lf, Treatment %in% c("C2T0","C2TF"))

phy_high_8
table(sample_data(phy_high_8)$Treatment)



lef_out_8<-run_lefse(phy_high_8, group = "Treatment", norm = "CPM",
                     taxa_rank = "Phylum",
                     
                     kw_cutoff = 0.05, lda_cutoff = 3)

lef_out_8
library(tidyverse)
library(knitr)

dat<-marker_table(lef_out_8) %>% data.frame() %>% select(1:5)

head(dat)
write.table(dat, file="lef_C1_C2_TF.txt", sep=",")

dat %>% filter(enrich_group=="ME") %>% head()
dat %>% kable(align = "c")
view(dat)
plot_ef_bar(lef_out_8)


phy_high_9<-subset_samples(psd5lf, Treatment %in% c("C1T0","C2TF"))

phy_high_9
table(sample_data(phy_high_9)$Treatment)



lef_out_9<-run_lefse(phy_high_9, group = "Treatment", norm = "CPM",
                     taxa_rank = "Phylum",
                     
                     kw_cutoff = 0.05, lda_cutoff = 3)

lef_out_9
library(tidyverse)
library(knitr)

dat1<-marker_table(lef_out_9) %>% data.frame() %>% select(1:5)

head(dat1)
write.table(dat1, file="lef_C_0_f.txt", sep=",")

dat %>% filter(enrich_group=="ME") %>% head()
dat %>% kable(align = "c")

view(dat1)
lefcontrol=plot_ef_bar(lef_out_9)
lefcontrol
phy_high_5<-subset_samples(psd5lf, Treatment %in% c("ME","BIO"))

phy_high_5
table(sample_data(phy_high_5)$Treatment)



lef_out_5<-run_lefse(phy_high_5, group = "Treatment", norm = "CPM",
                     taxa_rank = "Genus",
                     
                     kw_cutoff = 0.05, lda_cutoff = 3)

lef_out_5
library(tidyverse)
library(knitr)

dat<-marker_table(lef_out_5) %>% data.frame() %>% select(1:5)

head(dat)
write.table(dat, file="lef_BIO_ME.txt", sep=",")

dat %>% filter(enrich_group=="ME") %>% head()
dat %>% kable(align = "c")
view(dat)
plot_ef_bar(lef_out_5)

library("data.table")
library("phyloseq")
library("ggplot2")
library("grid")
library("tidyverse")
Sys.setlocale("LC_COLLATE", "C")
tax1 = lefse_1name_obj(psd5, sample_data(ps5)$Treatment)
lefse1 = lefse_obj(ps1)
lefse1 = rbind(tax1, lefse1)

pss_3 <-phy_high %>%
  tax_glom("Phylum")
phyloseq::taxa_names(pss_3) <- phyloseq::tax_table(pss_3)[, "Phylum"]
phyloseq::psmelt(pss_3) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "Time", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")
  


pss_3 <-phy_high_1 %>%
  tax_glom("Phylum")
phyloseq::taxa_names(pss_3) <- phyloseq::tax_table(pss_3)[, "Phylum"]
phyloseq::psmelt(pss_3) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "Time", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")


pss_3 <-phy_high_2 %>%
  tax_glom("Phylum")
phyloseq::taxa_names(pss_3) <- phyloseq::tax_table(pss_3)[, "Phylum"]
phyloseq::psmelt(pss_3) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "Time", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")
 
pss_3 <-phy_high_3 %>%
  tax_glom("Phylum")
phyloseq::taxa_names(pss_3) <- phyloseq::tax_table(pss_3)[, "Phylum"]
phyloseq::psmelt(pss_3) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "Time", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

pss_3 <-phy_high_4 %>%
  tax_glom("Phylum")
phyloseq::taxa_names(pss_3) <- phyloseq::tax_table(pss_3)[, "Phylum"]
phyloseq::psmelt(pss_3) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "Time", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

pss_3 <-phy_high_5 %>%
  tax_glom("Phylum")
phyloseq::taxa_names(pss_3) <- phyloseq::tax_table(pss_3)[, "Phylum"]
phyloseq::psmelt(pss_3) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "Time", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")
####################################################################################
exo = subset_taxa(phy_high_6, Family=="Geobacteraceae")

phyloseq::taxa_names(pss_3) <- phyloseq::tax_table(pss_3)[, "Family"]
phyloseq::psmelt(exo) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Treatment), height = 0, width = .2) +
  labs(x = "Treatment", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

#####################################################################################
pss_3 <-phy_high_6 %>%
  tax_glom("Family")
phyloseq::taxa_names(pss_3) <- phyloseq::tax_table(pss_3)[, "Family"]
phyloseq::psmelt(pss_3) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "Treatment", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

library(dplyr)
library(tidyr)

write.table(psd5 %>% tax_glom(taxrank = "Phylum") %>% 
              transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
              select(Phylum, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "ps.relative_abundance_p.csv", sep = "\t", quote = T, row.names = F, col.names = T)

write.table(psd5 %>% tax_glom(taxrank = "Family") %>% 
              transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
              select(Phylum, Class, Family, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "ps.relative_abundance.family_41.csv", sep = "\t", quote = T, row.names = F, col.names = T)

write.table(psd5 %>% tax_glom(taxrank = "Genus") %>% 
              transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
              select(Phylum, Class, Family, Genus, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "ps.relative_abundance.genus.csv", sep = "\t", quote = T, row.names = F, col.names = T)

phy_high_6= transform_sample_counts(phy_high_6, function(x) x / sum(x) * 100) 

dat <- phy_high_6 %>% tax_glom(taxrank = "Family") %>% psmelt()
Me = dat %>% group_by(name, Treatment) %>% 
  summarise(Geobacteraceae = Abundance[Family == "Geobacteraceae"], 
            Marinilabiliaceae = Abundance[Family == "Marinilabiliaceae"], 
            Vicinamibacterales = Abundance[Family == "Vicinamibacteraceae"]) %>% 
           # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  scale_fill_brewer(palette="Purples")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)


phy_high_6= transform_sample_counts(phy_high_6, function(x) x / sum(x) * 100) 

dat <- phy_high_6 %>% tax_glom(taxrank = "Family") %>% psmelt()
mag = dat %>% group_by(name, Treatment) %>% 
  summarise(Hungateiclostridiaceae  = Abundance[Family == "Hungateiclostridiaceae"], 
            Hydrogenophilaceae = Abundance[Family == "Hydrogenophilaceae"], 
            Anaerolineaceae = Abundance[Family == "Anaerolineaceae"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  scale_fill_brewer(palette="Blues")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)


phy_high_6= transform_sample_counts(phy_high_6, function(x) x / sum(x) * 100) 

dat <- phy_high_6 %>% tax_glom(taxrank = "Family") %>% psmelt()
dat %>% group_by(name, Treatment) %>% 
  summarise(Methylobacter  = Abundance[Family == "Methylomonadaceae"], 
             Methylobacterium = Abundance[Family == "Beijerinckiaceae"], 
            Methylomicrobium = Abundance[Family == "Methylomonadaceae"],
            Methyloligellaceae = Abundance [Family == "Methyloligellaceae"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  scale_fill_brewer(palette="Blues")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)



newSTorder = c("C2TF", "BIO", "ME", "MAG")
a_my_comparisons <- list(c("ME", "MAG"), c("MAG", "BIO"), c("BIO", "ME"), c("BIO", "C2TF"), c("ME", "C2TF"), c("C2TF", "MAG"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))                        
                         
#
palette <- c("#20DE8B","#CCDE8B", "#FFDE8B", "#FFA88B", "#FF6A8B", "#AD8CFF", "#90CFFF")
a_my_comparisons <- list(c("C2T0", "C2TF"), c("ME", "MAG"), c("C2TF", "MAG", c("MAG", "BIO"), c("BIO", "ME")))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
p2 <- plot_richness(psd5, x="Treatment", measures=c("ACE", "Shannon", "Simpson","Evenness"))+ geom_boxplot(aes(fill = Treatment), alpha=.6)+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))
scale_fill_manual(values=palette)+
  
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)


lef = data.frame(sample_data(phy_high_6))

#TUKEY SHAnNON

anova.sh = aov(Treatment ~ , data=lef)
summary(anova.sh)
tukey = TukeyHSD(anova.sh)
tukey
write.table(TukeyHSD, file="TukeyHSD.txt", sep=",")

library(agricolae)
tukey.test2 <- HSD.test(anova.sh, "sample_data(mrg.epi)$Treatment")
with(mrg.epi, tukey.test2)
tukey.test2

newSTorder = c("C2TF", "C1T0")
a_my_comparisons <- list(c("C2TF", "C1T0"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

phy_high_9= transform_sample_counts(phy_high_9, function(x) x / sum(x) * 100) 

dat <- phy_high_9 %>% tax_glom(taxrank = "Family") %>% psmelt()
dat %>% group_by(name, Treatment) %>% 
  summarise(Hydrogenophilaceae= Abundance[Family == "Hydrogenophilaceae"], 
            Anaerolineaceae = Abundance[Family == "Anaerolineaceae"], 
            HungateiclostridiaceaeAnaerolinea = Abundance[Family == "Hungateiclostridiaceae"],
            Methylomonadaceae = Abundance[Family == "Methylomonadaceae"],
            Microscillaceae = Abundance[Family == "Microscillaceae"],
            Spirochaetaceae = Abundance[Family == "Spirochaetaceae"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  scale_fill_brewer(palette="Reds")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)


newSTorder = c("C2TF", "C1T0")
a_my_comparisons <- list(c("C2TF", "C1T0"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

phy_high_9= transform_sample_counts(phy_high_9, function(x) x / sum(x) * 100) 

dat <- phy_high_9 %>% tax_glom(taxrank = "Phylum") %>% psmelt()
ControlA = dat %>% group_by(name, Treatment) %>% 
  summarise(Firmicutes= Abundance[Phylum == "Firmicutes"], 
            Chloroflexi = Abundance[Phylum == "Chloroflexi"], 
            Bacteroidota = Abundance[Phylum == "Bacteroidota"],
            Proteobacteria = Abundance[Phylum == "Proteobacteria"],
            Desulfobacterota = Abundance[Phylum == "Desulfobacterota"],
            Acidobacteriota = Abundance[Phylum == "Acidobacteriota"],
            Myxococcota = Abundance[Phylum == "Myxococcota"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Phylum", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Phylum, nrow = 1) +
  labs(y="Relative abundance")+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  theme(strip.text = element_text(face = "italic", hjust = 0.5, color = "black", size = 8.5))+
  scale_fill_brewer(palette="")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title=element_text(size=12,face="bold"))+    
  theme(axis.text = element_text(size = 10, face="bold"))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
ControlA
A=grid.arrange(lefcontrol,ControlA)


my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )


newSTorder = c("C2TF", "BIO", "ME", "MAG")
a_my_comparisons <- list(c("ME", "MAG"), c("MAG", "BIO"), c("BIO", "ME"), c("BIO", "C2TF"), c("ME", "C2TF"), c("C2TF", "MAG"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")) 

phy_high_6= transform_sample_counts(phy_high_6, function(x) x / sum(x) * 100) 
library(devtools)
install_github("umerijaz/microbiomeSeq") 
library(microbiomeSeq)
kw_sig <- kruskal_expression(phy_high_6, grouping_column = "Treatment")


dat <- phy_high_6 %>% tax_glom(taxrank = "Family") %>% psmelt()
control = dat %>% group_by(name, Treatment) %>% 
  summarise(Microscillaceae  = Abundance[Family == "Microscillaceae"], 
            Sandaracinaceae = Abundance[Family == "Sandaracinaceae"], 
            Desulfosarcinaceae = Abundance[Family == "Desulfosarcinaceae"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  scale_fill_brewer(palette="Greens")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  stat_compare_means(method = "kruskal.test", comparisons = my_comparisons)

  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
control



dat <- phy_high_6 %>% tax_glom(taxrank = "Family") %>% psmelt()
bio = dat %>% group_by(name, Treatment) %>% 
  summarise(Clostridiaceae = Abundance[Family == "Clostridiaceae"], 
            Christensenellaceae = Abundance[Family == "Christensenellaceae"], 
            BIrii41 = Abundance[Family == "BIrii41"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  scale_fill_brewer(palette="Reds")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

library(gridExtra)
grid.arrange(Me, mag, control, bio)
###################################################################################################
#duo
set.seed(67722)
library(ggplot2)
newSTorder = c("BIO", "C2TF")
a_my_comparisons <- list(c("BIO", "C2TF"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")) 
phy_high= transform_sample_counts(phy_high, function(x) x / sum(x) * 100)
dat <- phy_high %>% tax_glom(taxrank = "Order") %>% psmelt()
bio = dat %>% group_by(name, Treatment) %>% 
  summarise(Proteinivoracales = Abundance[Order == "Proteinivoracales"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  labs(y="Relative abundance")+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  theme(strip.text = element_text(face = "italic", hjust = 0.5, color = "black", size = 12))+
  scale_fill_brewer(palette="BrBG")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title=element_text(size=12,face="bold"))+    
  theme(axis.text = element_text(size = 12, face="bold"))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
bio

dat <- phy_high %>% tax_glom(taxrank = "Genus") %>% psmelt()
bio_1 = dat %>% group_by(name, Treatment) %>% 
  summarise(Sporacetigenium = Abundance[Genus == "Sporacetigenium"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  scale_fill_brewer(palette="BrBG")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
bio_1

dat <- phy_high %>% tax_glom(taxrank = "Family") %>% psmelt()
bio_2 = dat %>% group_by(name, Treatment) %>% 
  summarise(Alkaliphilus = Abundance[Family == "Alkaliphilus"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  labs(y="Relative abundance")+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  theme(strip.text = element_text(face = "italic", hjust = 0.5, color = "black", size = 12))+
  scale_fill_brewer(palette="BrBG")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title=element_text(size=12,face="bold"))+    
  theme(axis.text = element_text(size = 12, face="bold"))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
bio_2

bio = grid.arrange(bio, bio_2, ncol=2)

#######
a=grid.arrange(lefBIO, bio, ncol=2)
a=grid.arrange(lefBIO, bio)
###################################################################################################
#duo
newSTorder = c("MAG", "C2TF")
a_my_comparisons <- list(c("MAG", "C2TF"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")) 
#phy_high_2= transform_sample_counts(phy_high_2, function(x) x / sum(x) * 100)

dat <- phy_high_2 %>% tax_glom(taxrank = "Genus") %>% psmelt()
mag_1 = dat %>% group_by(name, Treatment) %>% 
  summarise(Thiobacillus = Abundance[Genus == "Thiobacillus"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  labs(y="Relative abundance")+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  theme(strip.text = element_text(face = "italic", hjust = 0.5, color = "black", size = 12))+
  scale_fill_brewer(palette="Reds")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title=element_text(size=12,face="bold"))+    
  theme(axis.text = element_text(size = 12, face="bold"))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
mag_1

mag_2 = dat %>% group_by(name, Treatment) %>% 
  summarise(Treponema = Abundance[Genus == "Treponema"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  labs(y="Relative abundance")+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  theme(strip.text = element_text(face = "italic", hjust = 0.5, color = "black", size = 12))+
  scale_fill_brewer(palette="Reds")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title=element_text(size=12,face="bold"))+    
  theme(axis.text = element_text(size = 12, face="bold"))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
mag_2

mag = grid.arrange(mag_1, mag_2, ncol=2)
b=grid.arrange(lefMag, mag_1, ncol=2)
b=grid.arrange(lefMag, mag)
b
############################################################
set.seed(9584)
newSTorder = c("ME", "C2TF")
a_my_comparisons <- list(c("ME", "C2TF"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")) 
#phy_high_1= transform_sample_counts(phy_high_1, function(x) x / sum(x) * 100)
phy_high_1
dat <- phy_high_1 %>% tax_glom(taxrank = "Family") %>% psmelt()
Mel = dat %>% group_by(name, Treatment) %>% 
  summarise(Geobacteraceae = Abundance[Family == "Geobacteraceae"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  labs(y="Relative abundance")+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  theme(strip.text = element_text(face = "italic", hjust = 0.5, color = "black", size = 12))+
  scale_fill_brewer(palette="Greens")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title=element_text(size=12,face="bold"))+    
  theme(axis.text = element_text(size = 12, face="bold"))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
Mel

dat <- phy_high_1 %>% tax_glom(taxrank = "Genus") %>% psmelt()
Mel_1 = dat %>% group_by(name, Treatment) %>% 
  summarise(Marmoricola = Abundance[Genus == "Marmoricola"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  labs(y="Relative abundance")+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  theme(strip.text = element_text(face = "italic", hjust = 0.5, color = "black", size = 12))+
  scale_fill_brewer(palette="Greens")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title=element_text(size=12,face="bold"))+    
  theme(axis.text = element_text(size = 12, face="bold"))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
Mel_1

mel = grid.arrange(Mel, Mel_1, ncol=2)
d=grid.arrange(lefMe, mel, ncol=2)
d=grid.arrange(lefMe, mel)

grid.arrange(a, b, d)
grid.arrange(lefBIO, bio, lefMe, Mel_1, lefMag, mag_1)
grid.arrange(lefBIO, lefMe,lefMag)
###############################################################

newSTorder = c("ME", "MAG")
a_my_comparisons <- list(c("ME", "MAG"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")) 
phy_high_4= transform_sample_counts(phy_high_4, function(x) x / sum(x) * 100)

dat <- phy_high_4 %>% tax_glom(taxrank = "Genus") %>% psmelt()
mel_mag = dat %>% group_by(name, Treatment) %>% 
  summarise(Cytophaga_xylanolytica = Abundance[Genus == "[Cytophaga] xylanolytica group"],
            Thiobacillus= Abundance[Genus == "Thiobacillus"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Genus", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Genus, nrow = 1) +
  labs(y="Relative abundance")+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  theme(strip.text = element_text(face = "italic", hjust = 0.5, color = "black", size = 8.5))+
  scale_fill_brewer(palette="YlGn")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title=element_text(size=12,face="bold"))+    
  theme(axis.text = element_text(size = 10, face="bold"))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
mel_mag

newSTorder = c("BIO", "MAG")
a_my_comparisons <- list(c("BIO", "MAG"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")) 
phy_high_3= transform_sample_counts(phy_high_3, function(x) x / sum(x) * 100)

dat <- phy_high_3 %>% tax_glom(taxrank = "Genus") %>% psmelt()
bio_mag = dat %>% group_by(name, Treatment) %>% 
  summarise(Clostridium = Abundance[Genus == "Clostridium sensu stricto 1"],
            Thiobacillus= Abundance[Genus == "Thiobacillus"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  labs(y="Relative abundance")+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  theme(strip.text = element_text(face = "italic", hjust = 0.5, color = "black", size = 8.5))+
  scale_fill_brewer(palette="Oranges")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title=element_text(size=12,face="bold"))+    
  theme(axis.text = element_text(size = 10, face="bold"))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
bio_mag

newSTorder = c("BIO", "ME")
a_my_comparisons <- list(c("BIO", "ME"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")) 
phy_high_5= transform_sample_counts(phy_high_5, function(x) x / sum(x) * 100)

dat <- phy_high_5%>% tax_glom(taxrank = "Family") %>% psmelt()
bio_mel = dat %>% group_by(name, Treatment) %>% 
  summarise(Anaerolineaceae = Abundance[Family == "Anaerolineaceae"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Family", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Family, nrow = 1) +
  labs(y="Relative abundance")+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  theme(strip.text = element_text(face = "italic", hjust = 0.5, color = "black", size = 8.5))+
  scale_fill_brewer(palette="BuPu")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title=element_text(size=12,face="bold"))+    
  theme(axis.text = element_text(size = 10, face="bold"))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
bio_mel

dat <- phy_high_5%>% tax_glom(taxrank = "Genus") %>% psmelt()
bio_mel_1 = dat %>% group_by(name, Treatment) %>% 
  summarise(Cytophaga_xylanolytica = Abundance[Genus == "[Cytophaga] xylanolytica group"]) %>% 
  # Other = sum(Abundance[!Family %in% c("Geobacteraceae","Pseudopelobacteraceae","Clostridiaceae")])) %>%
  pivot_longer(cols = -c(name, Treatment), names_to = "Genus", values_to = "Abundance") %>% 
  ggplot(aes(Treatment, Abundance, fill = Treatment)) + geom_boxplot() + facet_wrap(~ Genus, nrow = 1) +
  labs(y="Relative abundance")+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  theme_light()+
  theme(strip.text = element_text(face = "italic", hjust = 0.5, color = "black", size = 8.5))+
  scale_fill_brewer(palette="BuPu")+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title=element_text(size=12,face="bold"))+    
  theme(axis.text = element_text(size = 10, face="bold"))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

bio_mel_1
b_m=grid.arrange(bio_mel,bio_mel_1, ncol=2)

grid.arrange(mel_mag, bio_mag, b_m, ncol=3)


grid.arrange(a,b,c)
