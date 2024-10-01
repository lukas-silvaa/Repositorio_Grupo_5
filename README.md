# Repositorio_Grupo_5

#Se descargaron los archivos .txt "tax_table", "otu_table" y "df_table" (ubicados en github javi) y se abrieron en R como data.frames
#Se descargó el archivo "ps_arms.R" (ubicado en github javi). 

#librerias 
library(phyloseq) 
library (ggplot)
library (decontam)

#crear funciones utilizadas en phyloseq
#file_otu (en las filas los ASVs y en las columnas las muestras)

otu <- function(file_otu) {
  otus <- as.matrix(file_otu)
  otus_format <- otu_table(otus, taxa_are_rows = TRUE)
  return(otus_format)
}

#Crear un formato de tabla de taxonomía compatible con phyloseq
#file_taxa (en las filas los ASV y en las columnas los niveles taxonómicos)

taxa <- function(file_taxa) {
  taxas <- (as.matrix(file_taxa))
  taxa_format <- (tax_table(taxas))
  return(taxa_format)
}

#Crear un formato de datos de muestra compatible con phyloseq
#file_data (en las filas las muestras y en las columnas las variables)

data <- function(file_df) {
  df <- data.frame(file_df)
  df_format <- sample_data(df)
  return(df_format)
}

#PASO 1: Crear objeto phyloseq -> para ordenar datos de metabarcoding en un solo espacio

setwd("Poner ubicación donde se encuentran los archivos .txt")
project_path <- "Poner ubicación donde se encuentran los archivos .txt"
project_path
#abrir los archivos (también se pueden abrir vía import dataset en el enviroment -> Heading "yes")
asv_1 <- read.table("otu_table.txt", header = TRUE, sep = "\t", row.names = 1)
tax_1 <- read.table("tax_table.txt", header = TRUE, sep = "\t", row.names = 1)
tax_1[tax_1 == ""] <- "Unassigned"
df <- read.table("df_table.txt", header = TRUE, sep = "\t", row.names = 1)

otu_ps <- otu(asv_1)
tax_ps <- taxa(tax_1)
data_ps <- data(df)
ps <- phyloseq(otu_ps, tax_ps, data_ps)
ps <- subset_samples(ps, study == "arms")

#PASO 2: Descontaminar -> Remover contaminantes basados en control negativo  

df <- as.data.frame(sample_data(ps))
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

#Visualizar 

contam <- ggplot(data = df, aes(x = sample_name, y = LibrarySize,
        color = factor)) + geom_point() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 6))
contam
#ggsave("plot_papper/contam.png", contam, width = 6, height = 3) 

#identificación de controles negativos
sample_data(ps)$is.neg <- sample_data(ps)$type == "control"

#detección de contaminantes
contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
filter(contamdf.prev05, contaminant == TRUE)

#selection seq contaminants
contam <- head(which(contamdf.prev05$contaminant))
contam <- rownames(contamdf.prev05[c(contam),])

#remov contaminants
otu_table(ps) <- otu_table(ps)[!(row.names(otu_table(ps)) %in% contam),]
ps <- subset_samples(ps, type != "control")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

#AHORA EN  R
#Limpiar datos de especies que no son invertebrados marinos --> Data clean
ps = subset_taxa(ps, Specie !="Eurytemora foveola")#AGUA DULCE
ps = subset_taxa(ps, Specie !="Eurytemora herdmani")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Fridericia")
ps = subset_taxa(ps, Genus !="Dussartcyclops")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Neoergasilus")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Sinocalanus")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Bryodrilus")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Rhysida")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Potamothrix")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Baikalodrilus")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Fridericia")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Lumbriculus")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Lamprodrilus")#AGUA DULCE
ps = subset_taxa(ps, Family !="Diaptomidae")
ps = subset_taxa(ps, Family !="Paramelitidae")
ps = subset_taxa(ps, Family !="Potamidae")
ps = subset_taxa(ps, Family !="Ampullariidae")
ps = subset_taxa(ps, Family !="Orobdellidae")
ps = subset_taxa(ps, Family !="Trombiculidae")
ps = subset_taxa(ps, Family !="Buthidae")
ps = subset_taxa(ps, Family !="Phytoseiidae")
ps = subset_taxa(ps, Family !="Unionicolidae") 
ps = subset_taxa(ps, Order !="Crassiclitellata") #gusano terrestre
ps = subset_taxa(ps, Order !="Lumbriculida")#AGUA DULCE
ps = subset_taxa(ps, Order !="Anura")
ps = subset_taxa(ps, Order !="Artiodactyla") #angulados
ps = subset_taxa(ps, Order !="Chiroptera") #murcielago
ps = subset_taxa(ps, Order !="Eulipotyphla")
ps = subset_taxa(ps, Order !="Physariida")
ps = subset_taxa(ps, Order !="Mesostigmata")
ps = subset_taxa(ps, Order !="Araneae")
ps = subset_taxa(ps, Order !="Sarcoptiformes")
ps = subset_taxa(ps, Order !="Ixodida")
ps = subset_taxa(ps, Order !="Opiliones")
ps = subset_taxa(ps, Order !="Squamata")
ps = subset_taxa(ps, Class !="Eumycetozoa") #nolint
ps = subset_taxa(ps, Class !="Insecta") #nolint
ps = subset_taxa(ps, Class !="Collembola")
ps = subset_taxa(ps, Class !="Diplopoda")
ps = subset_taxa(ps, Class !="Actinopteri")
ps = subset_taxa(ps, Phylum !="Ascomycota")
ps = subset_taxa(ps, Phylum !="Basidiomycota")
ps = subset_taxa(ps, Phylum !="Mucoromycota")
ps = subset_taxa(ps, Phylum !="Chordata")
ps = subset_taxa(ps, Phylum !="Chlorophyta")
ps = subset_taxa(ps, Phylum !="Rhodophyta")
ps = subset_taxa(ps, Phylum !="Discosea")
ps = subset_taxa(ps, Phylum !="Bacillariophyta")
ps = subset_taxa(ps, Kingdom !="Bacteria")

#Contar número de secuencias, número de taxas 

ps <- prune_taxa(taxa_sums(ps) > 0,
                 ps)

count_total <- ntaxa(ps)

count_know <- ps %>%
  subset_taxa(Phylum != "Unassigned")  %>% 
  subset_taxa(Phylum != "NA")  %>% 
  ntaxa()

100 - ((count_know*100)/count_total)

#luego de la limpieza de datos se obtuvieron 10676 ASVs


#Borrar los primers con 0 lecturas

clean_zero_reads <- function(ps, taxonomic_level) {
  tax_glom_ps <- tax_glom(ps, taxonomic_level, NArm = FALSE) #agrupa por taxonomic leves eg. genus #nolint
  clean_zero_ps <- prune_taxa(taxa_sums(tax_glom_ps) > 0,
                              tax_glom_ps) #remove taxonommic_level=0
  return(clean_zero_ps)
}
ps <- clean_zero_reads(ps, "Specie")


#MERGE, juntar sub-replicas según profundidad

ps_sfs <- merge_samples(ps, group = "factor_r")
df_sfs <- data.frame(sample_data(ps_sfs))
df_sfs <- df_sfs  %>%
  mutate(study = rep("arms", nrow(.)),
         site = rep(c("Algarrobo", 8),
         ocean_depth= rep(c("30m", "60m"), each = 2, 2),
         factor = paste(site, ocean_depth, sep = "_"),
         factor_r = paste(site, ocean_depth, replica, sep = "_")) %>%
  dplyr::select(study, site, ocean_depth, factor, factor_r)
sample_data(ps_sfs) <- df_sfs

ps_sfs <- clean_zero_reads(ps_sfs, "Specie")


#Rarefacción: Hacer un submuestreo debido a que todas las muestras tienen un número de distinto de reads.

rarefaction <- function(clean_zero_ps) {
  sample_sums <- sample_sums(clean_zero_ps) #nolint
  min_reads_sample <- which.min(sample_sums) #nolint
  min_reads <- sample_sums[min_reads_sample] #nolint
  min_reads <- as.numeric(min_reads) #nolint
  rarefied_ps <- rarefy_even_depth(clean_zero_ps,
                                   sample.size = min_reads,
                                   replace = TRUE)
  return(rarefied_ps)
}

