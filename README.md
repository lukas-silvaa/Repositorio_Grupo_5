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


#Abrir los archivos (también se pueden abrir vía import dataset en el enviroment -> Heading "yes")--------------------

asv_1 <- read.table("otu_table.txt", header = TRUE, sep = "\t", row.names = 1)
tax_1 <- read.table("tax_table.txt", header = TRUE, sep = "\t", row.names = 1)
tax_1[tax_1 == ""] <- "Unassigned"
df <- read.table("df_table.txt", header = TRUE, sep = "\t", row.names = 1)

df <- df %>% 
mutate(factor =  paste(site, ocean_depth, sep = "_")) %>%
mutate(factor_r =  paste(site, ocean_depth, replica, sep = "_"))


#Crear objeto phyloseq----------------------------------------

otu <- as.matrix(asv_1)
otu_ps <- otu_table(otu, taxa_are_rows = TRUE)

taxa <- as.matrix(tax_1)
tax_ps <- tax_table(taxa)

df <- data.frame(df)
data_ps <- sample_data(df)

ps <- phyloseq(otu_ps, tax_ps, data_ps)
ps <- subset_samples(ps, study == "arms")


#PASO 2: Descontaminar -> Remover contaminantes basados en control negativo  

#Preparación de data
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

#seleccción de "seq contaminants"
contam <- head(which(contamdf.prev05$contaminant))
contam <- rownames(contamdf.prev05[c(contam),])

#remover contaminantes
otu_table(ps) <- otu_table(ps)[!(row.names(otu_table(ps)) %in% contam),]
ps <- subset_samples(ps, type != "control")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

#AHORA EN  R---------------------------------------------------------------
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


#MERGE, juntar sub-replicas según profundidad --------------------------------------

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

otu_table(ps_sfs) <- t(otu_table(ps_sfs))

ps_sfs <- clean_zero_reads(ps_sfs, "Specie")


#RAREFACCIÓN: Hacer un submuestreo debido a que todas las muestras tienen un número de distinto de reads-----------------------

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

sample_sums <- sample_sums(ps_sfs) #nolint
which.min(sample_sums)
r_ps_sfs <- rarefaction(ps_sfs)

##CURVAS DE RAREFACCION--------------------------------

calculate_rarefaction_curves2 <- function(psdata, measures, depths) {
  require("plyr") # ldply
  require("reshape2") # melt
  
estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
  rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
  alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
#Calcular equitabilidad (evenness)
evenness_values <- diversity(t(rarified_psdata@otu_table)) / log(specnumber(t(rarified_psdata@otu_table)))

#Añadir la equitabilidad a los resultados
evenness_df <- data.frame(Sample = rownames(alpha_diversity), Evenness = evenness_values)
alpha_diversity <- cbind(alpha_diversity, Evenness = evenness_values)
    
#Convertir los resultados a formato largo (melted)
molten_alpha_diversity <- melt(as.matrix(alpha_diversity),
                                   varnames = c('Sample', 'Measure'),
                                   value.name = 'Alpha_diversity')
    
  molten_alpha_diversity
  }
  
names(depths) <- depths # Esto habilita la adición automática de la profundidad al output por ldply
rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
#Convertir Depth de factor a numérico
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

#Resumen de los resultados obtenidos en las curvas de rarefacción
curve_summary_verbose <- function(rarefied_ps, r_curve_data) {
  r_curve_data_summary <- ddply(r_curve_data,
                                c("Depth", "Sample", "Measure"),
                                summarise,
                                Alpha_diversity_mean = mean(Alpha_diversity), #nolint
                                Alpha_diversity_sd = sd(Alpha_diversity)) #nolint
  r_curve_data_summary_verbose <- merge(r_curve_data_summary,
                                        data.frame(sample_data(rarefied_ps)),
                                        by.x = "Sample", by.y = "row.names")
  return(r_curve_data_summary_verbose)
}


r_curve2 <- calculate_rarefaction_curves2(r_ps_sfs,
                                        c("Observed", "Shannon", "Chao1", "Evenness"),
                                        rep(c(1:150 * 100), each = 5))

#Data summary curve rarefy
r_curve_summary <- curve_summary_verbose(r_ps_sfs, r_curve2)


##plot y:reads x:sample, agglomerated by taxonomic level----
a <- ggplot(data = data.frame(x = 1:length(sample_sums(ps_sfs)), #nolint
                              y = sort(sample_sums(ps_sfs), decreasing = TRUE)), #nolint
            aes(x = x, y = y)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(title = " Reads for Sample",
       x = "Samples",
       y = "Reads")
a
#plot whit sample rarefy
b <- ggplot(data = data.frame(x = 1:length(sample_sums(r_ps_sfs)), #nolint
                              y = sort(sample_sums(r_ps_sfs), decreasing = TRUE)),
            aes(x = x, y = y)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(title = "Rarefy",
       x = "Samples",
       y = "Reads")
b
plot_rarefy <- a | b
#ggsave ("Plot_papper/reads_rarefy.png", plot_rarefy)

#plot whit alpha diversity obs and shannon
c <- r_curve_summary  %>% 
  filter(Measure != "se.chao1")  %>% 
  ggplot( aes(
              x = Depth, #nolint
              y = Alpha_diversity_mean, #nolint
              ymin = Alpha_diversity_mean - Alpha_diversity_sd, #nolint
              ymax = Alpha_diversity_mean + Alpha_diversity_sd, #nolint
              colour = as.factor(factor), #nolint
              group = Sample)) + #nolint
  geom_line(linewidth = 0.2, linetype = "solid") +
  facet_wrap(facets = ~ Measure, scales = "free_y") +
  labs(x = "Reads",
       y = "Alpha diversity mean") +
  scale_colour_manual(values = c("#17becf", "#8263e4"), 
                    labels = c("Algarrobo 30m", "Algarrobo 60m") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5), 
        legend.title = element_blank())
c
#ggsave("plot_papper/rarefaction_curve.png", c, width = 8, height = 4)
