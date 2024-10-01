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

#PASO 3: visualizar 

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

#remove contaminants
otu_table(ps) <- otu_table(ps)[!(row.names(otu_table(ps)) %in% contam),]
ps <- subset_samples(ps, type != "control")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
