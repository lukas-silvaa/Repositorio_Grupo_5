# Repositorio_Grupo_5


#Se descargaron los archivos .txt "tax_table", "otu_table" y "df_table" (ubicados en github javi) y se abrieron en R como data.frames
#Se descargó el archivo "ps_arms.R" (ubicado en github javi). 

#librerias 
library(phyloseq) 
library (ggplot)
library (decontam)

#crear funciones utilizadas en phyloseq

otu <- function(file_otu) {
  otus <- as.matrix(file_otu)
  otus_format <- otu_table(otus, taxa_are_rows = TRUE)
  return(otus_format)
}

taxa <- function(file_taxa) {
  taxas <- (as.matrix(file_taxa))
  taxa_format <- (tax_table(taxas))
  return(taxa_format)
}

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

