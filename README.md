# Repositorio_Grupo_5
```ruby
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





#Barplot de Alfa diversidad--------------------------------------------------------------------

#ASV POR PHYLUM------------------
##contar numero de assigned by tax level-------
r_ps_sfs %>%
  clean_zero_reads(.,"Phylum" ) %>%
  subset_taxa(.,      Phylum != "Unassigned") %>%
  subset_taxa(.,      Phylum != "NA") %>%
  #subset_samples(., factor == "Algarrobo_30m") %>%
  ntaxa()

##contar numero de asv por phylum-----------
count_asv_phylum <- list()

r_ps_sfs_know <- subset_taxa(r_ps_sfs, Phylum != "Unassigned")
tax_phylum <- unique(tax_table(r_ps_sfs_know)[, "Phylum"])

for( i in tax_phylum) {
  resultado <- r_ps_sfs %>%
  subset_taxa(Phylum == i)  %>%
  prune_taxa(taxa_sums(.) > 0,
             .) %>%
  ntaxa()

 count_asv_phylum[[i]] <- resultado

 }

print(count_asv_phylum)

#ELEGIR PS A OCUPAR-------------------
r_ps_oc <- r_ps_sfs


#ABUNDNACIA POR PHYLUM------------------------
phylum_colors <- c("Annelida" = "#9467bd",# Lavanda
                   "Arthropoda" = "#ff7f0e",  # Naranja
                   "Bryozoa" = "#2ca02c",  # Verde
                   "Cnidaria" = "#bcbd22",  # Verde lima
                   "Echinodermata" = "#c5b0d5",  # Morado
                   "Mollusca" = "#d62728",  # Rojo
                   "Nematoda" = "#17becf",  # Cian
                   "Nemertea" = "#f7b6d2",  # Rosa claro
                   "Porifera" = "#1f77b4",  # Azul
                   "Phatyhelminthes" = "#62162f",  # Marrón
                   "Rotifera" = "#e377c2", 
                   "Chordata" = "#8c564b",
                   "Unassigned" = "#7f7f7f")  # Gris

###BAR PLOT ABSOLUTE RELATIVE#
taxonomic_level <- "Phylum"
sfs_rel <- r_ps_oc %>%
  tax_glom(taxrank = taxonomic_level) %>% 
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt() %>%                         # Melt to long format
  arrange(Phylum)

str(sfs_rel)
sfs_rel %>%
  group_by(factor) %>%
  summarise(total_abundance = sum(Abundance))

sfs_ps <- sfs_rel %>%
  dplyr::select(Phylum, Sample, Abundance,
         ocean_depth, site) %>%
  group_by(Phylum, ocean_depth, site) %>%
  filter(Abundance > 0.01)

bar_plot <- ggplot(sfs_ps) +
  geom_col(mapping = aes(
                         x = as.factor(ocean_depth),
                         y = Abundance,
                         fill = !!sym(taxonomic_level)),
           position = "fill",
           show.legend = TRUE) +
  theme_bw() +
  labs(x = "Ocean Depth (m)", y = "Proportion of Community") +
  facet_grid( ~ site, 
              labeller = as_labeller(c(`Algarrobo`= "Algarrobo", 
                                       `Las_Cruces` = "Las Cruces"))) +
  theme(axis.title = element_text(size = 13),  # Tamaño de los títulos de los ejes
        axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, hjust = 0.5),   # Tamaño de las etiquetas del eje X
        axis.text.y = element_text(size = 9), # Tamaño de las etiquetas del eje Y
        legend.title = element_text(size = 11), # Tamaño del título de la leyenda
        legend.text = element_text(size = 9.5),
        strip.text.x = element_text(size = 12),  # Tamaño del texto en el eje x
        strip.text.y = element_text(size = 12),
        text = element_text(family = "Arial")) +
  scale_fill_manual(values = phylum_colors) + 
  guides(fill = guide_legend(override.aes = list(size = 7.5)))
bar_plot
#ggsave("Plot_papper/bar_plot_relative.png", bar_plot, width = 5, height = 5)


#ALPHA DIVERSITy-----------------
#data frame con diversity
r_df_sfs <- data.frame(sample_data(r_ps_oc))
data_r_otu_sfs <- t(data.frame(otu_table(r_ps_oc)))
data_r_richness <- estimateR(data_r_otu_sfs)
S.evenness <- diversity(data_r_otu_sfs) / log(specnumber(data_r_otu_sfs))
S.shannon <- diversity(data_r_otu_sfs, index = "shannon")
sfs_alphadiv <- cbind(r_df_sfs, t(data_r_richness), S.shannon, S.evenness)
rm(r_df_sfs, data_r_otu_sfs, data_r_richness, S.evenness, S.shannon)

##grafico
D1 <- ggplot(sfs_alphadiv, aes(x = site, y = S.obs, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = NULL, y = "Richness", x =NULL, fill = "Depth",  tag = "A") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  scale_x_discrete(labels = c("Algarrobo" = "Algarrobo", "Las_Cruces" = "Las Cruces")) +
  theme_bw()
D2 <- ggplot(sfs_alphadiv, aes(x=site, y=S.chao1, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = NULL, y = "Chao1", x =NULL, fill = "Depth",  tag = "B") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  scale_x_discrete(labels = c("Algarrobo" = "Algarrobo", "Las_Cruces" = "Las Cruces")) +
  theme_bw()
D3 <- ggplot(sfs_alphadiv, aes(x=site, y=S.evenness, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = NULL, y = "Evenness", x ="Site", fill = "Depth",  tag = "C") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  scale_x_discrete(labels = c("Algarrobo" = "Algarrobo", "Las_Cruces" = "Las Cruces")) +
  theme_bw()
D4 <- ggplot(sfs_alphadiv, aes(x=site, y=S.shannon, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = NULL, y = "Shannon", x ="Site", fill = "Depth", tag = "D") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  scale_x_discrete(labels = c("Algarrobo" = "Algarrobo", "Las_Cruces" = "Las Cruces")) +
  theme_bw()

plot_div <-  D1 + D2 + D3 + D4 +
  plot_layout(guides = 'collect') & theme(
  plot.title = element_text(size = 24, hjust = 0.5, vjust = 0.5),  # Tamaño del título general
  axis.title = element_text(size = 22),  # Tamaño de los títulos de los ejes
  legend.title = element_text(size = 23), # Tamaño del título de la leyenda
  legend.text = element_text(size = 17),
  legend.key.size = unit(2, "lines"),
  axis.text.x = element_text(size = 17),   # Tamaño de las etiquetas del eje X
  axis.text.y = element_text(size = 17),
  axis.title.y = element_text(margin = margin(r = 5)) # Aumenta el margen derecho del título del eje Y
  
  ) & 
  guides(fill = guide_legend(override.aes = list(size = 5)))# Ajustar el tamaño de los símbolos
plot_div
#ggsave("Plot_papper/alpha_div.png", plot_div, width = 10, height = 8)


#GLM------------
###ver familias----------
hist(sfs_alphadiv$S.shannon)
hist(sfs_alphadiv$S.chao1)
hist(sfs_alphadiv$S.evenness)
hist(sfs_alphadiv$S.obs)

norm <- fitdist   (sfs_alphadiv$S.shannon, "norm")
lnormal <- fitdist(sfs_alphadiv$S.shannon, "lnorm")
exp <- fitdist    (sfs_alphadiv$S.shannon, "exp")
gamma <- fitdist  (sfs_alphadiv$S.shannon, "gamma")

# Cumulative distribution frequency, Gráficos de comparación CDF
(CDF.dens=cdfcomp(list(exp, norm, lnormal, gamma),
                  addlegend=T,main="",legendtext=c("Exp", "Normal", "lnorm", "gamma"),
                  plotstyle = "ggplot")+
    xlab("Densidad")+
    geom_line(size=0.8)+
    theme(axis.title=element_text(size=8), 
          axis.text = element_text(size=10), 
          legend.position = c(0.7,0.45),
          legend.text=element_text(size=6)))

#Gráficos de comparación QQ
(QQ.ZA.dens=qqcomp(list(exp, norm, lnormal, gamma),addlegend=F,main="",legendtext=c("Exp", "Normal", "lnorm", "gamma"),plotstyle = "ggplot")+
    theme_bw()+
    geom_jitter(size=2, height=0.2)+
    geom_line()+
    theme(axis.title=element_text(size=18), 
          axis.text = element_text(size=16), 
          title=element_blank(),
          legend.position = c(0.25,0.75),
          legend.text=element_text(size=14)))

###Homocedasticidad------------
car::leveneTest(S.shannon ~ site * ocean_depth,  data = sfs_alphadiv) # Usamos Levene porque es más robusto frente a distribuciones de variables aleatroias no normales. F max se basa en el supuesto de normalidad
car::leveneTest(S.evenness ~ site * ocean_depth, data = sfs_alphadiv) # Usamos Levene porque es más robusto frente a distribuciones de variables aleatroias no normales. F max se basa en el supuesto de normalidad
car::leveneTest(S.obs ~ site * ocean_depth,      data = sfs_alphadiv) # Usamos Levene porque es más robusto frente a distribuciones de variables aleatroias no normales. F max se basa en el supuesto de normalidad
car::leveneTest(S.chao1 ~ site * ocean_depth,    data = sfs_alphadiv) # Usamos Levene porque es más robusto frente a distribuciones de variables aleatroias no normales. F max se basa en el supuesto de normalidad


##GLM calculos-----
glm_shannon <- glm(S.shannon ~ site * ocean_depth,
                   family = gaussian(link = identity), data = sfs_alphadiv)

glm_evenness <- betareg(S.evenness ~ site + ocean_depth + site * ocean_depth,
                        data = sfs_alphadiv, link = "logit")

# Richness (Poisson or Negative Binomial regression)
# First check for overdispersion
poisson_model <- glm(S.obs ~ site * ocean_depth, 
                     family = poisson(link = "log"), data = sfs_alphadiv)
if (summary(poisson_model)$dispersion > 1) {
  # Overdispersion detected, use Negative Binomial
  glm_richness <- glm.nb(S.obs ~ site * ocean_depth, data = sfs_alphadiv)
} else {
  glm_richness <- poisson_model
}

# Chao1 (Gaussian or Gamma regression)
# Assess distribution of Chao1 values
if (skewness(sfs_alphadiv$S.chao1) > 1) {
  # Use Gamma if skewed
  glm_chao1 <- glm(S.chao1 ~ site * ocean_depth,
                   family = Gamma(link = "log"), data = sfs_alphadiv)
} else {
  # Use Gaussian if not skewed
  glm_chao1 <- glm(S.chao1 ~ site * ocean_depth,
                   family = gaussian(link = "identity"), data = sfs_alphadiv)
}

summary(glm_shannon)
summary(glm_evenness)
summary(glm_richness)
summary(glm_chao1)

``` 
