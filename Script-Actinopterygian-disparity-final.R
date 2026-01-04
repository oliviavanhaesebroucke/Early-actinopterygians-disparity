###############################################################################
############## Early actinopterygian morphological disparity ##################
################## Olivia Vanhaesebroucke - Octobre 2025 ######################
###############################################################################

getwd()
setwd("C:/Users/olivi/Documents/UQAR/DOC/Redaction-Articles-These/Article_1_Actinopterygian_disparity/Analyses-finales")

library(geomorph)
library(ggplot2)
library(Momocs)
library(ggfortify)
library(ggthemes)
library(ggdark)
library(Morpho)
library(readxl)
library(dplyr)
library(vegan)
library(divDyn)
library(dispRity)
library(scales)
library(ggpubr)
library(forcats)
library(readxl)
library(writexl)

########################################################################################################
############################################# COLOR CODES ##############################################
########################################################################################################

# Color for epochs #
cols_epoch <- c("Middle Devonian" = "#F1C868", "Upper Devonian" = "#F1E19D", "Mississippian" = "#678F66", "Pennsylvanian" = "#7EBCC6") #international chronostratigraphic chart colors

# Color for ages #
cols_age <- c("Eifelian" = "#F1D576", "Givetian" = "#F1E185", "Frasnian" = "#F2EDAD", "Famennian" = "#F2EDB3" , "Tournaisian" = "#8CB06C", "Viséan"= "#A6B96C", "Serpukhovian" = "#BFC26B", "Bashkirian" = "#99C2B5", "Moscovian" = "#B3CBB9", "Kasimovian" = "#BFD0C5", "Gzhelian" = "#CCD4C7") #color codes following the international chronostratigraphic chart 

########################################################################################################
################################# CONVEX HULL FUNCTION -- LAURENT HOULE ################################
########################################################################################################

# Creating convex hull: if there is only one observation for one category, the point is left alone in the graph
convex.hulls <- function(data, name.x, name.y, name.fill){
  
  library(tidyverse)
  library(grDevices)
  library(swaRm)
  
  data$X <- pca.class[[name.x]]
  data$Y <- data[[name.y]]
  data$FILL <- factor(data[[name.fill]])
  
  hull1 <- data %>%
    dplyr::slice(chull(X, Y))
  hull1 <- hull1[NULL,]
  SA <- c()
  ind <- 0
  for(j in 1:length(levels(data$FILL))){
    new <- filter(data, FILL == levels(data$FILL)[j])
    
    if(length(new$X) > 2){
      ind <- ind + 1
      hull.new <- new %>%
        dplyr::slice(chull(X, Y))
      
      Surf_area <- chull_area(hull.new$X,hull.new$Y)
      SA <- c(SA, Surf_area)
      names(SA)[ind] <- levels(data$FILL)[j]
      hull1 <- rbind(hull1,hull.new)
    }
    
  }
  
  return(list(table = hull1, surface_area = SA))
}
# end of the function

##############################################################################################################
##################################### SPECIES RICHNESS THROUGH TIME ##########################################
##############################################################################################################

data(stages) # loading dataset from divDyn package containing information about chronostratigraphic chart
stages_upd <- read_excel(file.choose(),1) #open stages-update.xlsx : the time for each stage from Ordovician, Silurian, Devonian and Carboniferous periods has been updated following the chronostratigraphic time chart 2024.

list_all <- read_excel(file.choose(),1) #open Matrix-all-species.xlsx
list_all <- list_all %>%
  mutate_if(is.character, as.factor)

list_actino <- read_excel(file.choose(),1) #open Matrix-coordinates.xlsx 
list_actino<- list_actino %>%
  mutate_if(is.character, as.factor)
list_actino <- list_actino[-c(7,35,40,59,63,91,100),]

list_all %>% mutate_at(c("max_ma", "min_ma"), as.numeric) #transform LAD/FAD column from character to numeric
list_all$me_ma <- apply(list_all[, c("max_ma", "min_ma")], 1, mean) # calculate the median age of each species

list_actino %>% mutate_at(c("max_ma", "min_ma"), as.numeric) #transform LAD/FAD column from character to numeric
list_actino$me_ma <- apply(list_actino[, c("max_ma", "min_ma")], 1, mean) # calculate the median age of each species

flDual_all <- fadlad(list_all, tax = "Species", age = c("max_ma", "min_ma")) # create a species x FAD/LAD matrix
flDual <- fadlad(list_actino, tax = "Species", age = c("max_ma", "min_ma")) # create a species x FAD/LAD matrix

list_all$mid <- stages$mid[list_all$Time_bin] #create a new column in the dataframe with the median age of the timebin
tsplot(stages_upd, shading="stage", boxes=c("short","series"),xlim=c(410,298.9), labels.args=list(cex=0.7), boxes.col=c("seriesCol", "systemCol")) #create background with the timechart
divDyn::ranges(list_all, tax="Species", bin=c("max_ma","min_ma"), labs=T, labels.args=list(cex=0.6, font = 3), occs=TRUE, filt="orig") #add duration of each species; font = 3 is italic

list_actino$mid <- stages$mid[list_actino$Time_bin] #create a new column in the dataframe with the median age of the timebin
tsplot(stages_upd, shading="stage", boxes=c("short","series"),xlim=c(410,298.9), labels.args=list(cex=0.7), boxes.col=c("seriesCol", "systemCol")) #create background with the timechart
divDyn::ranges(list_actino, tax="Species", bin=c("max_ma","min_ma"), labs=T, labels.args=list(cex=0.5, font = 3), occs=TRUE, filt="orig") #add duration of each species; font = 3 is italic

# Diversity through time #

actino_div <- divDyn(list_actino, bin = "Time_bin", tax = "Species")
all_div <- divDyn(list_all, bin = "Time_bin", tax = "Species")

tsplot(stages_upd, shading="stage", boxes=c("short", "series"), xlim=19:29, ylab="Richness (diversity)", ylim=c(0,80), labels.args = list(cex=0.7), boxes.col = c("seriesCol", "systemCol"))
lines(stages$mid[27:42], actino_div$divRT[27:42], col="blue", lwd = 3)
lines(stages$mid[27:42], all_div$divRT[27:42], col="black", lwd = 2)

##############################################################################################################
########################################### FULLBODY DISPARITY ############################################
##############################################################################################################

######################################### Import datasets ################################################

Actino_FB <- readland.tps("Fullbody_final.TPS",specID="imageID",negNA = TRUE,readcurves = TRUE,warnmsg = TRUE) # data for postcranial disparity

dim(Actino_FB)

name_FB <- dimnames(Actino_FB)[[3]]
name_FB

Period_Actino_FB <- read_excel(file.choose(), 1) #open Species-list-FB.xlsx
Period_Actino_FB

Period_Actino_FB<- Period_Actino_FB %>%
  mutate_if(is.character, as.factor)

# Convert curves into semi-landmarks #

sliders_FB <- read.table("Fullbody_final_sliders.nts")
sliders_FB <- as.matrix(sliders_FB)

########################################## GPA and PCA ####################################################

# Procrustes superimposition #
data.super_FB <- gpagen(Actino_FB, curves = sliders_FB, ProcD = FALSE)
attributes(data.super_FB)

plot(data.super_FB) 

# Principal component analyses #
pca_FB <- gm.prcomp(data.super_FB$coords)
pca_FB

plot_FB <- plot(pca_FB, axis1 = 1, axis2 = 2)
text(pca_FB[["x"]][,1], pca_FB[["x"]][,2], labels = name_FB) # PCA 1 vs 2

plot_PC2_PC3 <- plot(pca_FB, axis1=2, axis2=3)
text(pca_FB[["x"]][,2], pca_FB[["x"]][,3], labels = name_FB) # PCA 2 vs 3

# Saving PC scores #

PC.scores_FB <- pca_FB$x 
as.data.frame(PC.scores_FB) # Save PC scores as a data frame object

write.csv(PC.scores_FB,"PC.scores_FB.csv",row.names=TRUE) # Save PC scores as a csv file

# Scree plots #

var_expl_FB <- (pca_FB$sdev^2) / sum(pca_FB$sdev^2) * 100  

# Plot
plot(var_expl_FB, type = "b", pch = 19, col = "black",
     xlab = "Principal component",
     ylab = "Explained variance (%)",
     main = "Scree plot PC (Procrustes PCA)")

# threshold 5%
abline(h = 5, col = "red", lty = 2) # we will look at the first 5 axes

##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color epoch
ggplot(as.data.frame(PC.scores_FB), aes(x=PC.scores_FB[,1], y=PC.scores_FB[,2], label = name_FB, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Actino_FB$Epoque), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.3,0.5), y = c(-0.2,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 40.08 %", y = "PC2 = 16.27 %" ) 

#PC2 vs PC3 -- color epoch
ggplot(as.data.frame(PC.scores_FB), aes(x=PC.scores_FB[,2], y=PC.scores_FB[,3], label = name_FB, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Actino_FB$Epoque), color = "black", size = 6, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.2,0.3), y = c(-0.2,0.15)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 16.27 %", y = "PC3 = 8.33 %" )


#PC1 vs PC2 -- color age
ggplot(as.data.frame(PC.scores_FB), aes(x=PC.scores_FB[,1], y=PC.scores_FB[,2], label = name_FB, fontface = "italic")) +
  labs(fill="Age") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Actino_FB$Age), color = "black", size = 6, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_age)+
  lims(x=c(-0.3,0.5), y = c(-0.2,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 40.08 %", y = "PC2 = 16.27 %" ) 

#PC2 vs PC3 -- color age
ggplot(as.data.frame(PC.scores_FB), aes(x=PC.scores_FB[,2], y=PC.scores_FB[,3], label = name_FB, fontface = "italic")) +
  labs(fill="Age") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Actino_FB$Age), color = "black", size = 6, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_age)+
  lims(x=c(-0.2,0.3), y = c(-0.2,0.15)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 16.27 %", y = "PC3 = 8.33 %" )

################################################################################
########################## Warp Grids for each PC ##############################
################################################################################

# Use this to find which specimen is closest to the mean shape:
actino <- findMeanSpec(data.super_FB$coords) # turned out to be "Mesopoma pulchellum (51)"

# Estimate the mean shape for a set of aligned specimens
msh_actino <- mshape(data.super_FB$coords)

#Compare the minimum and maximum values to the global consensus:

# Full body 
plotRefToTarget(pca_FB$shapes$shapes.comp1$min, msh_actino) 
plotRefToTarget(msh_actino, pca_FB$shapes$shapes.comp1$max) # PC 1
plotRefToTarget(msh_actino, data.super_FB$coords[,,18], mag = 0.8) # target one specimen

plotRefToTarget(pca_FB$shapes$shapes.comp2$min, msh_actino) 
plotRefToTarget(msh_actino, pca_FB$shapes$shapes.comp2$max) # PC 2

plotRefToTarget(pca_FB$shapes$shapes.comp3$min, msh_actino) 
plotRefToTarget(msh_actino, pca_FB$shapes$shapes.comp3$max) # PC 3

par(mfrow = c(7, 12),
    mar = c(1, 1, 2, 1))

for (i in 1:84) {
  plotRefToTarget(
    msh_actino,
    data.super_FB$coords[,,i],
    mag = 0.4
  )
} 


#########################################################################################################
################################## MORPHOLOGICAL DISPARITY ##############################################
#########################################################################################################

################### Calculating convex hulls for the PC1 vs PC2 #######################

# Add classifier variables (group, age, habitats) by hand to .csv file
# loading .csv file with PC scores and classifiers:
pca.class <- read.csv(file.choose()) #Open PC.scores_FB.csv

# for the epochs
CHull_epoch <- convex.hulls(data = pca.class, name.x = "Comp4", name.y = "Comp5", name.fill = "Epoch")

CHull_epoch.table <- CHull_epoch$table
CHull_epoch$surface_area

# Weightened area
var_explained <- (pca_FB$sdev)^2 / sum(pca_FB$sdev^2)
var_explained # Explained variance
var_explained <- as.data.frame(var_explained)

# Explained variance for two axes - Simple additive weighting
saw_12 <- (sum(var_explained$var_explained[1:2]))/2
saw_13 <- (sum(var_explained$var_explained[1],var_explained$var_explained[3]))/2
saw_14 <- (sum(var_explained$var_explained[1],var_explained$var_explained[4]))/2
saw_15 <- (sum(var_explained$var_explained[1],var_explained$var_explained[5]))/2
saw_23 <- (sum(var_explained$var_explained[2],var_explained$var_explained[3]))/2
saw_24 <- (sum(var_explained$var_explained[2],var_explained$var_explained[4]))/2
saw_25 <- (sum(var_explained$var_explained[2],var_explained$var_explained[5]))/2
saw_34 <- (sum(var_explained$var_explained[3],var_explained$var_explained[4]))/2
saw_35 <- (sum(var_explained$var_explained[3],var_explained$var_explained[5]))/2
saw_45 <- (sum(var_explained$var_explained[4],var_explained$var_explained[5]))/2

# Explained variance for two axes - Weighted sum of pairwise area
wspa_12 <- var_explained$var_explained[1]*var_explained$var_explained[2]
wspa_13 <- var_explained$var_explained[1]*var_explained$var_explained[3]
wspa_14 <- var_explained$var_explained[1]*var_explained$var_explained[4]
wspa_15 <- var_explained$var_explained[1]*var_explained$var_explained[5]
wspa_23 <- var_explained$var_explained[2]*var_explained$var_explained[3]
wspa_24 <- var_explained$var_explained[2]*var_explained$var_explained[4]
wspa_25 <- var_explained$var_explained[2]*var_explained$var_explained[5]
wspa_34 <- var_explained$var_explained[3]*var_explained$var_explained[4]
wspa_35 <- var_explained$var_explained[3]*var_explained$var_explained[5]
wspa_45 <- var_explained$var_explained[4]*var_explained$var_explained[5]

# Create data frame with the PCs pairs weight
pc_pairs <- c("PC1 vs PC2", "PC1 vs PC3", "PC1 vs PC4", "PC1 vs PC5", "PC2 vs PC3", "PC2 vs PC4", "PC2 vs PC5", "PC3 vs PC4", "PC3 vs PC5", "PC4 vs PC5")

SAW <- c(saw_12, saw_13, saw_14, saw_15, saw_23, saw_24, saw_25, saw_34, saw_35, saw_45)
WSPA <- c(wspa_12, wspa_13, wspa_14, wspa_15, wspa_23, wspa_24, wspa_25, wspa_34, wspa_35, wspa_45)

weight <- data.frame(pc_pairs, SAW, WSPA)

print(weight)

################################# Morphological disparity analyses ###########################################

gdf_FB_epoch <- geomorph.data.frame(data.super_FB, species = Period_Actino_FB$Species, Time = Period_Actino_FB$Epoque) #gdf for the epochs

gdf_FB_age <- geomorph.data.frame(data.super_FB, species = Period_Actino_FB$Species, Time = Period_Actino_FB$Age) #gdf for the ages

# Sum of Procrustes variances : Epochs #########################

SOV_mean_FB_epoch <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_FB_epoch, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_FB_epoch)

SOV_FB_epoch <- morphol.disparity(coords~Time,groups=~Time, data = gdf_FB_epoch, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_FB_epoch)

# figure of the variation of the disparity through time

SoV_FB.epoch <- c("0.02016574", "0.02014144", "0.04327367", "0.03808941") # copy the SOV_FB_epoch results
Epoch <- c("Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian")
SoV_FB.epoch <- as.numeric(as.character(SoV_FB.epoch))
Epoch <- as.factor(as.character(Epoch))

DF_SoV.FB.epoch <- data.frame(Epoch = Epoch, SoV = SoV_FB.epoch)
DF_SoV.FB.epoch # Create data frame Epoch x SoV 

level_order <- factor(DF_SoV.FB.epoch$Epoch, level = c("Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian"))

ggplot(DF_SoV.FB.epoch, aes(x=level_order, y=SoV, group = 1)) +
  geom_line(color = "black", linewidth=1) +
  geom_point(color = cols_epoch, size=5)+
  scale_y_continuous(labels = label_number(accuracy = 0.001),limits = c(0,0.05))+
  xlab("Epoch") +
  ylab("Procrustes variance")  +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust=0.5, size = 14))


# Contribution of each time bin : Pie chart

PPV_PC.epoch <- c("0.03106806","0.05330455","0.64354710","0.27208029") # copy the SOV_mean_FB_epoch results
Epoch <- c("Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian")
PPV_PC.epoch <- as.numeric(as.character(PPV_PC.epoch))
Epoch <- as.factor(as.character(Epoch))

DF_PPV.epoch <- data.frame(Epoch = Epoch, PPV = PPV_PC.epoch)
DF_PPV.epoch # Create data frame Epoch x Proportion of variance

level_order <- factor(DF_PPV.epoch$Epoch, level = c("Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian"))

# Compute percentages
DF_PPV.epoch$fraction = DF_PPV.epoch$PPV / sum(DF_PPV.epoch$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.epoch$ymax = cumsum(DF_PPV.epoch$fraction)

# Compute the bottom of each rectangle
DF_PPV.epoch$ymin = c(0, head(DF_PPV.epoch$ymax, n=-1))

DF_PPV.epoch$labelPosition <- (DF_PPV.epoch$ymax + DF_PPV.epoch$ymin) / 2

# Compute a good label
DF_PPV.epoch$label <- paste0(DF_PPV.epoch$PPV)

# Make the plot
ggplot(DF_PPV.epoch, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Epoch)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_epoch) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.epoch, aes(x="", y=PPV, fill=Epoch)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_epoch) +
  theme_void() #Pie chart

# Sum of Procrustes variances : age ######################

SOV_mean_FB_age <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_FB_age, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_FB_age)

SOV_FB_age <- morphol.disparity(coords~Time,groups=~Time, data = gdf_FB_age, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_FB_age)

# Graphical representation #

# figure of the variation of the disparity through time
SoV_FB.age <- c("0.000000000", "0.008185613", "0.020368877", "0.007703818", "0.037190154", "0.039271170", "0.045587387", "0.011973400", "0.041385655", "0.036439432", "0.030227978") # copy the SOV_FB_age results
Age <- c("Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Viséan", "Serpukhovian", "Bashkirian", "Moscovian", "Kasimovian", "Gzhelian")
SoV_FB.age <- as.numeric(as.character(SoV_FB.age))
Age <- as.factor(as.character(Age))

DF_SoV.FB.age <- data.frame(Age = Age, SoV = SoV_FB.age)
DF_SoV.FB.age # Create data frame Age x SoV

level_order <- factor(DF_SoV.FB.age$Age, level = c("Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Viséan", "Serpukhovian", "Bashkirian", "Moscovian", "Kasimovian", "Gzhelian"))

ggplot(DF_SoV.FB.age, aes(x=level_order, y=SoV, group = 1)) +
  geom_line(color = "black", linewidth=1) +
  geom_point(color = cols_age, size=5)+
  scale_y_continuous(labels = label_number(accuracy = 0.001),limits = c(0,0.05))+
  xlab("Age") +
  ylab("Procrustes variance")  +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust=0.5, size = 14))

# Contribution of each time bin : Pie chart

PPV_FB.age <- c("0.01790052", "0.01316754", "0.03717931", "0.01612524", "0.08347281", "0.29468840", "0.26538589", "0.02239060", "0.15894886", "0.03925409", "0.05148675") # copy the SOV_mean_FB_age results
Age <- c("Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Viséan", "Serpukhovian", "Bashkirian", "Moscovian", "Kasimovian", "Gzhelian")
PPV_FB.age <- as.numeric(as.character(PPV_FB.age))
Age <- as.factor(as.character(Age))

DF_PPV.age <- data.frame(Age = Age, PPV = PPV_FB.age)
DF_PPV.age # Create data frame Age x Proportion of variance

level_order <- factor(DF_PPV.age$Age, level = c("Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Viséan", "Serpukhovian", "Bashkirian", "Moscovian", "Kasimovian", "Gzhelian"))

# Compute percentages
DF_PPV.age$fraction = DF_PPV.age$PPV / sum(DF_PPV.age$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.age$ymax = cumsum(DF_PPV.age$fraction)

# Compute the bottom of each rectangle
DF_PPV.age$ymin = c(0, head(DF_PPV.age$ymax, n=-1))

DF_PPV.age$labelPosition <- (DF_PPV.age$ymax + DF_PPV.age$ymin) / 2

# Compute a good label
DF_PPV.age$label <- paste0(DF_PPV.age$PPV)

# Make the plot
ggplot(DF_PPV.age, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Age)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_age) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.age, aes(x="", y=PPV, fill=Age)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_age) +
  theme_void() #Pie chart

##############################################################################################################
############################################## SKULL DISPARITY ###############################################
##############################################################################################################

######################################### Import datasets ################################################

Actino_SL <- readland.tps("Skull-lateral-R.TPS",specID="imageID",negNA = TRUE,readcurves = TRUE,warnmsg = TRUE) # data for skull disparity

dim(Actino_SL)

name_SL <- dimnames(Actino_SL)[[3]]
name_SL

Period_Actino_SL <- read_excel(file.choose(), 1) #open Species-list-SL.xlsx
Period_Actino_SL

Period_Actino_SL<- Period_Actino_SL %>%
  mutate_if(is.character, as.factor)

# Convert curves into semi-landmarks #

sliders_SL <- read.table("Skull-lateral-sliders.nts")
sliders_SL <- as.matrix(sliders_SL)

########################################## GPA and PCA ####################################################

# Procrustes superimposition #
data.super_SL <- gpagen(Actino_SL, curves = sliders_SL, ProcD = FALSE)
attributes(data.super_SL)

plot(data.super_SL) 

# Principal component analyses #
pca_SL <- gm.prcomp(data.super_SL$coords)
pca_SL

plot_SL <- plot(pca_SL, axis1 = 1, axis2 = 2)
text(pca_SL[["x"]][,1], pca_SL[["x"]][,2], labels = name_SL) # PCA 1 vs 2

plot_PC2_PC3 <- plot(pca_SL, axis1=2, axis2=3)
text(pca_SL[["x"]][,2], pca_SL[["x"]][,3], labels = name_SL) # PCA 2 vs 3

# Saving PC scores #

PC.scores_SL <- pca_SL$x 
as.data.frame(PC.scores_SL) # Save PC scores as a data frame object

write.csv(PC.scores_SL,"PC.scores_SL.csv",row.names=TRUE) # Save PC scores as a csv file

# Scree plots #

var_expl_SL <- (pca_SL$sdev^2) / sum(pca_SL$sdev^2) * 100  

# Plot
plot(var_expl_SL, type = "b", pch = 19, col = "black",
     xlab = "Principal component",
     ylab = "Explained variance (%)",
     main = "Scree plot PC (Procrustes PCA)")

# threshold 5%
abline(h = 5, col = "red", lty = 2) # we will look at the first 4 axes

##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color epoch
ggplot(as.data.frame(PC.scores_SL), aes(x=PC.scores_SL[,1], y=PC.scores_SL[,2], label = name_SL, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Actino_SL$Epoque), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.4,0.5), y = c(-0.2,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 48.15 %", y = "PC2 = 13.90 %" ) 

#PC2 vs PC3 -- color epoch
ggplot(as.data.frame(PC.scores_SL), aes(x=PC.scores_SL[,2], y=PC.scores_SL[,3], label = name_SL, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Actino_SL$Epoque), color = "black", size = 6, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.4,0.4), y = c(-0.15,0.35)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 13.90 %", y = "PC3 = 10.23 %" )


#PC1 vs PC2 -- color age
ggplot(as.data.frame(PC.scores_SL), aes(x=PC.scores_SL[,1], y=PC.scores_SL[,2], label = name_SL, fontface = "italic")) +
  labs(fill="Age") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Actino_SL$Age), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_age)+
  lims(x=c(-0.4,0.5), y = c(-0.2,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 48.15 %", y = "PC2 = 13.90 %" ) 

#PC2 vs PC3 -- color age
ggplot(as.data.frame(PC.scores_SL), aes(x=PC.scores_SL[,2], y=PC.scores_SL[,3], label = name_SL, fontface = "italic")) +
  labs(fill="Age") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Actino_SL$Age), color = "black", size = 6, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_age)+
  lims(x=c(-0.4,0.4), y = c(-0.15,0.35)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 13.90 %", y = "PC3 = 10.23 %" )

################################################################################
########################## Warp Grids for each PC ##############################
################################################################################

# Use this to find which specimen is closest to the mean shape:
skull <- findMeanSpec(data.super_SL$coords) # turned out to be "Yogoniscus gulo (86)"

# Estimate the mean shape for a set of aligned specimens
msh_skull <- mshape(data.super_SL$coords)

#Compare the minimum and maximum values to the global consensus:

# Skull
plotRefToTarget(pca_SL$shapes$shapes.comp1$min, msh_skull) 
plotRefToTarget(msh_skull, pca_SL$shapes$shapes.comp1$max) # PC 1
plotRefToTarget(msh_actino, data.super_SL$coords[,,18], mag = 0.8) # target one specimen

plotRefToTarget(pca_SL$shapes$shapes.comp2$min, msh_skull) 
plotRefToTarget(msh_skull, pca_SL$shapes$shapes.comp2$max) # PC 2

plotRefToTarget(pca_SL$shapes$shapes.comp3$min, msh_skull) 
plotRefToTarget(msh_skull, pca_SL$shapes$shapes.comp3$max) # PC 3

par(mfrow = c(8, 11),
    mar = c(1, 1, 2, 1))

for (i in 1:86) {
  plotRefToTarget(
    msh_skull,
    data.super_SL$coords[,,i],
    mag = 0.8
  )
} # All TPS grids

#########################################################################################################
################################## MORPHOLOGICAL DISPARITY ##############################################
#########################################################################################################

################### Calculating convex hulls for the PC1 vs PC2 #######################

# Add classifier variables (group, age, habitats) by hand to .csv file
# loading .csv file with PC scores and classifiers:
pca.class <- read.csv(file.choose()) #Open PC.scores_SL.csv

# for the epochs
CHull_epoch <- convex.hulls(data = pca.class, name.x = "Comp3", name.y = "Comp4", name.fill = "Age")

CHull_epoch.table <- CHull_epoch$table
CHull_epoch$surface_area

# Weightened area
var_explained <- (pca_SL$sdev)^2 / sum(pca_SL$sdev^2)
var_explained # Explained variance
var_explained <- as.data.frame(var_explained)

# Explained variance for two axes - Simple additive weighting
saw_12 <- (sum(var_explained$var_explained[1:2]))/2
saw_13 <- (sum(var_explained$var_explained[1],var_explained$var_explained[3]))/2
saw_14 <- (sum(var_explained$var_explained[1],var_explained$var_explained[4]))/2
saw_23 <- (sum(var_explained$var_explained[2],var_explained$var_explained[3]))/2
saw_24 <- (sum(var_explained$var_explained[2],var_explained$var_explained[4]))/2
saw_34 <- (sum(var_explained$var_explained[3],var_explained$var_explained[4]))/2

# Explained variance for two axes - Weighted sum of pairwise area
wspa_12 <- var_explained$var_explained[1]*var_explained$var_explained[2]
wspa_13 <- var_explained$var_explained[1]*var_explained$var_explained[3]
wspa_14 <- var_explained$var_explained[1]*var_explained$var_explained[4]
wspa_23 <- var_explained$var_explained[2]*var_explained$var_explained[3]
wspa_24 <- var_explained$var_explained[2]*var_explained$var_explained[4]
wspa_34 <- var_explained$var_explained[3]*var_explained$var_explained[4]

# Create data frame with the PCs pairs weight
pc_pairs <- c("PC1 vs PC2", "PC1 vs PC3", "PC1 vs PC4", "PC2 vs PC3", "PC2 vs PC4", "PC3 vs PC4")

SAW <- c(saw_12, saw_13, saw_14, saw_23, saw_24, saw_34)
WSPA <- c(wspa_12, wspa_13, wspa_14, wspa_23, wspa_24, wspa_34)

weight <- data.frame(pc_pairs, SAW, WSPA)

print(weight)

################################# Morphological disparity analyses ###########################################

gdf_SL_epoch <- geomorph.data.frame(data.super_SL, species = Period_Actino_SL$Species, Time = Period_Actino_SL$Epoque) #gdf for the epochs

gdf_SL_age <- geomorph.data.frame(data.super_SL, species = Period_Actino_SL$Species, Time = Period_Actino_SL$Age) #gdf for the ages

# Sum of Procrustes variances : Epochs #########################

SOV_mean_SL_epoch <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_SL_epoch, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_SL_epoch)

SOV_SL_epoch <- morphol.disparity(coords~Time,groups=~Time, data = gdf_SL_epoch, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_SL_epoch)

# figure of the variation of the disparity through time

SoV_SL.epoch <- c("0.009729606", "0.020682142", "0.038579207", "0.042323857") # copy the SOV_SL_epoch results
Epoch <- c("Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian")
SoV_SL.epoch <- as.numeric(as.character(SoV_SL.epoch))
Epoch <- as.factor(as.character(Epoch))

DF_SoV.SL.epoch <- data.frame(Epoch = Epoch, SoV = SoV_SL.epoch)
DF_SoV.SL.epoch # Create data frame Epoch x SoV 

level_order <- factor(DF_SoV.SL.epoch$Epoch, level = c("Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian"))

ggplot(DF_SoV.SL.epoch, aes(x=level_order, y=SoV, group = 1)) +
  geom_line(color = "black", linewidth=1) +
  geom_point(color = cols_epoch, size=5)+
  scale_y_continuous(labels = label_number(accuracy = 0.001),limits = c(0,0.05))+
  xlab("Epoch") +
  ylab("Procrustes variance")  +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust=0.5, size = 14))


# Contribution of each time bin : Pie chart

PPV_SL.epoch <- c("0.05537919","0.12536785","0.56836450","0.25088846") # copy the SOV_mean_SL_epoch results
Epoch <- c("Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian")
PPV_SL.epoch <- as.numeric(as.character(PPV_SL.epoch))
Epoch <- as.factor(as.character(Epoch))

DF_PPV.epoch <- data.frame(Epoch = Epoch, PPV = PPV_SL.epoch)
DF_PPV.epoch # Create data frame Epoch x Proportion of variance

level_order <- factor(DF_PPV.epoch$Epoch, level = c("Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian"))

# Compute percentages
DF_PPV.epoch$fraction = DF_PPV.epoch$PPV / sum(DF_PPV.epoch$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.epoch$ymax = cumsum(DF_PPV.epoch$fraction)

# Compute the bottom of each rectangle
DF_PPV.epoch$ymin = c(0, head(DF_PPV.epoch$ymax, n=-1))

DF_PPV.epoch$labelPosition <- (DF_PPV.epoch$ymax + DF_PPV.epoch$ymin) / 2

# Compute a good label
DF_PPV.epoch$label <- paste0(DF_PPV.epoch$PPV)

# Make the plot
ggplot(DF_PPV.epoch, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Epoch)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_epoch) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.epoch, aes(x="", y=PPV, fill=Epoch)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_epoch) +
  theme_void() #Pie chart

# Sum of Procrustes variances : age ######################

SOV_mean_SL_age <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_SL_age, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_SL_age)

SOV_SL_age <- morphol.disparity(coords~Time,groups=~Time, data = gdf_SL_age, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_SL_age)

# Graphical representation #

# figure of the variation of the disparity through time
SoV_SL.age <- c("0.000000000", "0.008667792", "0.010614208", "0.029955005", "0.017780439", "0.044414923", "0.027499534", "0.016938485", "0.033235288", "0.060665158", "0.014415969") # copy the SOV_SL_age results
Age <- c("Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Viséan", "Serpukhovian", "Bashkirian", "Moscovian", "Kasimovian", "Gzhelian")
SoV_SL.age <- as.numeric(as.character(SoV_SL.age))
Age <- as.factor(as.character(Age))

DF_SoV.SL.age <- data.frame(Age = Age, SoV = SoV_SL.age)
DF_SoV.SL.age # Create data frame Age x SoV

level_order <- factor(DF_SoV.SL.age$Age, level = c("Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Viséan", "Serpukhovian", "Bashkirian", "Moscovian", "Kasimovian", "Gzhelian"))

ggplot(DF_SoV.SL.age, aes(x=level_order, y=SoV, group = 1)) +
  geom_line(color = "black", linewidth=1) +
  geom_point(color = cols_age, size=5)+
  scale_y_continuous(labels = label_number(accuracy = 0.001),limits = c(0,0.065))+
  xlab("Age") +
  ylab("Procrustes variance")  +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust=0.5, size = 14))

# Contribution of each time bin : Pie chart

PPV_SL.age <- c("0.01556666", "0.03981253", "0.06604196", "0.05932589", "0.06873071", "0.36385291", "0.13578088", "0.02532063", "0.13789584", "0.07233764", "0.01533435") # copy the SOV_mean_SL_age results
Age <- c("Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Viséan", "Serpukhovian", "Bashkirian", "Moscovian", "Kasimovian", "Gzhelian")
PPV_SL.age <- as.numeric(as.character(PPV_SL.age))
Age <- as.factor(as.character(Age))

DF_PPV.age <- data.frame(Age = Age, PPV = PPV_SL.age)
DF_PPV.age # Create data frame Age x Proportion of variance

level_order <- factor(DF_PPV.age$Age, level = c("Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Viséan", "Serpukhovian", "Bashkirian", "Moscovian", "Kasimovian", "Gzhelian"))

# Compute percentages
DF_PPV.age$fraction = DF_PPV.age$PPV / sum(DF_PPV.age$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.age$ymax = cumsum(DF_PPV.age$fraction)

# Compute the bottom of each rectangle
DF_PPV.age$ymin = c(0, head(DF_PPV.age$ymax, n=-1))

DF_PPV.age$labelPosition <- (DF_PPV.age$ymax + DF_PPV.age$ymin) / 2

# Compute a good label
DF_PPV.age$label <- paste0(DF_PPV.age$PPV)

# Make the plot
ggplot(DF_PPV.age, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Age)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_age) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.age, aes(x="", y=PPV, fill=Age)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_age) +
  theme_void() #Pie chart

###############################################################################
############################ SOV through time #################################
###############################################################################

SOV <- data.frame(Time_bin = c(32:42), 
                  SOV_FB = c(0.000000000, 0.008185613, 0.020368877, 0.007703818, 0.037190154, 0.039271170, 0.045587387, 0.011973400, 0.041385655, 0.036439432, 0.030227978), 
                  SOV_SL = c(0.000000000, 0.008667792, 0.010614208, 0.029955005, 0.017780439, 0.044414923, 0.027499534, 0.016938485, 0.033235288, 0.060665158, 0.014415969))


tsplot(stages_upd, shading="stage", boxes=c("short", "series"), xlim=19:29, ylab="SOV", ylim=c(0,0.07), labels.args = list(cex=1), boxes.col = c("seriesCol", "systemCol"))

lines(stages_upd$mid[19:29], SOV$SOV_FB, col="blue", lwd = 4)
points(stages_upd$mid[19:29], SOV$SOV_FB, col="blue", pch = 15, cex = 2)
lines(stages_upd$mid[19:29], SOV$SOV_SL, col="#03be36", lwd = 4)
points(stages_upd$mid[19:29], SOV$SOV_SL, col="#03be36", pch = 17, cex = 2)


################################################################################
############################ Disparity for equal time bins #####################
################################################################################

# Full body disparity

Occurence_FB <- read_excel(file.choose(),1) #loading FAD and LAD - Open FB-Age.xlsx
Occurence_FB<-Occurence_FB %>%
  mutate_if(is.character, as.factor)

vec <- Occurence_FB$Species
Occurence_FB <- Occurence_FB[,-1]
rownames(Occurence_FB) <- vec

######### Find mismatch
which(!rownames(PC.scores_FB) %in% rownames(Occurence_FB))

## Issue - one mismatch in tip names. 
FADLAD.names <- rownames(Occurence_FB)
FADLAD.names
FADLAD.names[X] #X the number of the species mismatched

PC.names <- rownames(PC.scores_FB)
PC.names[X]

## First, convert FADLAD to data frame
Occurence_FB <- as.data.frame(Occurence_FB)

## Then correct name
FADLAD.names[X] <- PC.names[X]
rownames(Occurence_FB) <- FADLAD.names

####### Correction finished

chronosubset_FB <- chrono.subsets(data = PC.scores_FB, method = "discrete", time = 11, FADLAD = Occurence_FB) # subdivise species into 11 equal timebins
results_chrono <- dispRity(boot.matrix(chronosubset_FB, bootstraps = 100), metric =  function(X) return(sum(X^2)/nrow(X))) #sum of procrustes variance
results_chrono <- dispRity(boot.matrix(chronosubset_FB, bootstraps = 100), metric = c(sum,variances)) #sum of variance
 
summary(results_chrono)

get.disparity(results_chrono)
plot(results_chrono, type = "continuous")

# Skull disparity

Occurence_SL <- read_excel(file.choose(),1) #loading FAD and LAD - Open SL-Age.xlsx
Occurence_SL<-Occurence_SL %>%
  mutate_if(is.character, as.factor)

vec <- Occurence_SL$Species
Occurence_SL <- Occurence_SL[,-1]
rownames(Occurence_SL) <- vec

######### Find mismatch
which(!rownames(PC.scores_SL) %in% rownames(Occurence_SL))

## Issue - one mismatch in tip names. 
FADLAD.names <- rownames(Occurence_SL)
FADLAD.names
FADLAD.names[c(46,60)] #46 and 60 are mismatched

PC.names <- rownames(PC.scores_SL)
PC.names[c(46,60)]

## First, convert FADLAD to data frame
Occurence_SL <- as.data.frame(Occurence_SL)

## Then correct name
FADLAD.names[c(46,60)] <- PC.names[c(46,60)]
rownames(Occurence_SL) <- FADLAD.names

chronosubset_SL <- chrono.subsets(data = PC.scores_SL, method = "discrete", time = 8, FADLAD = Occurence_SL) # subdivise species into 11 equal timebins
results_chrono <- dispRity(boot.matrix(chronosubset_SL, bootstraps = 100), metric =  function(X) return(sum(X^2)/nrow(X))) #sum of procrustes variance
results_chrono <- dispRity(boot.matrix(chronosubset_SL, bootstraps = 100), metric = c(sum,variances)) #sum of variance

summary(results_chrono)

get.disparity(results_chrono)
plot(results_chrono, type = "continuous")

#########################################################################################################
######################################### PALAEOGEOGRAPHY ###############################################
#########################################################################################################

coordinates <- read_excel(file.choose(),1) # open file Matrix-coordinates.xlsx

coordinates <- coordinates %>%
  mutate_if(is.character, as.factor)
coordinates <- coordinates[,c(1,12,13)]

coordinates <- as.data.frame(coordinates)

rownames(coordinates) <- coordinates$Species
coordinates$Species <- NULL

# https://gplates.github.io/rgplates/

library(sf)
library(rgplates)
library(geojsonsf)

############### Transform actual coordinates into paleocoordinates ######################

# Eifelian (393 Mya)

paleocoords_eifelian <- rgplates::reconstruct(coordinates[18, c("lng", "lat")], 
                                              age = 393, model = "PALEOMAP")
sfObj_eifelian <- sf::st_as_sf(as.data.frame(paleocoords_eifelian), coords=c("paleolong", "paleolat"))
st_crs(sfObj_eifelian) <- "WGS84" # paleocoordinates are long-lat

# Givetian (388 Mya)

paleocoords_givetian <- rgplates::reconstruct(coordinates[c(17,18,29,72,73,99,100), c("lng", "lat")], 
                                              age = 388, model = "PALEOMAP")
sfObj_givetian <- sf::st_as_sf(as.data.frame(paleocoords_givetian), coords=c("paleolong", "paleolat"))
st_crs(sfObj_givetian) <- "WGS84" # paleocoordinates are long-lat

# Frasnian (383 Mya)

paleocoords_frasnian <- rgplates::reconstruct(coordinates[c(16,17,43,50,69,70,71,72,73,85), c("lng", "lat")], age = 383, model = "PALEOMAP")
sfObj_frasnian <- sf::st_as_sf(as.data.frame(paleocoords_frasnian), coords=c("paleolong", "paleolat"))
st_crs(sfObj_frasnian) <- "WGS84" # paleocoordinates are long-lat

# Famennian (372 Mya)

paleocoords_famennian <- rgplates::reconstruct(coordinates[c(25,34,35,54,76,106), c("lng", "lat")], 
                                               age = 372, model = "PALEOMAP")
sfObj_famennian <- sf::st_as_sf(as.data.frame(paleocoords_famennian), coords=c("paleolong", "paleolat"))
st_crs(sfObj_famennian) <- "WGS84" # paleocoordinates are long-lat

# Tournaisian (359 Mya)

paleocoords_tour <- rgplates::reconstruct(coordinates[c(6,7,30,34,35,38,42,46,47,53,57,75,84,97,101), c("lng", "lat")], age = 359, model = "PALEOMAP")
sfObj_tour <- sf::st_as_sf(as.data.frame(paleocoords_tour), coords=c("paleolong", "paleolat"))
st_crs(sfObj_tour) <- "WGS84" # paleocoordinates are long-lat

# Viséan (347 Mya) 

visean <- read_excel(file.choose(),1) # open matrix-coordinates-visean

visean <- visean %>%
  mutate_if(is.character, as.factor)
visean <- visean[,c(1,12,13)]

visean <- as.data.frame(visean)
rownames(visean) <- visean$Species
visean$Species <- NULL

paleocoords_visean1 <- rgplates::reconstruct(visean[c(1:15),c("lng","lat")], age = 347, model = "PALEOMAP")
sfObj_visean1 <- sf::st_as_sf(as.data.frame(paleocoords_visean1), coords=c("paleolong", "paleolat"))
st_crs(sfObj_visean1) <- "WGS84" # paleocoordinates are long-lat

paleocoords_visean2 <- rgplates::reconstruct(visean[c(16:34),c("lng","lat")], age = 347, model = "PALEOMAP")
sfObj_visean2 <- sf::st_as_sf(as.data.frame(paleocoords_visean2), coords=c("paleolong", "paleolat"))
st_crs(sfObj_visean2) <- "WGS84" # paleocoordinates are long-lat

# Serpukhovian (331 Mya)

paleocoords_serp <- rgplates::reconstruct(coordinates[c(4,9,24,27,28,39,40,41,45,51,52,55,56,58,59,62,63,81,87,88,108,109,112,113), c("lng", "lat")], age = 331, model = "PALEOMAP")
sfObj_serp <- sf::st_as_sf(as.data.frame(paleocoords_serp), coords=c("paleolong", "paleolat"))
st_crs(sfObj_serp) <- "WGS84" # paleocoordinates are long-lat

# Bashkirian (323 Mya)

paleocoords_bash <- rgplates::reconstruct(coordinates[c(20,22,107,110), c("lng", "lat")], age = 323, model = "PALEOMAP")
sfObj_bash <- sf::st_as_sf(as.data.frame(paleocoords_bash), coords=c("paleolong", "paleolat"))
st_crs(sfObj_bash) <- "WGS84" # paleocoordinates are long-lat

# Moscovian (315 Mya)

paleocoords_mosco <- rgplates::reconstruct(coordinates[c(19,20,31,33,48,61,67,68,77,82,90,91,92,96,107,110), c("lng", "lat")], age = 315, model = "PALEOMAP")
sfObj_mosco <- sf::st_as_sf(as.data.frame(paleocoords_mosco), coords=c("paleolong", "paleolat"))
st_crs(sfObj_mosco) <- "WGS84" # paleocoordinates are long-lat

# Kasimovian (307 Mya)

paleocoords_kasim <- rgplates::reconstruct(coordinates[c(3,12,78,92,96,104), c("lng", "lat")], age = 307, model = "PALEOMAP")
sfObj_kasim <- sf::st_as_sf(as.data.frame(paleocoords_kasim), coords=c("paleolong", "paleolat"))
st_crs(sfObj_kasim) <- "WGS84" # paleocoordinates are long-lat

# Gzhelian (304 Mya)

paleocoords_gzh <- rgplates::reconstruct(coordinates[c(1,21,32,89,98), c("lng", "lat")], age = 304, model = "PALEOMAP")
sfObj_gzh <- sf::st_as_sf(as.data.frame(paleocoords_gzh), coords=c("paleolong", "paleolat"))
st_crs(sfObj_gzh) <- "WGS84" # paleocoordinates are long-lat

################################## Construction of the map ########################################

coastlines <- rgplates::reconstruct("coastlines", age = 359, model = "PALEOMAP") #create coastlines at a given time

coastlines_sf <- sf::st_as_sf(coastlines) # convert into sf

epsg <- "ESRI:54030"  # projection Robinson
coastsRob <- sf::st_transform(coastlines_sf, crs = epsg)
sfObjRob  <- sf::st_transform(sfObj_gzh, crs = epsg)

# trace the map with the coordinates
plot(st_geometry(coastsRob), col = "black", border = "gray50", main = "304 Ma (PALEOMAP)")
plot(st_geometry(sfObjRob), col = "red", pch = 19, add = TRUE)


############################### Projection on paleomap from Scotese ####################################

library(jpeg)

map <- readJPEG("Gzhelian_300.jpg")

ggplot() +
  # backgroung : JPEG image 
  annotation_raster(map, xmin = -180, xmax = 180, ymin = -90, ymax = 90) +
  
  # paleocoordinates
  geom_sf(data = sfObj_gzh, color = "red", fill = "red", shape = 21, size = 3) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  theme_void() +
  labs(title = "300 Ma — Scotese's PALEOMAP")

################################################################################
############################## Geographic spreading ############################
################################################################################

# Spatial distribution of species 

df_eifelian <- as.data.frame(paleocoords_eifelian)
df_eifelian$Age <- "Eifelian"
df_givetian <- as.data.frame(paleocoords_givetian)
df_givetian$Age <- "Givetian"
df_frasnian <- as.data.frame(paleocoords_frasnian)
df_frasnian$Age <- "Frasnian"
df_famennian <- as.data.frame(paleocoords_famennian)
df_famennian$Age <- "Famennian"

df_tour <- as.data.frame(paleocoords_tour)
df_tour$Age <- "Tournaisian"
df_visean <- as.data.frame(paleocoords_visean1)
df_visean$Age <- "Viséan"
df_visean2 <- as.data.frame(paleocoords_visean2)
df_visean2$Age <- "Viséan"
df_serp <- as.data.frame(paleocoords_serp)
df_serp$Age <- "Serpukhovian"
df_bash <- as.data.frame(paleocoords_bash)
df_bash$Age <- "Bashkirian"
df_mosco <- as.data.frame(paleocoords_mosco)
df_mosco$Age <- "Moscovian"
df_kasim <- as.data.frame(paleocoords_kasim)
df_kasim$Age <- "Kasimovian"
df_gzh <- as.data.frame(paleocoords_gzh)
df_gzh$Age <- "Gzhelian"

data_all <- rbind(df_eifelian,df_givetian,df_frasnian,df_famennian,df_tour,df_visean,df_visean2,df_serp,df_bash,df_mosco,df_kasim,df_gzh)

data_all$paleolat <- as.numeric(data_all$paleolat)
data_all$paleolong <- as.numeric(data_all$paleolong)

ggplot(data_all, aes(x = paleolong, y = paleolat, fill = Age)) +
  geom_point(shape = 21, color = "black", size = 6, alpha = 1) +
  coord_fixed() +
  theme_minimal() +
  scale_fill_manual(
    name = "Âge géologique",
    values = cols_age
  ) +
  labs(
    x = "Paleo longitude (°)",
    y = "Paleo latitude (°)",
  )

ggplot() +
  geom_sf(
    data = coastlines_sf,
    color = "grey40",
    linewidth = 0.4,
    fill = NA
  ) +
  geom_point(
    data = data_all,
    aes(x = paleolong, y = paleolat, fill = Age),
    shape = 21,
    color = "black",
    size = 6
  ) +
  coord_sf(crs = 4326) +
  theme_minimal() +
  scale_fill_manual(
    name = "Âge géologique",
    values = cols_age
  ) +
  labs(
    x = "Paleo longitude (°)",
    y = "Paleo latitude (°)"
  )
