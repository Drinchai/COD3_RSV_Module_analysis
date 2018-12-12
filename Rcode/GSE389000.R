# Setup environment
rm(list=ls())
setwd("")
## install the dependencies (required packages)
### CRAN
required.packages <- c("gplots")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
### biocLite
source("http://bioconductor.org/biocLite.R")
required.biocLite.packages <- c("GEOquery")
missing.biocLite.packages <- required.biocLite.packages[!(required.biocLite.packages %in% installed.packages()[,"Package"])]
if(length(missing.biocLite.packages)) biocLite(missing.biocLite.packages)
## load packages into the memory
lapply(required.packages, library, character.only = TRUE)
lapply(required.biocLite.packages, library, character.only = TRUE)

#GET GSE soft file and matrix
GSE38900 <- getGEO("GSE38900", GSEMatrix=FALSE)

#sample names
names(GSMList(GSE38900 ))

#platforms used in this GSE
names(GPLList(GSE38900 ))
GSM.platforms <- lapply(GSMList(GSE38900 ),function(x) {Meta(x)$platform}) 
data.frame(GSM.platforms)

#example of an GSM experession vector 
Table(GSMList(GSE38900 )[[1]])[1:100,]

#example of gene anotation data from GPL of GSM 1
Probe.anotaion.table <- Table(GPLList(GSE38900)[[1]])[,c(1,2,4,6,11,12,10)]
rownames(Probe.anotaion.table) <- Probe.anotaion.table$ID
#Probeset extrated from GPL of GSM 1 
probesets <- as.character(Probe.anotaion.table$ID)
probesets <- Table(GPLList(GSE38900)[[1]])$ID
#creating the expression matrix ordered by the GPL order of probes
data.matrix <- do.call('cbind',lapply(GSMList(GSE38900 ),function(x) {
  tab <- Table(x)
  mymatch <- match(probesets,tab$ID_REF)
  return(tab$VALUE[mymatch])
}))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
# log2 trnasform the data # Data is already log2 transformed
#data.matrix.test <- log2(data.matrix)
# Add the probe names as row names

rownames(data.matrix) <- probesets
data.matrix <- data.matrix[complete.cases(data.matrix), ]
data.matrix[1:5,]
data.matrix[data.matrix<0]=1

###match probe and gene names
ProbeID.GSE38900 <- Probe.anotaion.table[which(rownames(Probe.anotaion.table)%in%rownames(data.matrix)),]
rownames(ProbeID.GSE38900)==rownames(data.matrix)
rownames(data.matrix) <- ProbeID.GSE38900$ILMN_Gene


#save to text file
write.csv (data.matrix,file="./GSE38900.epxr.matrix.txt")
write.csv (Probe.anotaion.table,file = "./GSE38900.probe.annot.table.txt")

#Anotation data
GSE38900.clinical  <- getGEO("GSE38900",GSEMatrix=TRUE)
Phenotipic.data <- pData(GSE38900.clinical[[2]])
Phenotipic.characteristics <- Phenotipic.data[,grep(x = colnames(Phenotipic.data), pattern = "characteristics")]
Phenotipic.characteristics <- data.frame(lapply(Phenotipic.characteristics, as.character), stringsAsFactors=FALSE)

#Fix colnames
fix.col.names <- t(as.data.frame(strsplit(x = as.character(Phenotipic.characteristics[1,]),split = ":")))
colnames(Phenotipic.characteristics) <- fix.col.names[,1]

#Fix data
for (i in 1:ncol(Phenotipic.characteristics)) {
  x<-colnames(Phenotipic.characteristics[i])
  x<- paste0(x,": ")
  print (x)
  x<-gsub(x=x ,pattern = "\\(",replacement = "\\\\(")
  x<-gsub(x=x ,pattern = "\\)",replacement = "\\\\)")
  Phenotipic.characteristics[,i] <- gsub(x = Phenotipic.characteristics[,i],pattern = x,replacement = "")
}

rownames(Phenotipic.characteristics) <- Phenotipic.data$geo_accession
sample.info1 <- Phenotipic.characteristics

sample.info2 <- Phenotipic.characteristics

sample.info <- rbind(sample.info1,sample.info2)

#Anotation data


save(data.matrix,sample.info, file = "./GSE38900_data.matrix.Rdata")

# preparing data ##
load("./GSE38900_data.matrix.Rdata")
GSE38900.dat.log2 <- log(data.matrix,2)  

rownames(data.matrix)


GSE38900.dat.log2.ann <- data.frame(Symbol = row.names(GSE38900.dat.log2), GSE38900.dat.log2)
GSE38900.dat.log2.ann <- data.frame(GSE38900.dat.log2)

### Prepare expression matrix with module list
Table.mod <- data.frame(GSE38900.dat.log2.ann[c(1:14168),])
Table.mod[,] <- NA
colnames(Table.mod) = gsub(colnames(Table.mod),pattern = "X",replacement = "")
New.mod.table <- as.data.frame(cbind(Modulelist_Gen3, Table.mod))

##Add expression value to moduleG2
module.genes <- unique(New.mod.table$Gene)
"%nin%" <- function (x, table) match(x, table, nomatch = 0L) == 0L # for exclude non-match RNAseq and Module
for (i in module.genes) {
  print(i)
  expression.gene <- as.numeric(GSE38900.dat.log2.ann[i,])
  samples = gsub(colnames(GSE38900.dat.log2.ann),pattern = "X",replacement = "")
  New.mod.table[New.mod.table$Gene == i,samples] <- expression.gene
}

GSE38900.mod <- New.mod.table[complete.cases(New.mod.table), ]    # exclused NA
rownames(GSE38900.mod) <- GSE38900.mod$Module_gene
GSE38900.mod.func.Gen3 <- GSE38900.mod[,c(1:6)]
GSE38900.mod.Gen3.log2 <- GSE38900.mod[,-c(1:5)]
rownames(GSE38900.mod.Gen3.log2) <- GSE38900.mod.Gen3.log2$Module_gene
GSE38900.mod.Gen3.log2$Module_gene <- NULL

save(GSE38900.mod.func.Gen3,GSE38900.mod.Gen3.log2,file = "./GSE38900_module_matrix.Rdata")

####Stat analysis####

## Clean up data
library("gtools")
#############################################
# Statistic analysis ##
############################################ 
GSE38900.mod.Gen3.log2 <- as.matrix(GSE38900.mod.Gen3.log2)
## prepare entry table 
tt_pval <- GSE38900.mod.Gen3.log2[,1,drop=FALSE]
colnames(tt_pval) <- c("pvalue")
tt_pval[,1] <- NA

for (k in 1:nrow(GSE38900.mod.Gen3.log2)) {
  signature = rownames(GSE38900.mod.Gen3.log2)[k]
  test.table <- sample.info 
  test.table$scores <- GSE38900.mod.Gen3.log2[k,]
  T2 <- test.table[test.table$`sample group` %in% c("Influenza_LRTI"),]
  T1 <- test.table[test.table$`sample group` %in% c("Healthy"),]
  if(T1$scores == T2$scores){
    tt_pval[k,] = 1
  }else{
    tt_pval[k,] <- t.test(x =T1$scores,y=T2$scores,paired = FALSE)$p.value
  }
}

pRSVvscontrol <- data.frame(tt_pval)


####calculate Fold change
GSE38900.mod.Gen3.raw <- 2^GSE38900.mod.Gen3.log2

############
library("gtools")
FC.group <- GSE38900.mod.Gen3.raw[,1,drop=FALSE]
colnames(FC.group) <- c("Foldchange")
FC.group[,1] <- NA

for (k in 1:nrow(GSE38900.mod.Gen3.raw)) {
  signature = rownames(GSE38900.mod.Gen3.raw)[k]
  test.table <- sample.info 
  test.table$scores <- GSE38900.mod.Gen3.raw[k,]
  T2 <- test.table[test.table$`sample group` %in% c("Influenza_LRTI"),]
  T1 <- test.table[test.table$`sample group` %in% c("Healthy"),]
  FC.group[k,] <- foldchange(mean(T2$scores),mean(T1$scores))
}

FCRSVvsControl <- data.frame(FC.group)

#### time point comparison ##### Group plot ###
#logical check ##
Group.up <- (FCRSVvsControl > 1)+(pRSVvscontrol < 0.05) == 2          # TRUE Up gene, Both TRUE

Group.down <- (FCRSVvsControl < -1) + (pRSVvscontrol < 0.05) == 2      # TRUE down gene, Both TRUE


################################################
Gene.matrix <- GSE38900.mod.func.Gen3[rownames(Group.up),]  
Gene.matrix$Module <- as.character(Gene.matrix$Module)

#####UP GENE######
up.mods.group <- data.frame(Module=NA,RSVvsControl=0,genes=0)

for (i in 1:length(unique(Gene.matrix$Module))){                                    # length of module
  module <- unique(Gene.matrix$Module)[i]                                           # look for only unique module
  sums <- colSums(Group.up[Gene.matrix$Module==module,1,drop=FALSE])              # sum upgene of each column by module 
  genes <- nrow(GSE38900.mod.func.Gen3[GSE38900.mod.func.Gen3$Module==module,])        # sum number of gene in each module
  up.mods.group <- rbind(up.mods.group,c(module,sums,genes))                               # paste result into a new fake table
}
up.mods.group <-up.mods.group[-1,]
rownames(up.mods.group) <- up.mods.group$Module
up.mods.group$Module <- NULL
up.mods.group.cal <- up.mods.group
up.mods.group <- as.data.frame(lapply(up.mods.group, as.numeric))                          # convert data frame to be numberic
up.mods.group <- (up.mods.group/up.mods.group$genes)*100 
rownames(up.mods.group) <-rownames(up.mods.group.cal)
up.mods.group <- up.mods.group[,-2,drop=FALSE]
#####DOWN GENE#######
down.mods.group <- data.frame (Module=NA,RSVvsControl=0,genes=0) # create a new blank table

for (i in 1:length(unique(Gene.matrix$Module))){
  module <- unique(Gene.matrix$Module)[i]
  sums <- colSums (Group.down[Gene.matrix$Module==module,1,drop=FALSE])
  genes <- nrow(GSE38900.mod.func.Gen3[GSE38900.mod.func.Gen3$Module==module,])
  down.mods.group <- rbind(down.mods.group,c(module,sums,genes))
}
down.mods.group<-down.mods.group[-1,]

rownames(down.mods.group) <- down.mods.group$Module
down.mods.group$Module <- NULL
down.mods.group.cal <- down.mods.group
down.mods.group <- as.data.frame(lapply(down.mods.group, as.numeric))
down.mods.group <- (down.mods.group/down.mods.group$genes)*100
rownames(down.mods.group) <- rownames(down.mods.group.cal)
down.mods.group <- down.mods.group[,-2,drop=FALSE]


########## DISPLAY DATA > 15 %
up.mods.group[up.mods.group < 15] <- 0
down.mods.group[down.mods.group < 15] <- 0

######## PREPARE DATA FOR PLOT ####
## Prepare data for ploting ## 
res.mods.group <- up.mods.group[,1,drop=FALSE]             ## prepare a new matrix for new data
res.mods.group[,1] <- NA                       ## Emtry matrix
colnames(res.mods.group) <- c("RSVvsControl")

for (i in 1: nrow(up.mods.group)){
  for (j in 1:ncol(up.mods.group)){
    #print (paste0 ("row:",i,"column:",j))
    up = up.mods.group[i,j]
    down = down.mods.group[i,j]
    if (up > down) {
      res = up
    }
    if (down > up){
      res = -down
    }
    if (up == down){
      res = 0
    }
    res.mods.group[i,j] = res
  }
}

################################################
## GRID PLOT
################################################

## prepared cluter position
ModuleG3_ann_grid <- ModuleG3_ann_grid[rownames(res.mods.group),]

rownames(res.mods.group)==rownames(ModuleG3_ann_grid)
rownames(res.mods.group) <- ModuleG3_ann_grid$Cluster_position


# creat new grid with all filtered cluster##
mod.group <- matrix(nrow=38,ncol=42)       
rownames (mod.group) <- paste0("A",c(1:38))
colnames (mod.group) <- paste0("",c(1:42))

## BY cluster
for (i in 1 : nrow(res.mods.group)){
  Mx <- as.numeric(gsub(x = strsplit (rownames(res.mods.group)[i],"\\.")[[1]][[1]],pattern = "A",replacement = ""))
  My <- as.numeric(strsplit (rownames(res.mods.group)[i],"\\.")[[1]][[2]])
  mod.group[Mx,My] <- res.mods.group$RSVvsControl[i] 
}

#cluster remove module that have only 1 module
mod.group <- mod.group[-c(9:14,19:23),]

library(ggplot2)
library(ggsignif)
library(reshape2)

melt_test <- melt(mod.group,id.var=c("row.names"))

################### MODULE MAP COLORs#########################################
dev.new()
ggplot(melt_test, aes(Var1, as.factor(Var2))) +
  geom_tile(color="#E6E6E6" , size = 0.2, fill= c(rep("white",27),
                                                  rep("white",27),
                                                  rep("white",9),"#E6E6E6","white","#E6E6E6",rep("white",15),
                                                  rep("white",9),"#E6E6E6","white","#E6E6E6",rep("white",15),
                                                  rep("white",9),"#E6E6E6","white","#E6E6E6",rep("white",14),"#E6E6E6",
                                                  rep("white",9),rep("#E6E6E6",3),rep("white",3),"#E6E6E6",rep("white",2),"#E6E6E6",rep("white",3),"#E6E6E6",rep("white",3),"#E6E6E6",
                                                  rep("white",6),"#E6E6E6",rep("white",2),rep("#E6E6E6",3),rep("white",3),rep("#E6E6E6",2),"white","#E6E6E6",rep("white",3),"#E6E6E6","white","#E6E6E6","white","#E6E6E6",
                                                  rep("white",6),"#E6E6E6",rep("white",2),rep("#E6E6E6",3),rep("white",3),rep("#E6E6E6",2),"white","#E6E6E6",rep("white",3),"#E6E6E6","white","#E6E6E6","white","#E6E6E6",
                                                  rep("white",6),rep("#E6E6E6",2),"white",rep("#E6E6E6",3),rep("white",3),rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white","#E6E6E6","white","#E6E6E6",
                                                  rep("white",6),rep("#E6E6E6",2),"white",rep("#E6E6E6",3),rep("white",3),rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white","#E6E6E6","white","#E6E6E6",
                                                  rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",2),"white",rep("#E6E6E6",3),rep("white",3),rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white","#E6E6E6","white","#E6E6E6",
                                                  rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",7),rep("white",2),rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
                                                  rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",8),"white",rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
                                                  rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",8),"white",rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
                                                  rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",8),"white",rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
                                                  rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
                                                  rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
                                                  rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
                                                  rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
                                                  rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
                                                  rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
                                                  rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",15),"white",rep("#E6E6E6",5),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",15),"white",rep("#E6E6E6",5),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",21),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",21),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",21),
                                                  rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",21),
                                                  rep("white",2),rep("#E6E6E6",25),
                                                  rep("white",2),rep("#E6E6E6",25),
                                                  "#E6E6E6","white",rep("#E6E6E6",25),
                                                  "#E6E6E6","white",rep("#E6E6E6",25),
                                                  "#E6E6E6","white",rep("#E6E6E6",25),
                                                  "#E6E6E6","white",rep("#E6E6E6",25)))+
  geom_point(aes(colour=value,size=1))+   
  ylab("") +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0))+
  scale_color_gradient2(low = "blue", mid="white", high = "red", limits=c(-100,100),na.value = "#E6E6E6", guide = "colourbar")+
  theme_light() +
  theme(panel.grid.minor = element_line(colour="black", size=0.9))+
  coord_flip() + scale_x_discrete(limits = rev(levels(melt_test$Var1))) +
  theme(panel.border = element_rect(color = "black",size = 0.5),
        axis.text.x = element_text(colour="black",size=9,angle=0,hjust=0.5,vjust=2,face="plain"),
        axis.text.y = element_text(colour="black",size=9,angle=0,hjust=0.5,vjust=0.5,face="plain"))







