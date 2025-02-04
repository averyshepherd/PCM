## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(phylocompmeth)

## -----------------------------------------------------------------------------
library(ape)

## -----------------------------------------------------------------------------
f <- "https://raw.githubusercontent.com/averyshepherd/phylocompmethods/main/Potential%20Nexus%20Files/T6506.nex"
tree <- read.nexus(f)

## -----------------------------------------------------------------------------
class(tree)
str(tree)

## -----------------------------------------------------------------------------
head(tree$edge)

## -----------------------------------------------------------------------------
head(tree$tip.label)

## -----------------------------------------------------------------------------
head(tree$Nnode)

## ----fig1, fig.height = 6, fig.width = 5--------------------------------------
plot(tree, cex = 0.6, type = "cladogram")

## -----------------------------------------------------------------------------
library(phytools)

## -----------------------------------------------------------------------------
anole.resids<-cbind(anole.data[,1],
    phyl.resid(anoletree,anole.data[,1,drop=FALSE],
    anole.data[,2:ncol(anole.data)])$resid)
colnames(anole.resids)[1]<-"SVL"

## -----------------------------------------------------------------------------
phylo.heatmap(anoletree,anole.resids,
    split=c(0.7,0.3),fsize=c(0.4,0.8,0.8),
    standardize=TRUE)

## -----------------------------------------------------------------------------
anole.pca <- pca(anoletree, anole.data)

## -----------------------------------------------------------------------------
print(anole.pca)
summary(anole.pca)

## -----------------------------------------------------------------------------
plot(anole.pca)

## -----------------------------------------------------------------------------
biplot(anole.pca)

## -----------------------------------------------------------------------------
# 10 species, total length of the tree is 100
World.tree<-pbtree(n=10,scale=100)
# setting labels for branches
World.tree$tip.label<-replicate(Ntip(World.tree),
    paste(sample(LETTERS,1),".",
    paste(sample(letters,round(runif(n=1,min=3,max=10))),
    collapse=""),
    sep=""))
# simulating latitudes and longitudes
lat<-fastBM(World.tree,sig2=10,bounds=c(-90,90))
long<-fastBM(World.tree,sig2=80,bounds=c(-180,180))
# group data
World<-cbind(lat,long)
for(i in 1:Ntip(World.tree)){
    ni<-sample(0:2,1)
    for(j in 1:ni){ 
        World<-rbind(World,c(World[i,1]+rnorm(n=1,sd=4),
            World[i,2]+rnorm(n=1,sd=4)))
        rownames(World)[nrow(World)]<-rownames(World)[i]
    }
}


## -----------------------------------------------------------------------------
head(World, 20)

## -----------------------------------------------------------------------------
obj <- phylo.to.map(World.tree, World, plot=FALSE)
# default settings
plot(obj, ftype="i")

# colors, font size, point size, etc. can be customized
cols<-setNames(sample(rainbow(n=Ntip(World.tree))),
    World.tree$tip.label)
plot(obj,colors=cols,ftype="i",fsize=0.8,cex.points=c(0.7,1.2))

## -----------------------------------------------------------------------------
# simulate data for the analysis
# phylogeny with 100 species
tree<-pbtree(n=100)
# transition rates
Q<-matrix(c(-2,1,1,
    1,-2,1,
    1,1,-2),3,3)
rownames(Q)<-colnames(Q)<-letters[1:3]
x<-as.factor(sim.history(tree,Q)$states)
# simulating phenotypic trait values for species for two different scenarios
y1<-fastBM(tree)
y2<-fastBM(tree,sig2=0.5)+as.numeric(x)

## -----------------------------------------------------------------------------
# y1
phenogram(tree,y1,ftype="off",col=make.transparent("blue",0.5))
tiplabels(pie=to.matrix(x,letters[1:3]),
    piecol=colorRampPalette(c("blue", "yellow"))(3),
    cex=0.4)

# y2
phenogram(tree,y2,ftype="off",col=make.transparent("blue",0.5))
tiplabels(pie=to.matrix(x,letters[1:3]),
    piecol=colorRampPalette(c("blue", "yellow"))(3),
    cex=0.4)

## -----------------------------------------------------------------------------
phyloANOVA(tree,x,y1)

phyloANOVA(tree,x,y2)

