#Initialization and restore data from check point
load("../code/FaroeQTL_savepoint_2.RData")
library(qtl)
require(car)
require(gridExtra)
require(ggplot2)
library(QTLRel)
library(abind)

#Sourcing mapping functions from lgsmfear in Parker et al., 2014, Genetics. Refer to:
#https://github.com/pcarbo/lgsmfear

Palmer <- dir<-"/external_data/Palmer_LG_J_SM_J"
source(paste0(Palmer <- dir,"lgsmfear/code/misc.R"))
source(paste0(Palmer <- dir,"lgsmfear/code/read.data.R"))
source(paste0(Palmer <- dir,"lgsmfear/code/data.manip.R"))
source(paste0(Palmer <- dir,"lgsmfear/code/mapping.tools.R"))

#Cleaning up and labelling data columns
colnames(pruned.new.data$pheno)[10]<-"GrowthModel_µ/MaxRate"
colnames(pruned.new.data$pheno)[32]<-'powerTransformed_GrowthModel_mu_MaxRate'
colnames(pruned.new.data$pheno)[51]<-'zScore_GrowthModel_mu_MaxRate'

#The genetic map needs to be reformatted for QTLrel
map<-pull.map(pruned.new.data)

amap<-matrix(NA, nrow=sum(nmar(pruned.new.data))-94-5, ncol=6)
names(amap)<-c("snp","chr","dist","refSNP","pos","distold")
amap<-data.frame(amap)
colnames(amap)<-c("snp","chr","dist","refSNP","pos","distold")

sum<-1
astr<-vector()
asnp<-vector()
adist<-vector()
apos<-vector()
chromosomes<-names(map)
for (chr in names(map)[1:19]) {
    print(chr);
    a<-rep(chr, length(map[[chr]]))
    astr<-c(astr, a)
    asnp<-c(asnp,names(map[[chr]]))
    adist<-c(adist, as.vector(map[[chr]]))
    apos<-c(apos, as.numeric(sapply(strsplit(names(map[[chr]]), "[:]"), "[[", 2)))
}

amap[,2]<-astr
amap[,4]<-asnp
amap[,1]<-asnp
amap[,3]<-adist
amap[,6]<-adist
amap[,5]<-apos

amap$chr<-as.factor(amap$chr)
amap$dist<-as.numeric(amap$dist)
str(amap$chr)
str(amap$dist)
G<-pruned.new.data$geno

#Reformatting the genotype files
ageno<-matrix(nrow=nind(pruned.new.data));
sum<-1
for (chr in names(map)[1:19]) {
    print(chr);
    atemp<-G[[chr]]$data
    print(str(atemp))
    atemp<-as.matrix(atemp)
    ageno<-cbind(ageno, atemp)
   sum<-sum+length(map[[chr]])
}

# Initialize storage for QTL mapping results (gwscan), parameter
# estimates of variance components from QTLRel analysis (vcparams),
# parameter estimates of additive QTL effects for individual markers
# (additive) and dominance QTL effects (dominance), and permutation
# tests (perms).
gwscan    <- list(F2.qtl   = NULL,
                  F2.rel   = empty.scanone(amap))
r         <- list(F2       = empty.scanone(amap))
additive  <- r
dominance <- r
pve       <- r
vcparams  <- list(F2 = NULL)
perms     <- list(F2 = NULL,F2.rel = NULL)

relatedness<-"markers"
covariates<-"sex"
num.perm<-1000

ind_retain<-!(is.na(pheno[,phenotype]))
pheno<-pheno[ind_retain,]
ageno<-ageno[ind_retain,]
cat("Calculating probabilities of missing genotypes.\n")
G<-ageno
G<-zero.na(G)
gp <- genoProb(G,amap,method = 'Kosambi', gr=2)

pheno <- transform(pheno,sex    = factor2integer(sex) - 1)
markers<-seq(1,length(amap$snp),1)
rows<-seq(1,length(pheno[,phenotype]), 1)

  cat("(a) QTL mapping with F2 cross using QTLRel\n")
#  F2.rows <- which(pheno$generation == "F2" &
#                   none.missing.row(pheno[cols]))
if (relatedness == "markers") {
    out.qtlrel <- map.cross.rr(pheno,G,amap,phenotype,covariates,gp, 
                        num.perm,rows=rows, markers=markers, verbose = TRUE)
}
  # Get the parameter estimates, QTL mapping results and estimated
  # distribution of maximum LOD scores under the null.
  gwscan$F2.rel[[phenotype]] <- out.qtlrel$gwscan$lod
  perms$F2.rel[[phenotype]]  <- out.qtlrel$perms
  additive$F2[[phenotype]]   <- out.qtlrel$gwscan$additive
  dominance$F2[[phenotype]]  <- out.qtlrel$gwscan$dominance
  pve$F2[[phenotype]]        <- out.qtlrel$gwscan$pve
  vcparams$F2                <- rbind(vcparams$F2,out.qtlrel$vcparams)


  #save.image("FaroeQTL_savepoint_3.RData")
