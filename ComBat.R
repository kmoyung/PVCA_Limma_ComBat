# Modified by Kevin Moyung
# Removes a known batch effect using the ComBat function in the SVA package (Bioconductor)

library(sva)
library(bladderbatch)
setwd("C:/Kevin/LUAD")
#setwd("C:/Justin")

# Read the CPM file and batch effects file
edata <- read.table("LUAD_CPM.txt", sep="\t", header=TRUE, row.names=1 )
pheno <- read.table("LUAD_BatchEffects_ComBat.txt", sep="\t", header=TRUE)
# Justin's data
#edata <- read.table("PDAC_CPM_ComBat.txt", sep="\t", header=TRUE, row.names=1 )
#pheno <- read.table("PAAD_Technical_Factors_ComBat.txt", sep="\t", header=TRUE)

# Select batch effect to be removed - CHANGE as needed
batch = pheno$Source_Center
#batch = pheno$The_28s_18s


# Create model matrix based on the batch effects file
mod = model.matrix(~1, data=pheno)

### Run ComBat ###

# parametric adjustment
combat_edata1 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)

# non-parametric adjustment, mean-only version
combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)

# reference-batch version, with covariates
combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3)

# Save new expression data batch removed
write.table(combat_edata1, file = "C:/Kevin/LUAD/LUAD_CPMs_SourceCenter_Adjusted.txt", append = TRUE, sep = "\t", row.names=TRUE, col.names=TRUE)
#write.table(combat_edata1, file = "C:/Justin/PDAC_CPM_28s_18s_Adjusted.txt", append = TRUE, sep = "\t", row.names=TRUE, col.names=TRUE)
