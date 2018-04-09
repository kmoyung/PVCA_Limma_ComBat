# Kevin Moyung

library(edgeR)
#setwd("C:/Users/ksriram/Documents/rwkd")
setwd("C:/Kevin/LUAD")

mats <- read.table( file = "LUAD_BatchEffects_Limma_SourceCenter.txt" )
#mats[1] <- NULL

#plate.id <- mats[-c(1, 2), ]
#plates <- as.matrix(mats[-c(1), ])
#plateval <- as.vector(as.matrix(mats[-c(1, 3), ]))	
sourcecenter <- as.vector(as.matrix(mats[-c(1), ]))

raw.data <- read.table( file = "LUAD_CPM_Limma.txt" )
#raw.data2 <- raw.data[-c(1:88) ]

batch_var <- sourcecenter


counts <- raw.data[ , ]
#cds <- DGEList( counts )
#head(cds$counts)
#cds$samples

################ filter out the genes with >1 read per million in at least 3 samples

#cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
#dim( cds )    ########### the matrix will now be smaller due to filtered genes
#cds <- calcNormFactors( cds )
##### cds <- calcNormFactors( cds, c )
#cds$samples
# effective library sizes
#cds$samples$lib.size*cds$samples$norm.factors

#batch2 <- plates
#batch2 <- plateval

#logCPM <- cpm(cds, log=TRUE, prior.count=5)
logCPMc <- removeBatchEffect(counts, batch_var)

# Appends a 'c' to the end of the sample names for the corrected file
colnames(logCPMc) <- paste(colnames(logCPMc), "c", sep="")

#write.table(logCPM, file = "OV_logCPM1", sep="\t", col.names=NA)        ##### Rename this as needed
write.table(logCPMc, file = "LUAD_logCPM_SourceCenter_Adjusted_Limma.txt", sep="\t", col.names=NA)      ##### Rename this as needed



mdssave <- plotMDS(logCPM , col=as.numeric(cds$samples$group))
qq <- unlist(mdssave)
write.table(qq, file = "OV_uncorrected_MDS", sep="\t")      ##### Rename this as needed


mdssave <- plotMDS(logCPMc , col=as.numeric(cds$samples$group))
qq <- unlist(mdssave)
write.table(qq, file = "OV_corrected_MDS", sep="\t")      ##### Rename this as needed


raw.data2 <- cbind(logCPM, logCPMc) #### correct this syntax
#write.table(raw.data2, file = "OV_logCPMmerged", sep="\t", col.names=NA)
#raw.data2 = as.data.table(raw.data2)
counts2 <- 2^raw.data2

counts2 <- as.data.frame(counts2)

cds_new <- DGEList(counts2)
head(cds_new$counts2)
cds_new$samples


group <- c(rep("uncorrected", 404) , rep("corrected", 404))
mdssave <- plotMDS( cds_new , col=as.numeric(cds_new$samples$group))       


