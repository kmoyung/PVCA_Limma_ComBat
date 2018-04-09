#Programmer:  Pierre R. Bushel / Modified by Krishna Sriram, Justin Chou and Kevin Moyung 
#Location: NIEHS / UCSD Dept of Pharmacology
#email  Bushel@niehs.nih.gov / ksriram@ucsd.edu
#Code: R
#Program name: pvca.R
#Date: May 26, 2009 / November, 2017

# Justin's version Jan 18th 2018 

#setwd("C:/Users/ksriram/Documents/rwkd") #change as appropriate
#setwd("C:/Users/yogi.sriram/Documents/rwkd")
#setwd("C:/Users/ksriram/Documents/Kevin Files/PVCA Test on ACC")
#setwd("C:/Kevin/PVCA and Combat on OV CPMs")
setwd("C:/Kevin/LUAD")


########## load libraries ##########
library(lme4)
library(RSpectra)
####### Edit these variables according to user defined parameters and the path to your data and data files names ##########

#myPath <- "../../Data/"
theGene_expression_file <- "LUAD_CPMs_SourceCenter_Adjusted_Limma.txt"
theExperiment_data_file <- "LUAD_BatchEffects.txt"
pct_threshold = .9 # Amount of variability desired to be explained by the principal components.  Set to match the results in book chapter and SAS code.  User can adjust this to a higher (>= 0.8) number but < 1.0

### In addition, be sure to modify the mixed linear model by adding the appropriate random effects terms in the model

################################################

theGEDFilePath = paste(theGene_expression_file, sep="")
theExpDataFilePath = paste(theExperiment_data_file, sep="")

########## Load data ##########

theDataMatrix <- read.delim(theGEDFilePath, row.names = 1, header = TRUE, sep = "\t")
dataRowN <- nrow(theDataMatrix)
dataColN <- ncol(theDataMatrix)

########## Center the data (center rows) ##########

theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, ncol = dataColN)
theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
theDataMatrixCentered = t(theDataMatrixCentered_transposed)
exp_design <- read.delim(theExpDataFilePath, sep = "\t", header = TRUE, row.names = 1)
expDesignRowN <- nrow(exp_design)
expDesignColN <- ncol(exp_design)
myColNames <- names(exp_design)

########## Compute correlation matrix ##########

theDataCor <- cor(theDataMatrixCentered)

########## Obtain eigenvalues ##########

# Use RSpectra's eigen decomposition for more efficiency. Be sure to change '12' to the
# correct dimensions of the correlation matrix when using other files.
#eigenData <- eigen(theDataCor)
eigenData <- eigs_sym(theDataCor, 41, which = "LM")

eigenValues = eigenData$values
ev_n <- length(eigenValues)
eigenVectorsMatrix = eigenData$vectors
eigenValuesSum = sum(eigenValues)
percents_PCs = eigenValues /eigenValuesSum 

########## Merge experimental file and eigenvectors for n components ##########

my_counter_2 = 0
my_sum_2 = 1
for (i in ev_n:1){
  my_sum_2  = my_sum_2 - percents_PCs[i]
  if ((my_sum_2) <= pct_threshold ){
    my_counter_2 = my_counter_2 + 1
  }
  
}
if (my_counter_2 < 3){
  pc_n  = 3
  
}else {
  pc_n = my_counter_2 
}

# pc_n is the number of principal components to model

pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
mycounter = 0
for (i in 1:pc_n){
  for (j in 1:expDesignRowN){
    mycounter <- mycounter + 1
    pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]
    
  }
}

AAA <- exp_design[rep(1:expDesignRowN,pc_n),]

Data <- cbind(AAA,pc_data_matrix)

####### Edit these variables according to your factors #######

Data$Plate_ID <- as.factor(Data$Plate_ID)
Data$Source_Center <- as.factor(Data$Source_Center)
#Data$Center_ID <- as.factor(Data$Center_ID)
Data$The_28s_18s <- as.factor(Data$The_28s_18s)
Data$RIN_Score <- as.factor(Data$RIN_Score)


########## Mixed linear model ##########
### MIGHT NEED TO ALTER THE MATH FOR effects_n IF USING DIFFERENT # OF VARIABLES IN EXP DATA ###
op <- options(warn = (-1)) 
effects_n = (expDesignColN - 2) + ((expDesignColN - 2)*(((expDesignColN - 2)-1)))/2 + 1
#effects_n = (expDesignColN) + 1 #hardcoding
randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)


for (i in 1:pc_n){
  y = (((i-1)*expDesignRowN)+1)
  #randomEffects <- (summary(Rm1ML <- lmer(pc_data_matrix ~ (1|Time) + (1|Treatment) + (1|Batch) + (1|Time:Treatment) + (1|Time:Batch) + (1|Treatment:Batch), Data[y:(((i-1)*expDesignRowN)+expDesignRowN),], REML = TRUE, control=lmerControl(maxIter = 1000000, msMaxIter=1000000, singular.ok=TRUE,tolerance=1e-4, returnObject=TRUE),verbose = FALSE, na.action = na.omit))@REmat)
  #randomEffects <- (summary(Rm1ML <- lmer(pc_data_matrix ~ (1|Time) + (1|Treatment) + (1|Batch) + (1|Time:Treatment) + (1|Time:Batch) + (1|Treatment:Batch), Data[y:(((i-1)*expDesignRowN)+expDesignRowN),], REML = TRUE, verbose = FALSE, na.action = na.omit))@REmat)
  #randomEffects <- (VarCorr(Rm1ML <- lmer(pc_data_matrix ~ (1|RIN_Score) + (1|Gender) + (1|rnaRatio) + (1|RIN_Score:Gender) + (1|RIN_Score:rnaRatio) + (1|Gender:rnaRatio), Data[y:(((i-1)*expDesignRowN)+expDesignRowN),], REML = TRUE, verbose = FALSE, na.action = na.omit)))
  randomEffects <- (VarCorr(Rm1ML <- lmer(pc_data_matrix ~ (1|Plate_ID) + (1|The_28s_18s) + (1|RIN_Score) + (1|Source_Center) + (1|The_28s_18s:RIN_Score) + (1|Plate_ID:The_28s_18s) + (1|Plate_ID:RIN_Score) + (1|Source_Center:Plate_ID) + (1|Source_Center:The_28s_18s) + (1|Source_Center:RIN_Score), Data[y:(((i-1)*expDesignRowN)+expDesignRowN),], REML = TRUE, verbose = FALSE, na.action = na.omit)))
  randomEffectsDF <- as.data.frame(randomEffects)
  
  for (j in 1:effects_n){
    randomEffectsMatrix[i,j] = as.numeric(randomEffectsDF[j,5])
  }
  
}

effectsNames <- randomEffectsDF[,1]

########## Standardize Variance ##########

randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
for (i in 1:pc_n){
  mySum = sum(randomEffectsMatrix[i,])
  for (j in 1:effects_n){
    randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum	
  }
}

########## Compute Weighted Proportions ##########

randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
for (i in 1:pc_n){
  weight = eigenValues[i]/eigenValuesSum
  for (j in 1:effects_n){
    randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight
  }
}

########## Compute Weighted Ave Proportions ##########

randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
randomEffectsSums <-colSums(randomEffectsMatrixWtProp)
totalSum = sum(randomEffectsSums)
randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, ncol = effects_n)

for (j in 1:effects_n){
  randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum 	
  
}
bp <- barplot(randomEffectsMatrixWtAveProp,  xlab = "Effects", ylab = "Weighted average proportion variance", ylim= c(0,1.1),col = c("blue"), las=2)

axis(1, at = bp, labels = effectsNames, xlab = "Effects", cex.axis = 0.5, las=2)
values = randomEffectsMatrixWtAveProp
new_values = round(values , 3)
text(bp,randomEffectsMatrixWtAveProp,labels = new_values, pos=3, cex = 0.8) # place numbers on top of bars 




