# The intended purpose is to calculate the area between curves from Hydrogen/deuterium
#   exchange (HDX) data for two experimental states. The following steps should occur:
# 1. Enter experimental details: protein, labels for State1 and State2, timepoints, number of observations
# Note that the present calculation assumes the same number of observations for all time points for State 1 and State 2 and that there are no missing values.
# 2. Read in peptide info: PeptideID, Start, End, z, Sequence.
# 3. Read in deuterium incorporation % values and convert to remaining hydrogen fraction.
# 4. Read in associated standard deviations (SDs) and convert to hydrogen fraction SDs.
# 5. Calculate the difference of the means (Y-X), the pooled std.dev.(Sxy), and p values
# 6. Calculate the area (Abc),the standard deviation (Sbc), and p values
# 7. Output results as a dataframe.
# 8. make plots?

#' Area between exchange curves.
#'
#' @param proteinName       A String.
#' @param state0Name        A String.
#' @param state1Name        A String.
#' @param observations      A Number or a List of numbers.
#' @param timeList          A List of numbers. Each is the number of seconds a batch has been exposed to deuterium.
#' @param peptidesFile      A String.
#' @param state0DFracFile   A String.
#' @param state1DFracFile   A String.
#' @param state0SDFile      A String.
#' @param state1SDFile      A String.
#' @return myOutput         meanX, meanY, DifMeans, Sdm, pvaluedm, Areabc, Sbc, pvaluebc
#'          meanX           mean remaining hydrogen fractions of state 0.
#'          meanY           mean remaining hydrogen fractions of state 1.
#'          DifMeans        Difference of the mean hydrogen remaining between the two states.
#'          Sdm             The standard deviation of the difference of means.
#'          pvaluedm        The p-value of the difference of means.
#'          Areabc          Area between the exchange curves.
#'          Sbc             Standard deviation of the curve area.
#'          pvaluebc        The p-value of the area between the curves.
#' @examples
#' abec("PPARg", "apo", "Rosi", 4, c(1, 30, 60, 900, 3600), "peptides.txt", "apo.txt", "Rosi.txt",
#'      "apoSD.txt", "RosiSD.txt")
#' @export
abec <- function(proteinName, state0Name, state1Name, observations, timeList, peptidesFile, state0DFracFile, state1DFracFile, state0SDFile, state1SDFile) {

# 1. Enter experimental details
# Enter the name of the protein
myprotein <- proteinName
# Enter the labels for the two states
myLabel1 <- state0Name
myLabel2 <- state1Name
# Enter the number of observations per timepoint.
# In general it could vary for each timepoint. Needed improvement here.
nobs <- observations
# Enter the timepoints in seconds
times <- timeList
ntimes <- length(times)
ntimesp <- ntimes+1 # useful
# TODO: %s/(ntimes+1)/ntimesp/g

# 2. Read in Peptide info
myPeptides <- utils::read.table(peptidesFile,header=TRUE)
mynrow <- nrow(myPeptides)
PeptideID <- myPeptides[["PeptideID"]]

# outputfile
myOutput <- myPeptides
#TODO: replace with creating a new, blank file?

# 3. Read in D% values and convert to remaining hydrogen fractions.
myDup1 <- utils::read.table(state0DFracFile,header=TRUE)
myDup2 <- utils::read.table(state1DFracFile,header=TRUE)

# Setting up dataframes for remaining hydrogen fractions in State X and State Y
# Probably it would be better just to have these as matrices
myX <- myDup1; myY <- myDup2
# TODO: Need to set up automated generation of headings for myX and myY DFs
colnames(myX) <- c("PeptideID","apo_1s","apo_30s","apo_60s","apo_900s","apo_3600s")
myY <- myX;colnames(myY) <- c("PeptideID","Rosi_1s","Rosi_30s","Rosi_60s","Rosi_900s","Rosi_3600s")
myX[,2:6] <- 1-.01*myDup1[,2:6];myY[,2:6] <- 1-.01*myDup2[,2:6] # Calculating remaining hydrogen fractions

# 4. Read in SDs and convert to hydrogen fraction SDs.
myD_SD1 <- utils::read.table(state0SDFile,header=TRUE)
myD_SD2 <- utils::read.table(state1SDFile,header=TRUE)
mySx <- myD_SD1; mySy <- myD_SD2 #Creating DFs for Std. deviations as hydrogen fractions
# Probably it would be better just to have these as matrices
mySx[,2:6] <- .01*mySx[,2:6];mySy[,2:6] <- .01*mySy[,2:6]
colnames(mySx) <- c("PeptideID","Sx_1s","Sx_30s","Sx_60s","Sx_900s","Sx_3600s")
colnames(mySy) <- c("PeptideID","Sy_1s","Sy_30s","Sy_60s","Sy_900s","Sy_3600s")

# 5. Calculate the means, difference of the means (Y-X), the pooled std.dev.(Sdm), and p values
myXdata <- as.matrix(myX[,2:(ntimes+1)]); myYdata <- as.matrix(myY[,2:(ntimes+1)])
myOutput$meanX <- rowMeans(myXdata); myOutput$meanY <- rowMeans(myYdata)
myOutput$DifMeans <- myOutput$meanY-myOutput$meanX # Difference of means

# values currently are standard deviations.
# Square to get variances.
myXvar <- as.matrix(mySx[,2:ntimesp]);myYvar <- as.matrix(mySy[,2:ntimesp])
myXvar <- myXvar^2;myYvar <- myYvar^2;
colnames(myXvar) <- c("Varx_1s","Varx_30s","Varx_60s","Varx_900s","Varx_3600s")
colnames(myYvar) <- c("Vary_1s","Vary_30s","Vary_60s","Vary_900s","Vary_3600s")
myOutput$Sdm <- sqrt(0.5*(rowMeans(myXvar)+rowMeans(myYvar)))

# Calculate p values for difference of means
myOutput$pvaluedm <- 2*stats::pt(abs(nobs*myOutput$DifMeans/(2*myOutput$Sdm)),2*nobs-2,lower.tail = FALSE, log.p = FALSE)

# 6. Calculate the area (Abc),the standard deviation (Sbc), and p values
# set up a matrix of weights derived from log time ratios. First a vector.
logtimeweights <- rep(0.0,ntimes)
logtimeweights[1] <- log(times[2]/times[1],10)
for (i in 2:(ntimes-1)) {
  logtimeweights[i] <- log(times[i+1]/times[i-1],10)
}
logtimeweights[ntimes] <- log(times[ntimes]/times[ntimes-1],10)

# Now the matrix
myweights <- rep(logtimeweights,each=mynrow); dim(myweights) <- c(mynrow,ntimes)

# Matrices with weighted differences at each timepoint, or weighted standard deviations
mywtdiffs <- (myYdata-myXdata)*myweights;
colnames(mywtdiffs) <- c("wtdif_1s","wtdif_30s","wtdif_60s","wtdif_900s","wtdif_3600s")
mywtvars <- (myYvar+myXvar)*myweights;
colnames(mywtvars) <- c("wtvar_1s","wtvar_30s","wtvar_60s","wtvar_900s","wtvar_3600s")

# Calculate Areabc and standard deviations from row sums
Areabc <- rep(0.0,mynrow); Sbc2 <- rep(0.0,mynrow)
myOutput$Areabc <- 0.5*rowSums(mywtdiffs) # Area between curves
myOutput$Sbc <- sqrt(0.5*rowSums(mywtvars)) # Calculated SD of Area between curves
myOutput$pvaluebc <- 2*stats::pt(abs(nobs*myOutput$Areabc/(2*myOutput$Sbc)),2*nobs-2,lower.tail = FALSE, log.p = FALSE)
myOutput
}
