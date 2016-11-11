# The intended purpose is to calculate the area between curves from Hydrogen/deuterium
#   exchange (HDX) data for two experimental states. The following steps should occur:
# 1. Enter experimental details: protein, labels for State1 and State2, timepoints, number of repetitions
# Note that the present calculation assumes the same number of repetitions for all time points for State 1 and State 2 and that there are no missing values.
# 2. Read in peptide info: PeptideID, Start, End, z, Sequence.
# 3. Read in deuterium incorporation % values and convert to remaining hydrogen fraction.
# 4. Read in associated standard deviations (SDs) and convert to hydrogen fraction SDs.
# 5. Calculate the difference of the means (Y-X), the pooled std.dev.(Sxy), and p values
# 6. Calculate the area (Abc), the standard deviation (Sbc), and p values
# 7. Output results as a dataframe.
# 8. make plots?

#' Area between exchange curves.
#'
#' @param proteinName       A String.
#' @param state0Name        A String.
#' @param state1Name        A String.
#' @param repetitions       A Number or an array.
#' @param timeList          A List of numbers. Each is the number of seconds a batch has been exposed to deuterium.
#' @param peptidesFile      A String.
#' @param state0DFracFile   A String.
#' @param state1DFracFile   A String.
#' @param state0SDFile      A String.
#' @param state1SDFile      A String.
#' @return result         state0mean, state1mean, mean_difs, Sdm, pvaluedm, Areabc, Sbc, pvaluebc
#'          state0mean           mean remaining hydrogen fractions of state 0.
#'          state1mean           mean remaining hydrogen fractions of state 1.
#'          mean_difs        Difference of the mean hydrogen remaining between the two states.
#'          Sdm             The standard deviation of the difference of means.
#'          pvaluedm        The p-value of the difference of means.
#'          Areabc          Area between the exchange curves.
#'          Sbc             Standard deviation of the curve area.
#'          pvaluebc        The p-value of the area between the curves.
#' @examples
#' abec("PPARg", "data/Peptides.txt", 4, c(1, 30, 60, 900, 3600), "apo", "Rosi",
#'  "data/apo.txt", "data/Rosi.txt", "data/apoSD.txt", "data/RosiSD.txt")
#' @export
abec <- function(protein_name, peptides_file,
                 repetitions, times_list,
                 state0name, state1name,
                 state0dfracfile, state1dfracfile,
                 state0sdfile, state1sdfile) {
times_length <- length(times_list)

# 2. Read in Peptide info
peptides_table <- utils::read.table(peptides_file, header = TRUE)
num_peptides <- nrow(peptides_table)

# outputfile
result <- peptides_table

# 3. Read in D% values and convert to remaining hydrogen fractions.
# Setting up dataframes for remaining hydrogen fractions in State X and State Y
# Probably it would be better just to have these as matrices
state0dfrac <- utils::read.table(state0dfracfile, header = TRUE)
state1dfrac <- utils::read.table(state1dfracfile, header = TRUE)

# Calculating remaining hydrogen fractions
state0dfrac[, 2:6] <- 1 - .01 * state0dfrac[, 2:6]
state1dfrac[, 2:6] <- 1 - .01 * state1dfrac[, 2:6]

# 4. Read in SDs and convert to hydrogen fraction SDs.
#Creating DFs for Std. deviations as hydrogen fractions
# Probably it would be better just to have these as matrices
state0sd <- utils::read.table(state0sdfile, header = TRUE)
state1sd <- utils::read.table(state1sdfile, header = TRUE)
state0sd[, 2:6] <- .01 * state0sd[, 2:6]
state1sd[, 2:6] <- .01 * state1sd[, 2:6]

# 5. Calculate the means, difference of the means (Y-X), the pooled std.dev.(Sdm), and p values
state0dfracdata <- as.matrix(state0dfrac[, 2:(times_length + 1)])
state1dfracdata <- as.matrix(state1dfrac[, 2:(times_length + 1)])
result$state0mean <- rowMeans(state0dfracdata)
result$state1mean <- rowMeans(state1dfracdata)
result$mean_difs <- result$state1mean - result$state0mean # Difference of means

# values currently are standard deviations.
# Square to get variances.
state0dfracvar <- as.matrix(state0sd[, 2:(times_length + 1)])
state1dfracvar <- as.matrix(state1sd[, 2:(times_length + 1)])
state0dfracvar <- state0dfracvar ^ 2
state1dfracvar <- state1dfracvar ^ 2
result$Sdm <- sqrt(0.5 * (rowMeans(state0dfracvar) + rowMeans(state1dfracvar)))

# Calculate p values for difference of means
result$pvaluedm <- 2 * stats::pt(abs(repetitions * result$mean_difs / (2 * result$Sdm)), 2 * repetitions - 2, lower.tail = FALSE, log.p = FALSE)

# 6. Calculate the area (Abc), the standard deviation (Sbc), and p values
# set up a matrix of weights derived from log time ratios. First a vector.
logtimeweights <- rep(0.0, times_length)
logtimeweights[1] <- log(times_list[2] / times_list[1], 10)
for (i in 2:(times_length - 1)) {
  logtimeweights[i] <- log(times_list[i + 1] / times_list[i - 1], 10)
}
logtimeweights[times_length] <- log(times_list[times_length] / times_list[times_length - 1], 10)

# Now the matrix
myweights <- rep(logtimeweights, each = num_peptides)
dim(myweights) <- c(num_peptides, times_length)

# Matrices with weighted differences at each timepoint, or weighted standard deviations
mywtdiffs <- (state1dfracdata - state0dfracdata) * myweights
mywtvars <- (state1dfracvar + state0dfracvar) * myweights


# Calculate Areabc and standard deviations from row sums
# Area between curves
result$Areabc <- 0.5 * rowSums(mywtdiffs)
# Calculated SD of Area between curves
result$Sbc <- sqrt(0.5 * rowSums(mywtvars))
result$pvaluebc <- 2 * stats::pt(abs(repetitions * result$Areabc / (2 * result$Sbc)), 2 * repetitions - 2, lower.tail = FALSE, log.p = FALSE)
result
}

# Generate a list of the form prefix_times (e.g. wtdif_1s, wtdif_30s...)
generate_time_names <- function(prefix, times_list) {
  result_list <- paste(prefix, "_", times_list, "s", sep = "")
  result_list
}
