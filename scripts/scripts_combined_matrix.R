library(tximport)
library(tidytree)

# List only featureCounts files (exclude .summary)
files <- list.files(path = "/Users/hellenhiza/Desktop/Scribles/RNA_snake/results/featurecounts", pattern = "\\.txt$", full.names = TRUE)

###Remove the txt.summary files 
files <- files[!grepl("\\.summary$", files)]

# Read counts
count_list <- lapply(files, function(f) {
  df <- read.table(f, header = TRUE, row.names = 1)
  df <- df[, ncol(df), drop = FALSE]  # keep only counts column
  return(df)
})

# Merge into one count matrix
count_matrix <- do.call(cbind, count_list)
colnames(count_matrix) <- gsub(".txt", "", basename(files))

#Now count_matrix is ready in memory
head(count_matrix)

###can also write it out but i am using analysis so will keep it 
write.csv(count_matrix, "count_matrix.csv") ##don't forget to add path 

#######################################
############CALCULATE TRASCRIPTS PER MILLION (TPM)
###Here using gene length from feature counts 
# Read first file to grab gene IDs and Length
first <- read.table(files[1], header = TRUE, row.names = 1)
gene_length <- first$Length / 1000   # convert bp → kb

#######RPK= raw counts/(gene_length_in_kb) then TPM= (RPK/sum(rpk)) x1e6
tpm_matrix <- apply(count_matrix, 2, function(x) {
  rpk <- x / gene_length
  rpk / sum(rpk) * 1e6
})

# Both matrices are in memory
head(count_matrix)
head(tpm_matrix)


##################

