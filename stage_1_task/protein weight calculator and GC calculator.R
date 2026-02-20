gc_calculator <- function(sequence) {
  
  g_count <- 0
  c_count <- 0
  sequence <- toupper(sequence)
  bases <- strsplit(sequence, "")[[1]]
  seq_length <- length(bases)
  for (base in bases) {
    if( base == 'C') {
      c_count <- c_count + 1
    } else if (base == 'G') {
      g_count <- g_count + 1
    }
  }
  gc_percent <- (g_count + c_count) * 100 /seq_length
  return(paste0('The percentage of GC is ', gc_percent,'%'))  
  
}


gc_calculator('GCAgTTAT')



protein_weight_calculator <- function(protein) {
# A function that calculates the weight of a protein in Kilodalton (Kda).
#  It returns 0 once an invalid amino acid is found in the sequence.

# Requesting for a valid protein sequence from a user  
protein <- readline("Kindly enter a valid protein sequence :   ")

# Initializing the weight of the protein
protein_weight_da <- 0

# Converting the protein sequence to an upper case
protein  <- toupper(protein)

# Splitting the protein string into individual characters(amino acid)
amino_acids <- strsplit(protein, '')[[1]]

# Creating a vector (lookup table) of amino acids and it's molecular weights(da)
data_dictionary <- c('A' = 89.09, 'R' =174.20,
                     'N' = 132.12, 'D' =133.10,
                     'C' =121.15, 'E' = 147.13,
                     'Q' = 146.15, 'G' = 75.07,
                     'H' = 155.16, 'I' = 131.18,
                     'L' = 131.18, 'K' = 146.19,
                     'M' = 149.21, 'F' = 165.19,
                     'P' = 115.13, 'S' = 105.09,
                     'T' = 119.12, 'W' = 204.23,
                     'Y' = 181.19, 'V' = 117.15)

# Loop through each amino acid in the protein sequence
for (amino_acid in amino_acids) {
  
# Check if the amino acid exists in the dictionary
  if (amino_acid %in% names(data_dictionary)) {
    
# Add the amino acid's weight to the total protein weight
    protein_weight_da <- protein_weight_da + data_dictionary[amino_acid]
    
# Convert the protein weight from dalton(da) to kilodalton(Kda) 
    protein_weight_kda <- protein_weight_da / 1000
    #  If an invalid character is found    
  } else {

# Stop the function and return an error message
    return('The protein is invalid, hence the weight is 0')
  }
  
}
# Return the final protein weight as a readable message
return(
  paste0('The value of the protein provided is ', protein_weight_kda, ' KiloDalton'))   
}


# Testing the function 
protein_weight_calculator()
chamo
