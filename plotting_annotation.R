# List of package names
packages_to_install <- c("ggplot2", "ggrepel")

# Loop through the packages
for (package_name in packages_to_install) {
  # Check if the package is installed
  if (!require(package_name, character.only = TRUE)) {
    # If not installed, install the package
    install.packages(package_name, dependencies = TRUE)
    
    # Load the package after installation
    library(package_name, character.only = TRUE)
  } else {
    # If already installed, just load the package
    library(package_name, character.only = TRUE)
  }
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if an argument is provided
if (length(args) == 0) {
  stop("No directory path provided.")
}

# Extract the directory path
print(args[1])
print(args[2])
input_file <- args[1]
file_name <- paste(args[2], "pdf", sep=".")

gene_hits = read.table(input_file, stringsAsFactors=FALSE, header=TRUE)

myplot <- ggplot(gene_hits, aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene)) +
   geom_point(alpha=0.5) +
   geom_text_repel(aes(size=15),  show.legend = FALSE, colour='black', alpha = 1) +
   scale_size("Number of k-mers", range=c(1,10)) +
   scale_colour_gradient('Average MAF') +
   theme_bw(base_size=8) +
   ggtitle("sccmec") +
   xlab("Average effect size") +
   ylab("Maximum -log10(p-value)")+
   ylim(0,NA)+
   xlim(0,NA)
   
# Split the string based on '\'
split_string <- strsplit(input_file, "/")[[1]]

# Remove the last element
split_string <- split_string[-length(split_string)]

# Combine the remaining elements back into a string
result_string <- paste(split_string, collapse = "/")
print(result_string)
print(file_name)

# Allocate result file name
result_file <- file.path(result_string, file_name)
  
pdf(result_file)
print(myplot) 
dev.off() 