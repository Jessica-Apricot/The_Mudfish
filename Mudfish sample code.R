#Packages
##########install.packages(c("betapart", "ade4", "labdsv", "ape", "ggplot2", "vegan"))
library(betapart)
library(ade4)
library(labdsv)
library(ape)
library(ggplot2)
library(vegan)


sample_data = read.csv("presence_absence_matrix.csv")
# Remove first column (assumed taxa names or IDs) this is a test

dat <- sample_data[,-1]

# Transpose data so rows = sites, columns = taxa
dat <- t(dat)
colnames(dat) <- sample_data$Taxa

# Convert to data frame for easier manipulation
datfull <- as.data.frame(dat)
#Add forestcover vector
datfull$forestcover <- c(rep("Forested", 15), rep("Reforested", 15))

# Add site names as a column
datfull$site <- rownames(datfull)

# Make sure all columns (taxa) are numeric
for(i in 1:(ncol(datfull)-2)){
  datfull[, i] <- as.numeric(as.character(datfull[, i]))
}

# Convert back to matrix for beta diversity functions
dat_matrix <- as.matrix(datfull[, 1:(ncol(datfull)-2)])

# Calculate multi-site beta diversity metrics
beta_multi <- beta.multi(dat_matrix)
print(beta_multi)

# Calculate pairwise Sorensen dissimilarity matrix
soren <- beta.pair(dat_matrix, index.family = "sorensen")
Bsor <- soren$beta.sor

# Histogram of dissimilarities (optional)
hist(Bsor, breaks = 30, main = "Histogram of Sorensen distances")

# Principal Coordinates Analysis (PCoA)
pcoa <- pcoa(Bsor)

forestcover <- c(rep("Forested", 15), rep("Reforested", 15))

# Create data frame for plotting PCoA results

datpcoa <- data.frame(
  p1 = pcoa$vectors[, 1],
  p2 = pcoa$vectors[, 2],
  p3 = pcoa$vectors[, 3],
  p4 = pcoa$vectors[, 4],
  forestcover = forestcover
)

datpcoa$label <- rownames(datpcoa)  # Add labels for sites

#PCoA plot with labels

ggplot(datpcoa, aes(x = p1, y = p2)) +
  geom_point(aes(fill = forestcover), shape = 21, size = 6) +
  geom_text(aes(label = label), vjust = -1, size = 3) +
  theme_classic() +
  labs(x = "PCoA 1", y = "PCoA 2", fill = "") +
  scale_fill_manual(
    labels = c("Forested", "Reforested"),
    values = c("Forested" = "forestgreen", "Reforested" = "goldenrod")
  ) +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black")
  )

# PERMANOVA test for differences between forestcover groups

permanova_result <- adonis2(Bsor ~ forestcover, data = datpcoa, permutations = 9999)

print(permanova_result)





### Testing for species specific shifts ###

dat2=Sample_data[,-1]

#no. forested and deforested sites in which each taxa occurs

dat2$occurforest=rep(NA,nrow(dat2))

dat2$occurreforest=rep(NA,nrow(dat2))

for(i in 1:nrow(dat2)){
  dat2$occurforest[i]=sum(dat2[i,1:15])  #Adjust this based on the number of you have of each category##
  50
  dat2$occurreforest[i]=sum(dat2[i,16:30]) } #Adjust this based on the number of you have of each category##

dim(dat2)

datz=dat2[,31:32]

datz$totaloccurrence=datz$occurforest+datz$occurreforest

str(datz)

rownames(datz) = data$Taxa

datz=datz[which(datz$totaloccurrence>5),] **#May need to lower this number if you don't have many sites##**
  
  #p-value vector
  
  datz$pval=rep(NA,nrow(datz))

#Fisher’s exact test

for(i in 1:nrow(datz)){
  
  test=fisher.test(matrix(c(datz$occurforest[i], 15-datz$occurforest[i], #Adjust this based on the number of you have of each category##**
                            
 datz$occurreforest[i], 15-datz$occurreforest[i]), ncol=2)) **#Adjust this based on the number of you have of each category##**
    
    datz$pval[i]=test$p.value
    
}

#number of significant taxa (unadjusted p-value)

length(which(datz$pval<0.1))

datz[which(datz$pval<0.1),]



#Benjamini-Hochberg correction

datz$adjustedpval=p.adjust(datz$pval, method="fdr")

#number of significant taxa after p-value adjustment for multiple comparisons

length(which(datz$adjustedpval<0.05))


datz[which(datz$adjustedpval<0.05),]






###to combine files to binary presence/absence
# Load libraries
library(dplyr)
library(tidyr)
library(stringr)

# Set the folder with your CSV files
data_folder <- "/Users/jessdarnley/Library/CloudStorage/OneDrive-UniversityofOtago/Honours\ 2025/Lamprey/The_Mudfish/"

# List all CSV files in the folder
files <- list.files(data_folder, pattern = "\\.csv$", full.names = TRUE)

# Function to read each file and return a dataframe with taxa + location
read_location <- function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  # Assume first column is taxa names
  taxa_col <- names(df)[1]
  
  # Get location name from file (remove folder + extension)
  location <- str_remove(basename(file), [/.csv$]\\.csv$)
  
  # Create dataframe with taxa and presence (1)
  data.frame(Taxa = df[[taxa_col]],
             Location = location,
             Presence = 1)
}

# Read all files
all_data <- lapply(files, read_location) %>% bind_rows()

# Spread into presence/absence matrix
presence_absence <- all_data %>%
  pivot_wider(names_from = Location, values_from = Presence, values_fill = 0) %>%
  arrange(Taxa)

# View result
print(presence_absence)

# Save to CSV
write.csv(presence_absence, "presence_absence_matrix.csv", row.names = FALSE)
