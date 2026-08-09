#******************************************************************************
#### 0. LIBRARIES ####
#******************************************************************************

library(dplyr)
library(openxlsx)
library(readxl)
library(tidyr)
library(stringr)
library(readr)
library(ggrepel)
library(circlize)
library(ggplot2)

#******************************************************************************
#### 1. DATA INPUT ####
#******************************************************************************

# NOTE: Uncomment the specific path that matches your current OS environment

# Dataset containing all publications related to glaucoma
all.gwas <- read.delim("~/data/combined_gwas.tsv")

# Metadata of the articles included in the screening process (Rayyan)
included.articles <- read.csv("~/data/bibliometric/rayyan_included.csv")

# GWAS data included for analysis
Snps.GPAA <- read_excel("~/data/combined_gwas_included.xlsx")
# Snps.GPAA <- read_excel("C:/Users/mloel/OneDrive/Desktop/glaucoma-gwas/data/combined_gwas_included.xlsx")

#******************************************************************************
#### 2. DATA DEPURATION ####
#******************************************************************************

# Retain only the GWAS data from the included articles
all.gwas.included <- all.gwas %>% 
  filter(PUBMEDID %in% included.articles$pubmed_id)

# Export the filtered GWAS dataset (Optional)
# write.xlsx(all.gwas.included, "~/data/combined_gwas_included.xlsx")


# Clean and filter the SNPs dataset using a single concise pipeline
Snps.GPAA <- Snps.GPAA %>%
  # 1. Separate the limits into two new columns with English names
  extract(
    col = X95..CI..TEXT., 
    into = c("lower_limit", "upper_limit"), 
    regex = "^\\[(-?[0-9.]+)-(-?[0-9.]+)\\].*", 
    remove = FALSE,  # Keeps original column intact
    convert = TRUE   # Converts extracted text to numbers automatically
  ) %>%
  # 2. Remove rows with missing values in effect size or confidence limits
  filter(!is.na(OR.or.BETA) & !is.na(lower_limit) & !is.na(upper_limit)) %>%
  # 3. Significance criterion: Exclude intervals crossing 1
  # Both values must be strictly above 1, or strictly below 1
  filter((lower_limit > 1 & upper_limit > 1) | (lower_limit < 1 & upper_limit < 1)) %>%
  # 4. Filter by strict SNP structure:
  # ^rs      : Starts with "rs"
  # \\d+     : Followed by one or more numbers
  # -        : Followed by a hyphen
  # [A-Z]$   : Ends EXACTLY with a single uppercase letter
  filter(str_detect(STRONGEST.SNP.RISK.ALLELE, "^rs\\d+-[A-Z]$"))

#******************************************************************************
#### 2.3 EXPORTING SIGNIFICANT SNPs ####
#******************************************************************************

# Extract unique SNPs and convert to a data frame
Filtered_SNPS <- as.data.frame(unique(Snps.GPAA$SNPS))

# Export the filtered SNPs (Optional)
# write_csv(Filtered_SNPS, file = "~/Desktop/glaucoma-gwas/data/Filtered_SNPS.csv")

# Load SNP frequencies (Uncomment path based on OS)
SNPs.FQs <- read.csv("~/data/SNPs_FQs.csv")
# SNPs.FQs <- read.csv("C:/Users/mloel/OneDrive/Desktop/glaucoma-gwas/data/SNPs_FQs.csv")

#******************************************************************************
#### 2.4 MERGING FREQUENCIES ####
#******************************************************************************

# Create a new column containing the specific allele
Snps.GPAA <- Snps.GPAA %>%
  mutate(Allele = str_extract(STRONGEST.SNP.RISK.ALLELE, "(?<=-)[A-Z]"))

# 1. Initialize empty population columns in the main dataset
Snps.GPAA$PLQ    <- NA
Snps.GPAA$CHG    <- NA
Snps.GPAA$CLM    <- NA
Snps.GPAA$ATQCES <- NA
Snps.GPAA$ATQPGC <- NA

# 2. Define the target population codes
poblaciones <- c("PLQ", "CHG", "CLM", "ATQCES", "ATQPGC")

# 3. Loop through the dataset to match and populate allele frequencies
for (i in 1:nrow(Snps.GPAA)) {
  
  snp_actual   <- Snps.GPAA$SNPS[i]
  alelo_actual <- Snps.GPAA$Allele[i] # The allele extracted in the previous step
  
  # Skip row to avoid errors if the allele is missing
  if (is.na(alelo_actual)) next
  
  for (pop in poblaciones) {
    if (is.na(Snps.GPAA[i, pop])) {
      
      # Find the exact match in the frequencies table
      match_idx <- which(SNPs.FQs$rs == snp_actual & SNPs.FQs$Population.Group == pop)
      
      if (length(match_idx) > 0) {
        idx <- match_idx[1]
        
        # Extract the frequency value for the specific allele (A, C, T, or G)
        valor_frecuencia <- SNPs.FQs[[alelo_actual]][idx]
        
        # Save the retrieved frequency into the main dataset
        Snps.GPAA[i, pop] <- valor_frecuencia
      }
    }
  }
}

# 4. Filter rows keeping only those with at least 2 valid population frequencies
Snps.GPAA.final <- Snps.GPAA[rowSums(!is.na(Snps.GPAA[, poblaciones])) >= 2, ]


#******************************************************************************
#### 3. DATA ANALYSIS ####
#******************************************************************************

#******************************************************************************
##### 3.1 Correlation Analysis #####
#******************************************************************************

# 1. Define the population groups we are going to compare
poblaciones <- c("PLQ", "CHG", "CLM", "ATQCES", "ATQPGC")

# 2. Get all possible combinations of 2 (pairs) without repetition
pares <- combn(poblaciones, 2, simplify = FALSE)

# 3. Create an empty data.frame to store the results neatly
resultados_cor <- data.frame(
  Grupo_1 = character(),
  Grupo_2 = character(),
  Shapiro_p_Grupo1 = numeric(),
  Shapiro_p_Grupo2 = numeric(),
  Metodo = character(),
  Coeficiente = numeric(),
  P_Valor_Cor = numeric(),
  stringsAsFactors = FALSE
)

# 4. Start the loop to analyze each pair of populations
for (par in pares) {
  pop1 <- par[1]
  pop2 <- par[2]
  
  # Extract only the two columns of interest from the main dataset
  datos_par <- Snps.GPAA[, c(pop1, pop2)]
  
  # Remove rows that have NA in either of the two populations
  # This is CRUCIAL because correlation requires paired observations
  datos_par <- na.omit(datos_par)
  
  # Convert to numeric vectors
  var1 <- as.numeric(datos_par[[pop1]])
  var2 <- as.numeric(datos_par[[pop2]])
  
  n_datos <- length(var1)
  
  # Validate that we have enough data (Shapiro requires at least 3 data points)
  if (n_datos >= 3) {
    
    # IMPORTANT: The shapiro.test() function in R has a maximum limit of 5000 data points.
    # If the table has more than 5000 SNPs, take a random sample of 5000 for the test.
    if (n_datos > 5000) {
      set.seed(123) # For reproducible sampling
      var1_shapiro <- sample(var1, 5000)
      var2_shapiro <- sample(var2, 5000)
    } else {
      var1_shapiro <- var1
      var2_shapiro <- var2
    }
    
    # Perform the Shapiro-Wilk normality test
    # (A p-value < 0.05 indicates that the data is NOT normally distributed)
    shapiro1 <- shapiro.test(var1_shapiro)$p.value
    shapiro2 <- shapiro.test(var2_shapiro)$p.value
    
    # 5. Decision logic: Pearson vs Spearman
    # If BOTH p-values are > 0.05, we assume normal distribution -> Pearson
    # If AT LEAST ONE is < 0.05, we cannot assume joint normality -> Spearman
    if (shapiro1 > 0.05 && shapiro2 > 0.05) {
      metodo_elegido <- "pearson"
    } else {
      metodo_elegido <- "spearman"
    }
    
    # 6. Execute the correlation test with the chosen method
    # (Use exact = FALSE to avoid warnings with ties in Spearman)
    test_correlacion <- cor.test(var1, var2, method = metodo_elegido, exact = FALSE)
    
    # 7. Save the results in a new row
    nueva_fila <- data.frame(
      Grupo_1 = pop1,
      Grupo_2 = pop2,
      Shapiro_p_Grupo1 = round(shapiro1, 4),
      Shapiro_p_Grupo2 = round(shapiro2, 4),
      Metodo = metodo_elegido,
      Coeficiente = round(test_correlacion$estimate, 4), # rho or r
      P_Valor_Cor = signif(test_correlacion$p.value, 4)
    )
    
    # Bind the new row to our main results table
    resultados_cor <- rbind(resultados_cor, nueva_fila)
    
  } else {
    # Warning message if a pair has almost no SNPs in common
    message(paste("Insufficient data to compare", pop1, "and", pop2))
  }
}

# 8. View the final table with all results
print(resultados_cor)

#******************************************************************************
######## 3.1.1 Plotting Correlation #########
#******************************************************************************

# 1. Prepare the data (using the results table from the previous step)
# Converting to factors ensures the matrix maintains the correct shape and order
datos_triangulo <- resultados_cor %>%
  mutate(
    Grupo_1 = factor(Grupo_1, levels = poblaciones),
    # Reverse the levels for the Y-axis so the triangle builds from top to bottom
    Grupo_2 = factor(Grupo_2, levels = rev(poblaciones)) 
  )

# 2. Generate the minimalist triangular plot
ggplot(datos_triangulo, aes(x = Grupo_1, y = Grupo_2, fill = Coeficiente)) +
  # Draw the correlation tiles with a clean white border
  geom_tile(color = "white", linewidth = 1) +
  # Apply a smooth, professional color gradient (Red -> White -> Blue)
  scale_fill_gradient2(
    high = "#ef233c", mid = "white", low = "#0582ca", 
    midpoint = 0, limit = c(-1, 1), space = "Lab", 
    name = "Correlation"
  ) +
  # Add ONLY the coefficient numbers, rounded to 2 decimals
  geom_text(aes(label = sprintf("%.2f", Coeficiente)), color = "black", size = 4.5) +
  # Apply an ultra-minimalist theme
  theme_minimal() + 
  # Remove all unnecessary visual clutter
  theme(
    axis.text.x = element_text(size = 11, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.x = element_blank(),       
    axis.title.y = element_blank(),       
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank(),   
    panel.background = element_blank(),   
    plot.background = element_blank(),
    legend.position = "right"             
  ) +
  # Force the tiles to be perfect squares
  coord_fixed()

#******************************************************************************
######## 3.1.2 PCA Analysis ######## 
#******************************************************************************

# 1. Extract the columns of interest BY NAME from the dataset
poblaciones_pca <- c("PLQ", "CHG", "CLM", "ATQCES", "ATQPGC")
data_pca <- Snps.GPAA[, poblaciones_pca]

# 2. Force everything to be numeric
data_pca[] <- lapply(data_pca, as.numeric)

# 3. Strict cleaning: Classic PCA does not support NAs (missing data). 
# Keep only the SNPs that have frequency data for ALL 5 populations.
data_pca_completo <- na.omit(data_pca)

# 4. Calculate variance and remove SNPs identical across all populations (variance = 0)
snp_variances <- apply(data_pca_completo, 1, var)
data_pca_clean <- data_pca_completo[snp_variances > 0, ]

# Print a message to see how much data survived the cleaning process
cat("From", nrow(data_pca), "original SNPs, there are", nrow(data_pca_clean), "complete and variant SNPs left for the PCA.\n")

# 5. Transpose the matrix (Populations as rows, SNPs as columns)
transposed_data <- t(data_pca_clean)

# 6. Run R's classic PCA
pca_result <- prcomp(transposed_data, center = TRUE, scale. = TRUE)

# 7. Extract the results for plotting
pca_scores <- as.data.frame(pca_result$x) 
pca_scores$Population <- rownames(pca_scores)

# Calculate the % of variance explained by each Principal Component
variance_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

# 8. Create the PCA plot
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 8, alpha = 0.9) +
  geom_text_repel(
    aes(label = Population),
    color = "black",
    box.padding   = 0.8,
    point.padding = 0.5,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20,
    force = 20,
    show.legend = FALSE
  ) +
  theme_bw(base_size = 14) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = sprintf("PC1 (%.1f%% variance)", variance_explained[1]),
    y = sprintf("PC2 (%.1f%% variance)", variance_explained[2])
  ) +
  # Define the limits from -10 to 10 for both axes
  coord_cartesian(xlim = c(-10, 9), ylim = c(-10, 9)) + 
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

#******************************************************************************
#### 3.2 DESCRIPTIVE ANALYSIS ####
#******************************************************************************

#******************************************************************************
########### 3.2.1 Top 15 Genes (Adapted for DISEASE.TRAIT) ########### 
#******************************************************************************

# 1. Find the top 15 most frequent genes (filtering out NAs and empty strings)
top_15_genes <- Snps.GPAA %>%
  filter(!is.na(MAPPED_GENE) & MAPPED_GENE != "") %>%  
  count(MAPPED_GENE, name = "Total") %>%               
  top_n(15, Total) %>%                                 
  pull(MAPPED_GENE)                                    

# 2. Filter the data to keep only those top 15 genes, 
# and count combinations of MAPPED_GENE and DISEASE.TRAIT
plot_data <- Snps.GPAA %>%
  filter(MAPPED_GENE %in% top_15_genes) %>%
  filter(!is.na(DISEASE.TRAIT) & DISEASE.TRAIT != "") %>% 
  count(MAPPED_GENE, DISEASE.TRAIT, name = "Frequency")  

# 3. Create the bar chart
ggplot(plot_data, aes(x = reorder(MAPPED_GENE, Frequency, sum), 
                      y = Frequency, 
                      fill = MAPPED_GENE)) +
  geom_bar(stat = "identity") +                
  coord_flip() +                               
  scale_fill_manual(values = c(
    "#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", 
    "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", 
    "#499894", "#86BCB6", "#E15759", "#FF9D9A", 
    "#79706E", "#BAB0AC", "#D37295", "#FABFD2"
  )) +
  labs(
    x = "Genes",
    y = "Number of Associations"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.position = "none"
  )

#******************************************************************************
###### 3.2.2 Jitter Plot (Adapted for GWAS: Ancestry Disparities) ######
#******************************************************************************

# 1. Calculate the mean for admixed/European populations
Snps.GPAA <- Snps.GPAA %>%
  mutate(
    mean.admixed = rowMeans(select(., CLM, ATQCES, ATQPGC), na.rm = TRUE)
  )

# 2. Filter and classify the data (Strict limits at <0.40 and >0.60)
plot_data_traits <- Snps.GPAA %>%
  mutate(
    Trait_Group = case_when(
      # High frequency in Palenque (>=60%) and low in Admixed (<=40%)
      PLQ >= 0.60 & mean.admixed <= 0.40 ~ "PLQ Freq >= 60%",
      # High frequency in Admixed (>=60%) and low in Palenque (<=40%)
      mean.admixed >= 0.60 & PLQ <= 0.40 ~ "Admixed Freq >= 60%",
      # All other data points
      TRUE ~ "Other"
    )
  ) %>%
  # Sort to ensure colored points are drawn on top of the grey ones
  arrange(desc(Trait_Group == "Other"))

# 3. Restore the original colors
colors_traits <- c(
  "PLQ Freq >= 60%"     = "#3a86ff", # Vibrant blue
  "Admixed Freq >= 60%" = "#ff006e", # Neon pink
  "Other"               = "#d9d9d9"  # Original light grey
)

# 4. Generate the plot
ggplot(plot_data_traits, aes(x = PLQ, y = mean.admixed, color = Trait_Group)) +
  # Exact limit lines at 0.4 and 0.6
  geom_hline(yintercept = c(0.4, 0.6), linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = c(0.4, 0.6), linetype = "dashed", color = "gray50", linewidth = 0.5) +
  # Point size increased for visibility
  geom_jitter(size = 2.5, alpha = 0.7, width = 0.02, height = 0.02) +
  scale_color_manual(values = colors_traits) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  labs(
    x = "African (San Basilio de Palenque)",
    y = "European",
    color = "Highlighted Traits"
  ) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.text = element_text(color = "black"),
    legend.position = "none" 
  )

#******************************************************************************
######## 3.3 Number of Associations per Chromosome ####### 
#******************************************************************************

plot_data_chr <- Snps.GPAA %>%
  # 1. Remove explicit NAs, empty strings, and the word "NA" if read as text
  filter(!is.na(CHR_ID) & CHR_ID != "" & CHR_ID != "NA") %>%
  # Extract the "p" or "q" arm from the REGION column
  mutate(Arm = str_extract(REGION, "[pq]")) %>%
  # Remove rows where the Arm couldn't be determined
  filter(!is.na(Arm)) %>%
  # Count the number of associations per Chromosome and per Arm
  count(CHR_ID, Arm, name = "Associations") %>%
  # Ensure no NAs remain in the final dataframe
  drop_na()

# 2. Fix the Chromosome Order (Biological order: 1 to 22, then X and Y)
chr_order <- c(as.character(1:22), "X", "Y")

# Use rev() so Chromosome 1 appears at the TOP when coordinates are flipped
plot_data_chr$CHR_ID <- factor(plot_data_chr$CHR_ID, levels = rev(chr_order))

# Discard any NAs introduced if there were chromosomes outside chr_order (e.g., mitochondrial)
plot_data_chr <- plot_data_chr %>% filter(!is.na(CHR_ID))

# 3. Create the Figure
ggplot(plot_data_chr, aes(y = CHR_ID, x = Associations, fill = Arm)) +
  # Draw the bars (mapping Y = CHR_ID draws them horizontally automatically)
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("p" = "#ff006e", "q" = "#3a86ff")) +
  # Add 'n = 4' to pretty() to force fewer numbers and space them out horizontally
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, max(x, na.rm = TRUE)), n = 4)))) +
  labs(
    y = "Chromosome",
    x = "N. of Associations",
    fill = "Arm"
  ) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = "gray80", linetype = "dashed", linewidth = 0.6),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(), # No minor lines
    axis.line = element_line(color = "black", linewidth = 0.8), # Thick black axis lines
    axis.title = element_text(size = 14, face = "bold", color = "black"), # Bold axis titles
    axis.text = element_text(size = 12, color = "black"), # General axis text formatting
    axis.text.x = element_text(margin = margin(t = 12)),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

#******************************************************************************
##### 3.4 Ancestry Analysis (Adapted for GWAS Catalog Text Mining) #####
#******************************************************************************

final_table <- data.frame(text = unique(Snps.GPAA$INITIAL.SAMPLE.SIZE), stringsAsFactors = FALSE) %>%
  # 1. Text cleaning and separation (removes commas inside numbers)
  mutate(clean_text = str_replace_all(text, "(?<=\\d),(?=\\d)", "")) %>%
  separate_rows(clean_text, sep = ",\\s*") %>%
  # 2. Extraction and classification
  mutate(
    Count = as.numeric(str_extract(clean_text, "\\d+")),
    Type = case_when(
      str_detect(clean_text, "(?i)cases") ~ "cases",
      str_detect(clean_text, "(?i)controls") ~ "controls",
      str_detect(clean_text, "(?i)individuals") ~ "individuals",
      TRUE ~ "other"
    ),
    Clean_Ancestry = case_when(
      str_detect(clean_text, "(?i)European|Erasmus|Orcadian") ~ "European",
      str_detect(clean_text, "(?i)Asian|Japanese|Chinese|Indian") ~ "Asian",
      str_detect(clean_text, "(?i)African|Sub-Saharan") ~ "African",
      str_detect(clean_text, "(?i)Hispanic|Latino|Latin American") ~ "Latino",
      TRUE ~ "Other" # Catches mixed ancestry and anything unrecognized
    )
  ) %>%
  # Filter out missing numbers and keep only cases/controls
  filter(!is.na(Count), Type %in% c("cases", "controls")) %>%
  # 3. Summarization and pivoting
  group_by(Clean_Ancestry, Type) %>%
  summarise(Total = sum(Count, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Type, values_from = Total, values_fill = list(Total = 0)) %>%
  # 4. Calculation of totals and percentages per row
  mutate(
    Total_Individuals = cases + controls,
    Percent_Cases = round((cases / sum(cases)) * 100, 2),
    Percent_Controls = round((controls / sum(controls)) * 100, 2),
    Percent_Total = round((Total_Individuals / sum(Total_Individuals)) * 100, 2)
  ) %>%
  # 5. Add the GRAND TOTAL row dynamically
  {
    bind_rows(., summarise(., 
                           Clean_Ancestry = "GRAND TOTAL",
                           cases = sum(cases),
                           controls = sum(controls),
                           Total_Individuals = sum(Total_Individuals),
                           Percent_Cases = 100,
                           Percent_Controls = 100,
                           Percent_Total = 100
    ))
  }

#******************************************************************************
##### 3.5 Chord Plot #####
#******************************************************************************

ancestry_totals <- final_table %>%
  # Remove percentage columns
  select(-Percent_Cases, -Percent_Controls, -Percent_Total) %>%
  # Filter out the GRAND TOTAL row for the plot
  filter(Clean_Ancestry != "GRAND TOTAL") %>%
  # Rename columns to match links dataframe structure
  rename(Ancestry_Unified = Clean_Ancestry) %>%
  rename(Total_Cases = cases, Total_Controls = controls) %>%
  rename(Grand_Total = Total_Individuals)

# Order the ancestry groups by the total number of individuals (descending)
ancestry_totals <- ancestry_totals %>% arrange(desc(Grand_Total))

# Create a links dataframe for the chord diagram
links <- data.frame(
  from = rep(ancestry_totals$Ancestry_Unified, 2),  # Each ancestry appears twice (cases and controls)
  to = c(rep("casos", nrow(ancestry_totals)), rep("controles", nrow(ancestry_totals))),
  value = c(ancestry_totals$Total_Cases, ancestry_totals$Total_Controls)  # Corresponding values
)

# Remove connections with a value of zero to avoid unnecessary lines in the diagram
links <- links[links$value > 0, ]

# Manually define the order of the ancestry groups
order_manual <- c("European", "Asian", "African", "Latino", "Other")

# Define the fixed order for the 'cases' and 'controls' sectors
order_casos_controles <- c("casos", "controles")

# Apply the desired order to the factors in the links dataframe
links$from <- factor(links$from, levels = order_manual)
links$to <- factor(links$to, levels = order_casos_controles)

# Define custom colors for each ancestry group and the 'cases'/'controls'
grid.col <- c(
  "European"  = "#6CB4EE", 
  "Asian"     = "#1B4D3E", 
  "African"   = "#FFBF00",
  "Latino"    = "#6C3082",
  "Other"     = "#3AB09E",
  "casos"     = "black", 
  "controles" = "darkgray"
)

# Clear any previous circlize plots
circos.clear()

# Define reversed order to display the ancestry sectors counterclockwise
orden <- rev(c("European", "Asian", "African", "Latino", "Other"))

# Update the 'from' column with the reversed order
links$from <- factor(links$from, levels = orden)

# Sort the links dataframe according to the updated ancestry order
links <- links[order(links$from), ]

# Set custom transparency for each ancestry group
# A value of 0 = completely opaque, 1 = completely transparent
transparency_values <- c(
  "European" = 0.2, 
  "Asian"    = 0.2, 
  "African"  = 0.2, 
  "Latino"   = 0.2, 
  "Other"    = 0.2
)

# Generate the chord diagram
chordDiagram(
  links,
  grid.col = grid.col, # Use the custom color mapping
  transparency = transparency_values[as.character(links$from)], # Apply specific transparency
  annotationTrack = "grid", # Add grid lines for the sectors
  preAllocateTracks = list(track.height = 0.1) # Pre-allocate track height
)