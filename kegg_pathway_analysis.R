setwd("C:/B.TECH/COVID19_VS_HEALTHY")

library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)

# Load KEGG file
data <- read_excel("COVID19_Upregulated_KEGG_DAVID.xlsx")

# Clean columns
colnames(data) <- trimws(colnames(data))

# Fix Benjamini
data$Benjamini <- as.character(data$Benjamini)
data$Benjamini <- gsub(",", "", data$Benjamini)
data$Benjamini <- as.numeric(data$Benjamini)

# Filter significant
kegg_sig <- data %>%
  filter(Benjamini < 0.05)

# Top pathways
top_kegg <- kegg_sig %>%
  arrange(Benjamini) %>%
  head(8)

# Create logP
top_kegg$logP <- -log10(top_kegg$Benjamini)

# Clean names
top_kegg$Term <- str_wrap(top_kegg$Term, width = 40)

# Dotplot
ggplot(top_kegg, aes(x = logP, y = reorder(Term, logP))) +
  geom_point(aes(size = Count, color = logP)) +
  scale_color_gradient(low = "blue", high = "red") +
  xlab("-log10 (Adjusted P-value)") +
  ylab("KEGG Pathways") +
  ggtitle("KEGG Pathway Enrichment (Upregulated Genes)") +
  theme_minimal()

ggsave("KEGG_dotplot_upregulated.png", width = 12, height = 8)




# ================================
# LOAD KEGG FILE (DOWNREGULATED)
# ================================
data <- read_excel("COVID19_DOWNREGULATED_KEGG_DAVID.xlsx")

# ================================
# CLEAN COLUMN NAMES
# ================================
colnames(data) <- trimws(colnames(data))

# ================================
# FIX BENJAMINI COLUMN
# ================================
data$Benjamini <- as.character(data$Benjamini)
data$Benjamini <- gsub(",", "", data$Benjamini)
data$Benjamini <- as.numeric(data$Benjamini)

# ================================
# FILTER SIGNIFICANT KEGG PATHWAYS
# ================================
kegg_sig <- data %>%
  filter(Benjamini < 0.05)

# ???? CHECK (VERY IMPORTANT)
print(nrow(kegg_sig))

# ================================
# SELECT TOP PATHWAYS
# ================================
top_kegg <- kegg_sig %>%
  arrange(Benjamini) %>%
  head(8)

# ================================
# CREATE -log10 VALUES
# ================================
top_kegg$logP <- -log10(top_kegg$Benjamini)

# ================================
# CLEAN LONG PATHWAY NAMES
# ================================
top_kegg$Term <- str_wrap(top_kegg$Term, width = 40)

# ================================
# DOTPLOT
# ================================
ggplot(top_kegg, aes(x = logP, y = reorder(Term, logP))) +
  geom_point(aes(size = Count, color = logP)) +
  scale_color_gradient(low = "blue", high = "red") +
  xlab("-log10 (Adjusted P-value)") +
  ylab("KEGG Pathways") +
  ggtitle("KEGG Pathway Enrichment (Downregulated Genes)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

ggsave("KEGG_dotplot_downregulated.png", width = 12, height = 8)