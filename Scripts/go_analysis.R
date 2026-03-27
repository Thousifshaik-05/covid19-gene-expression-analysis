setwd("C:/B.TECH/COVID19_VS_HEALTHY")

# Install packages (run once)
install.packages("readxl")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("stringr")

# Load libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)

# GO VISUALIZATION - UPREGULATED GENES

# ================================
# LOAD DAVID FILE
# ================================
data <- read_excel("COVID19_Upregulated_GO_DAVID.xlsx")

# ================================
# Clean columns
colnames(data) <- trimws(colnames(data))

# Fix Benjamini first
data$Benjamini <- as.character(data$Benjamini)
data$Benjamini <- gsub(",", "", data$Benjamini)
data$Benjamini <- as.numeric(data$Benjamini)

# Filter BP terms (flexible)
bp_data <- data %>%
  filter(grepl("GOTERM_BP", Category))

# Filter significant
bp_sig <- bp_data %>%
  filter(Benjamini < 0.05)

# Check rows
nrow(bp_sig)

# Take top terms
top_bp <- bp_sig %>%
  arrange(Benjamini) %>%
  head(8)
# ================================
# CREATE -log10 VALUES
# ================================
top_bp$logP <- -log10(top_bp$Benjamini)

# ================================
# CLEAN LONG GO TERMS
# ================================
top_bp$Term <- str_wrap(top_bp$Term, width = 40)

# ================================
# BARPLOT
# ================================
ggplot(top_bp, aes(x = reorder(Term, logP), y = logP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  xlab("GO Terms") +
  ylab("-log10 (Adjusted P-value)") +
  ggtitle("Top GO Biological Processes (Upregulated Genes)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

ggsave("GO_barplot_upregulated.png", width = 10, height = 7)

# ================================
# DOTPLOT
# ================================
ggplot(top_bp, aes(x = logP, y = reorder(Term, logP))) +
  geom_point(aes(size = Count, color = logP)) +
  scale_color_gradient(low = "blue", high = "red") +
  xlab("-log10 (Adjusted P-value)") +
  ylab("GO Terms") +
  ggtitle("GO Enrichment Dotplot (Upregulated Genes)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

ggsave("GO_dotplot_upregulated.png", width = 12, height = 8)





# ================================
# LOAD DAVID FILE (DOWNREGULATED)
# ================================
data <- read_excel("COVID19_Downregulated_GO_DAVID.xlsx")

# ================================
# CLEAN COLUMN NAMES
# ================================
colnames(data) <- trimws(colnames(data))

# ================================
# FIX BENJAMINI COLUMN (IMPORTANT)
# ================================
data$Benjamini <- as.character(data$Benjamini)
data$Benjamini <- gsub(",", "", data$Benjamini)
data$Benjamini <- as.numeric(data$Benjamini)

# ================================
# FILTER GO BIOLOGICAL PROCESS
# ================================
bp_data <- data %>%
  filter(grepl("GOTERM_BP", Category))

# ================================
# FILTER SIGNIFICANT TERMS
# ================================
bp_sig <- bp_data %>%
  filter(Benjamini < 0.05)

# ???? CHECK (IMPORTANT)
print(nrow(bp_sig))

# ================================
# SELECT TOP TERMS
# ================================
top_bp <- bp_sig %>%
  arrange(Benjamini) %>%
  head(8)

# ================================
# CREATE -log10 VALUES
# ================================
top_bp$logP <- -log10(top_bp$Benjamini)

# ================================
# CLEAN LONG TERMS
# ================================
top_bp$Term <- str_wrap(top_bp$Term, width = 40)

# ================================
# BARPLOT
# ================================
ggplot(top_bp, aes(x = reorder(Term, logP), y = logP)) +
  geom_bar(stat = "identity", fill = "darkred") +
  coord_flip() +
  xlab("GO Terms") +
  ylab("-log10 (Adjusted P-value)") +
  ggtitle("Top GO Biological Processes (Downregulated Genes)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

ggsave("GO_barplot_downregulated.png", width = 10, height = 7)

# ================================
# DOTPLOT
# ================================
ggplot(top_bp, aes(x = logP, y = reorder(Term, logP))) +
  geom_point(aes(size = Count, color = logP)) +
  scale_color_gradient(low = "blue", high = "red") +
  xlab("-log10 (Adjusted P-value)") +
  ylab("GO Terms") +
  ggtitle("GO Enrichment Dotplot (Downregulated Genes)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

ggsave("GO_dotplot_downregulated.png", width = 12, height = 8)




# ================================
# LOAD FILES
# ================================
up <- read_excel("COVID19_Upregulated_GO_DAVID.xlsx")
down <- read_excel("COVID19_Downregulated_GO_DAVID.xlsx")

# ================================
# CLEAN COLUMN NAMES
# ================================
# Function to clean columns consistently
clean_data <- function(df) {
  colnames(df) <- trimws(colnames(df))
  
  # Fix Benjamini
  df$Benjamini <- as.character(df$Benjamini)
  df$Benjamini <- gsub(",", "", df$Benjamini)
  df$Benjamini <- as.numeric(df$Benjamini)
  
  # Fix FDR if present
  if("FDR" %in% colnames(df)) {
    df$FDR <- as.character(df$FDR)
    df$FDR <- gsub(",", "", df$FDR)
    df$FDR <- as.numeric(df$FDR)
  }
  
  return(df)
}

# Apply to both
up <- clean_data(up)
down <- clean_data(down)

# ================================
# FIX BENJAMINI
# ================================
fix_benjamini <- function(df) {
  df$Benjamini <- as.character(df$Benjamini)
  df$Benjamini <- gsub(",", "", df$Benjamini)
  df$Benjamini <- as.numeric(df$Benjamini)
  return(df)
}

up <- fix_benjamini(up)
down <- fix_benjamini(down)

# ================================
# FILTER BP + SIGNIFICANT
# ================================
up_bp <- up %>%
  filter(grepl("GOTERM_BP", Category)) %>%
  filter(Benjamini < 0.05)

down_bp <- down %>%
  filter(grepl("GOTERM_BP", Category)) %>%
  filter(Benjamini < 0.05)

# ================================
# SELECT TOP TERMS
# ================================
up_top <- up_bp %>%
  arrange(Benjamini) %>%
  head(6) %>%
  mutate(Regulation = "Upregulated")

down_top <- down_bp %>%
  arrange(Benjamini) %>%
  head(6) %>%
  mutate(Regulation = "Downregulated")

# ================================
# COMBINE
# ================================
combined <- bind_rows(up_top, down_top)

# ================================
# CREATE logP
# ================================
combined$logP <- -log10(combined$Benjamini)

# ================================
# CLEAN TERMS
# ================================
combined$Term <- str_wrap(combined$Term, width = 40)

# ================================
# COMBINED DOTPLOT
# ================================
ggplot(combined, aes(x = logP, y = reorder(Term, logP))) +
  geom_point(aes(size = Count, color = Regulation)) +
  xlab("-log10 (Adjusted P-value)") +
  ylab("GO Terms") +
  ggtitle("GO Enrichment: Up vs Downregulated Genes") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

ggsave("GO_combined_dotplot.png", width = 12, height = 8)