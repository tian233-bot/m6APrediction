# 🧬 m6APrediction

### Predict m6A RNA Methylation Sites Using Random Forest Models

---

## 📘 Overview

**m6APrediction** is an R package developed to predict **N6-methyladenosine (m6A)** RNA methylation sites in RNA sequences.  
It provides user-friendly functions for:

- 🧩 **Encoding** DNA 5-mer sequences into position-based categorical features  
- 🔍 **Predicting** m6A modification probabilities using a trained **Random Forest** model  
- ⚡ Supporting both **single** and **multiple** sample predictions  

This package was created for the **BIO215 coursework** to demonstrate the application of machine learning in RNA modification analysis.

---

## ⚙️ Installation

You can install the package directly from GitHub using **devtools** or **remotes**:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install m6APrediction from GitHub
devtools::install_github("tian233-bot/m6APrediction")

# Or alternatively, using remotes
# remotes::install_github("tian233-bot/m6APrediction")
---

## 🧬 Example Usage

### 💡 Example 1: Encode DNA 5-mer Sequences
```{r}
library(m6APrediction)

# Example: Encode 5-mer sequences into positional factor features
dna_encoding(c("ATCGA", "GGGTT"))
```

###🔁 Example 2: Predict Multiple m6A Sites
```{r}
library(m6APrediction)
library(randomForest)

# Example feature table (minimal viable columns)
feature_df <- data.frame(
  gc_content = c(0.45, 0.62),
  RNA_type   = c("mRNA", "lncRNA"),
  RNA_region = c("CDS", "3'UTR"),
  exon_length = c(1500, 800),
  distance_to_junction = c(120, 50),
  evolutionary_conservation = c(0.8, 0.3),
  DNA_5mer = c("ATCGA", "GGGTT"),
  stringsAsFactors = FALSE
)

# A small toy random forest model for demonstration
set.seed(123)
ml_fit <- randomForest::randomForest(
  x = dna_encoding(feature_df$DNA_5mer),
  y = factor(c("Positive", "Negative"))
)

# Run batch prediction
preds <- prediction_multiple(ml_fit, feature_df)
head(preds)
```

###🎯 Example 3: Predict a Single m6A Site
```{r}
# Predict m6A status for one example sequence
prediction_single(
  ml_fit,
  gc_content = 0.55,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 1200,
  distance_to_junction = 80,
  evolutionary_conservation = 0.7,
  DNA_5mer = "ATCGA"
)
```

##📈 Model Performance Visualization

###To showcase the model’s predictive power, include the ROC and PRC curve images from my Practical 4 results.
These figures should be placed under man/figures/ within your package directory.

![](figures/ROC_curve.png)
![](figures/PRC_curve.png)

