#' Encode equal-length DNA strings into positional factor features
#'
#' @description
#' Given a set of DNA sequences (e.g., 5-mers) composed of A/T/C/G with the same length,
#' this function splits them by position and encodes each base as a factor
#' (levels = A, T, C, G). The resulting columns are named nt_pos1 to nt_posN.
#'
#' @param dna_strings A character vector of DNA strings composed only of A, T, C, and G,
#'   all of the same length.
#'
#' @return A data frame where:
#'   - Columns are named nt_pos1, nt_pos2, …, nt_posN
#'   - Each column is a factor with levels c("A","T","C","G")
#'
#' @details
#' The function does not check the validity of the bases automatically;
#' please ensure all input strings use uppercase A/T/C/G and have identical lengths.
#'
#' @examples
#' # Minimal example for encoding 5-mers
#' dna_encoding(c("ATCGA", "GGGTT"))
#'
#' @seealso [prediction_multiple()], [prediction_single()]
#' @export
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A","T","C","G"))
  seq_df
}

#' Predict m6A status and probability for multiple candidate sites
#'
#' @description
#' Uses a trained randomForest classification model to predict the m6A probability
#' and class label for multiple candidate sites provided in feature_df.
#' If the table does not already contain nt_pos1–nt_pos5, these columns are
#' automatically generated using [dna_encoding()] from the DNA_5mer column.
#'
#' @param ml_fit A fitted randomForest model object (with the positive class labeled "Positive").
#' @param feature_df A data frame containing at least the following columns:
#'   gc_content, RNA_type, RNA_region, exon_length,
#'   distance_to_junction, evolutionary_conservation, and DNA_5mer.
#'
#'   Optionally, may include pre-encoded nt_pos1–nt_pos5 factor columns.
#' @param positive_threshold Numeric. The probability threshold above which
#'   predictions are labeled "Positive". Default = 0.5.
#'
#' @return A data frame identical to feature_df but with two additional columns:
#'   - predicted_m6A_prob: numeric probability of "Positive" prediction
#'   - predicted_m6A_status: predicted class label ("Positive" or "Negative")
#'
#' @details
#' The function:
#' 1. Validates required feature columns
#' 2. Converts categorical features into fixed-level factors
#' 3. Converts numeric columns properly
#' 4. Encodes 5-mer nucleotides if missing
#' 5. Predicts probabilities using predict(..., type = "prob")
#'
#' @import randomForest
#' @importFrom stats predict
#' @export
#'
#' @examples
#' # ---- Working example that uses package-bundled data (inst/extdata) ----
#' # Load fitted RF model (.rds) and example features (.csv) from the installed package
#' rf_path  <- system.file("extdata", "rf_fit.rds", package = "m6APrediction")
#' csv_path <- system.file("extdata", "m6A_input_example.csv", package = "m6APrediction")
#' ml_fit   <- readRDS(rf_path)
#' feature_df <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
#'
#' # Run batch prediction; should add prob/label columns without errors
#' preds <- prediction_multiple(ml_fit, feature_df, positive_threshold = 0.5)
#' head(preds)
#'
#' @seealso [dna_encoding()], [prediction_single()], [randomForest::randomForest()]
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content","RNA_type","RNA_region","exon_length",
                  "distance_to_junction","evolutionary_conservation","DNA_5mer")
                %in% colnames(feature_df)))

  feature_df$RNA_type   <- factor(feature_df$RNA_type,
                                  levels = c("mRNA","lincRNA","lncRNA","pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS","intron","3'UTR","5'UTR"))

  num_cols <- c("gc_content","exon_length","distance_to_junction","evolutionary_conservation")
  feature_df[num_cols] <- lapply(feature_df[num_cols], function(x) as.numeric(as.character(x)))

  nt_cols <- paste0("nt_pos", 1:5)
  if (!all(nt_cols %in% names(feature_df))) {
    dna_df <- dna_encoding(feature_df$DNA_5mer)
    feature_df <- cbind(feature_df, dna_df)
  } else {
    for (k in nt_cols) feature_df[[k]] <- factor(feature_df[[k]], levels = c("A","T","C","G"))
  }

  prob <- predict(ml_fit, newdata = feature_df, type = "prob")[, "Positive"]

  feature_df$predicted_m6A_prob   <- as.numeric(prob)
  feature_df$predicted_m6A_status <- ifelse(feature_df$predicted_m6A_prob > positive_threshold,
                                            "Positive","Negative")
  feature_df
}

#' Predict m6A status for a single sample (convenience wrapper)
#'
#' @description
#' A simple wrapper function for predicting a single record’s m6A probability and label.
#' Internally constructs a one-row feature data frame and calls [prediction_multiple()].
#'
#' @param ml_fit A fitted randomForest model object.
#' @param gc_content Numeric GC content (fraction or between 0–1).
#' @param RNA_type Character/factor. One of "mRNA", "lincRNA", "lncRNA", "pseudogene".
#' @param RNA_region Character/factor. One of "CDS", "intron", "3'UTR", "5'UTR".
#' @param exon_length Numeric exon length.
#' @param distance_to_junction Numeric distance to the splice junction.
#' @param evolutionary_conservation Numeric conservation score.
#' @param DNA_5mer Character string for the 5-mer DNA sequence (e.g. "ATCGA").
#' @param positive_threshold Numeric threshold for labeling "Positive" (default = 0.5).
#'
#' @return A named vector with:
#' - predicted_m6A_prob: predicted probability of "Positive"
#' - predicted_m6A_status: predicted label ("Positive"/"Negative")
#'
#' @export
#'
#' @examples
#' # ---- Working example using the same bundled data as above ----
#' rf_path  <- system.file("extdata", "rf_fit.rds", package = "m6APrediction")
#' csv_path <- system.file("extdata", "m6A_input_example.csv", package = "m6APrediction")
#' ml_fit   <- readRDS(rf_path)
#' feature_df <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
#'
#' # Take the first row to form a single-record prediction
#' res <- prediction_single(
#'   ml_fit,
#'   gc_content = feature_df$gc_content[1],
#'   RNA_type   = feature_df$RNA_type[1],
#'   RNA_region = feature_df$RNA_region[1],
#'   exon_length = feature_df$exon_length[1],
#'   distance_to_junction = feature_df$distance_to_junction[1],
#'   evolutionary_conservation = feature_df$evolutionary_conservation[1],
#'   DNA_5mer  = feature_df$DNA_5mer[1],
#'   positive_threshold = 0.5
#' )
#' res
#'
#' @seealso [prediction_multiple()], [dna_encoding()]
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction, evolutionary_conservation,
                              DNA_5mer, positive_threshold = 0.5){
  feature_df <- data.frame(
    gc_content = as.numeric(gc_content),
    RNA_type = factor(RNA_type, levels = c("mRNA","lincRNA","lncRNA","pseudogene")),
    RNA_region = factor(RNA_region, levels = c("CDS","intron","3'UTR","5'UTR")),
    exon_length = as.numeric(exon_length),
    distance_to_junction = as.numeric(distance_to_junction),
    evolutionary_conservation = as.numeric(evolutionary_conservation),
    DNA_5mer = as.character(DNA_5mer),
    stringsAsFactors = FALSE
  )
  pred_df <- prediction_multiple(ml_fit, feature_df, positive_threshold)
  c(predicted_m6A_prob   = as.numeric(pred_df$predicted_m6A_prob[1]),
    predicted_m6A_status = as.character(pred_df$predicted_m6A_status[1]))
}
