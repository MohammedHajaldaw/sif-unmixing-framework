# Function to calculate class-average SIF emission and select fraction thresholds
calculate_class_avg_sif <- function(sif_cont_pred, frac,
                                    candidate_thresholds = seq(.01, .1, .005)) {
  
  # Overview:
  #   Converts predicted SIF contribution to class-average SIF emission,
  #   evaluates candidate fraction thresholds using modified z-scores,
  #   and selects the threshold with the lowest outlier rate.
  
  # Requires:
  #   - sif_cont_pred (SpatRaster): Predicted SIF contribution per vegetation class.
  #   - frac (SpatRaster): Raster stack containing vegetation class fractions.
  #   - candidate_thresholds (numeric): Candidate minimum fraction thresholds.
  
  # Effects:
  #   - Calculates class-average SIF emission from SIF contribution and fraction.
  #   - Evaluates outlier rates for candidate fraction thresholds.
  #   - Produces diagnostic threshold plots.
  
  # Returns:
  #   - A SpatRaster containing class-average SIF emission after fraction filtering.
  
  if (nlyr(sif_cont_pred) != nlyr(frac))
    stop("sif_cont_pred and frac must have the same number of layers.")
  
  n <- nlyr(frac); selected <- numeric(n)
  selected_outlier <- numeric(n); selected_retained <- numeric(n)
  
  # Set plotting layout for the three vegetation classes.
  op <- par(mfrow = c(3, 1), mar = c(4.5, 5.2, 2.5, 1),
            mgp = c(3.7, .8, 0), las = 1, bty = "l")
  on.exit(par(op))
  
  cols <- c(Crops = "#984EA3", Mixed_Vegetation = "#a6d854", Trees = "#1C7255")
  
  # Evaluate fraction thresholds separately for each vegetation class.
  for (i in seq_len(n)) {
    
    # Extract class name, fraction values, and predicted SIF contribution.
    class_name <- sub(" fractions$", "", names(frac)[i], ignore.case = TRUE)
    f <- values(frac[[i]], mat = FALSE)
    sif_cont <- values(sif_cont_pred[[i]], mat = FALSE)
    
    # Retain finite pixels with a positive class fraction and calculate
    # class-average SIF emission as SIF contribution divided by fraction.
    ok <- is.finite(f) & is.finite(sif_cont) & f > 0
    f <- f[ok]
    sif_avg <- sif_cont[ok] / f
    total_valid <- length(f)
    
    # Store outlier rate and number of retained pixels for each threshold.
    z <- data.frame(
      threshold = candidate_thresholds,
      outlier_percent = NA_real_,
      retained_n = NA_integer_,
      retained_percent = NA_real_
    )
    
    # Evaluate each candidate minimum fraction threshold.
    for (k in seq_along(candidate_thresholds)) {
      frac_val <- candidate_thresholds[k]
      x <- sif_avg[f >= frac_val]
      
      z$retained_n[k] <- length(x)
      z$retained_percent[k] <- 100 * length(x) / total_valid
      
      if (length(x) < 2) next
      
      # Calculate the median absolute deviation and modified z-scores.
      med <- median(x, na.rm = TRUE)
      mad_value <- mad(x, med, constant = 1, na.rm = TRUE)
      
      if (!is.finite(mad_value) || mad_value == 0) {
        z$outlier_percent[k] <- 0
      } else {
        modified_z <- .6745 * (x - med) / mad_value
        z$outlier_percent[k] <- mean(abs(modified_z) > 3.5, na.rm = TRUE) * 100
      }
    }
    
    # Calculate the reduction in outlier rate between consecutive thresholds.
    z$reduction <- NA_real_
    if (nrow(z) > 1)
      z$reduction[-1] <- z$outlier_percent[-nrow(z)] - z$outlier_percent[-1]
    
    # Select the threshold with the minimum outlier rate.
    min_outlier <- min(z$outlier_percent, na.rm = TRUE)
    min_indices <- which(abs(z$outlier_percent - min_outlier) < 1e-12)
    selected_index <- min_indices[1]
    
    # If several thresholds have the same minimum outlier rate,
    # select the first threshold for which an earlier outlier rate was higher.
    if (length(min_indices) > 1) {
      for (j in min_indices) {
        if (j == 1) next
        if (any(z$outlier_percent[1:(j - 1)] > min_outlier + 1e-12, na.rm = TRUE)) {
          selected_index <- j
          break
        }
      }
    }
    
    selected[i] <- z$threshold[selected_index]
    selected_outlier[i] <- z$outlier_percent[selected_index]
    selected_retained[i] <- z$retained_percent[selected_index]
    
    # Plot outlier rate as a function of the minimum class fraction.
    x <- z$threshold
    y <- z$outlier_percent
    col <- unname(cols[gsub(" ", "_", class_name)])
    if (is.na(col)) col <- "black"
    
    plot(x, y, type = "o", pch = 16, col = col, lwd = 1.5, cex = 1.3,
         xlab = if (i == n) "Minimum class fraction" else "",
         ylab = "Outlier rate (%)", main = class_name, xaxt = "n",
         cex.main = 1.5, cex.lab = 1.35, cex.axis = 1.2)
    axis(1, x)
    abline(h = pretty(y), lty = 3, col = "grey85")
    
    # Mark all thresholds with the minimum outlier rate.
    best <- which(abs(z$outlier_percent - min_outlier) < 1e-12)
    points(x[best], y[best], pch = 21, bg = col, col = "white",
           lwd = 1.5, cex = 1.7)
    
    # Highlight the selected threshold.
    j <- selected_index
    points(x[j], y[j], pch = 21, bg = col, col = "black",
           lwd = 2, cex = 2.4)
    
    # Position the threshold annotation so it does not overlap the plot edge.
    usr <- par("usr")
    dy <- diff(usr[3:4])
    gap <- .3 * dy
    left <- x[j] <= median(x)
    text_x <- x[j] + if (left) .006 else -.006
    text_pos <- if (left) 4 else 2
    
    if (y[j] < usr[3] + .35 * dy) {
      text_y <- y[j] + gap
      arrow_y <- text_y - .04 * dy
    } else {
      text_y <- y[j] - gap
      arrow_y <- text_y + .04 * dy
    }
    
    text(text_x, text_y, paste0("Best threshold = ", x[j]),
         pos = text_pos, font = 2, cex = .95)
    arrows(text_x, arrow_y, x[j], y[j], length = .08, lwd = 1.5)
  }
  
  # Convert predicted SIF contribution to class-average SIF emission
  sif_avg_pred <- sif_cont_pred
  
  for (i in seq_len(n)) {
    f <- frac[[i]]
    sif_avg <- sif_cont_pred[[i]] / f
    sif_avg[f < selected[i]] <- 0 # Set to zero for visualization; treated as NA in analysis.
    sif_avg_pred[[i]] <- sif_avg
  }
  
  # Crop and mask the output to the extent of the predicted SIF contribution.
  sif_avg_pred <- crop(sif_avg_pred, sif_cont_pred, touches = FALSE, mask = TRUE)
  names(sif_avg_pred) <- c("Crops", "Mixed_Vegetation", "Trees")
  
  return(sif_avg_pred)
}
