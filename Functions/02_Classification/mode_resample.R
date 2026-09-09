# Custom function to calculate mode while excluding class 3 (Non-Fluorescent) if other classes are present
mode_resample <- function(x) {
  
  # Overview:
  #   Calculates the most frequent class value while excluding class 3
  #   when other vegetation classes are present.
  
  # Requires:
  #   - x (numeric vector): Class values, possibly containing NA values.
  
  # Effects:
  #   - Removes NA values.
  #   - Excludes class 3 when other classes are present.
  #   - Returns the most frequent remaining class or NA if none remain.
  
  # Returns:
  #   - A numeric value representing the modal class.
  
  # Remove NA values
  x <- x[!is.na(x)]
  
  # Exclude Non-Fluorescent class when other classes are present.
  if (sum(x == 3) > 0 && length(unique(x)) > 1) {
    x <- x[x != 3]
  }
  
  # Return NA if no classes remain.
  if (length(x) == 0) {
    return(NA)
  } else {
    # Return the most frequent class.
    uniq_x <- unique(x)
    return(uniq_x[which.max(tabulate(match(x, uniq_x)))])
  }
}