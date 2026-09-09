# Function to perform supervised classification using Random Forest
rf_classification <- function(raster_data, classes_shapefile, study_area, split_percentage = 0.75, num_trees = 100, output_path = NULL) {
  
  # Overview:
  #   Performs supervised classification of a raster using Random Forest.
  
  # Requires:
  #   - raster_data (SpatRaster): Multi-band raster to classify.
  #   - classes_shapefile (sf): Training sample points with class labels.
  #   - study_area (sf or SpatVector): Area used to crop and mask the classification.
  #   - split_percentage (numeric): Proportion of samples used for training.
  #   - num_trees (integer): Number of trees in the Random Forest.
  #   - output_path (character, optional): Output file path.
  
  # Effects:
  #   - Extracts raster values at training sample locations.
  #   - Trains and evaluates a Random Forest classifier.
  #   - Classifies the full raster and fills gaps using a modal filter.
  #   - Crops and masks the classification to the study area.
  #   - Optionally saves the classified raster.
  
  # Returns:
  #   - A SpatRaster containing the classified vegetation classes.
  
  # Ensure the shapefile does not contain a Z-dimension
  classes_shapefile <- st_zm(classes_shapefile, drop = TRUE)
  
  # Convert class labels to factors for categorical classification
  classes_shapefile$Label <- as.factor(classes_shapefile$Label)
  
  # Extract raster values at training sample locations.
  extracted_values <- terra::extract(raster_data, classes_shapefile)
  classification_data <- data.frame(extracted_values)
  classification_data$class <- classes_shapefile$Label[classification_data$ID]  # Link to class labels
  classification_data$ID <- NULL  # Remove ID column
  classification_data <- na.omit(classification_data)  # Remove rows with missing values
  
  # Split samples into training and test sets.
  set.seed(111)  # Set seed for reproducibility
  train_index <- createDataPartition(classification_data$class, p = split_percentage, list = FALSE)
  train_data <- classification_data[train_index, ]
  test_data <- classification_data[-train_index, ]
  
  # Train the Random Forest model.
  model_rf <- ranger(x = train_data[, 1:(ncol(train_data) - 1)],
                     y = train_data$class,
                     num.trees = num_trees,
                     importance = "permutation",
                     seed = 111)
  
  # Predict test dataset and compute accuracy.
  test_predictions <- predict(model_rf, data = test_data[, 1:(ncol(test_data) - 1)])$predictions
  confusion_matrix <- table(test_data$class, test_predictions)
  test_accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  print(paste("Test Accuracy: ", round(test_accuracy * 100, 2), "%"))
  
  # Convert confusion matrix to data frame for visualization.
  confusion_df <- as.data.frame(confusion_matrix)
  colnames(confusion_df)[1:2] <- c("Var1", "Var2")
  
  # Define vegetation class labels (modify based on dataset).
  class_labels <- c("Crops", "Mixed_Vegetation", "Non-Fluorescent", "Trees")
  
  # Plot the confusion matrix.
  confusion_plot <- ggplot(confusion_df, aes(Var1, Var2, fill = Freq)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "lightgray", high = "steelblue") +
    geom_text(aes(label = Freq, color = Freq > (max(Freq) / 2)), size = 4) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
    labs(title = "Test Data Confusion Matrix", x = "Actual Class", y = "Predicted Class") +
    scale_x_discrete(labels = class_labels) +
    scale_y_discrete(labels = class_labels) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(confusion_plot)  # Display confusion matrix.
  
  # Classify the full raster using the trained model.
  classified_raster <- predict(raster_data, model_rf, fun = function(...) predict(...)$predictions, na.rm = TRUE)
  
  # Fill missing classification values using a modal filter.
  classified_raster_filled <- focal(classified_raster, w = matrix(1, 15, 15), fun = "modal", na.rm = TRUE, cores = 14)
  classified_raster[is.na(classified_raster)] <- classified_raster_filled[is.na(classified_raster)]
  
  # Crop and mask the classification to the study area.
  classified_raster <- crop(classified_raster, study_area, touches = TRUE, mask = TRUE)
  
  # Set class levels for visualization.
  levels(classified_raster)[[1]]$class <- class_labels
  
  # Plot the final classified raster.
  plot(classified_raster, main = "Classified Raster")
  
  # Optionally save the classified raster.
  if (!is.null(output_path)) {
    writeRaster(classified_raster, output_path, overwrite = TRUE)
  }
  
  return(classified_raster)
}
