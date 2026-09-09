# Function to plot a ternary colour map or colour scale
plot_ternary <- function(tern_values = NULL, type = c("map", "mesh"), res = 10,
                         hex1 = "#984EA3", hex2 = "#a6d854", hex3 = "#1C7255",
                         hexC = "#000000", blendPow = 500) {
  
  # Overview:
  #   Creates either a ternary RGB map from vegetation fractions or a
  #   ternary colour mesh showing the mixture of three vegetation classes.
  
  # Requires:
  #   - tern_values (SpatRaster): Three-layer raster containing vegetation fractions.
  #   - type (character): "map" for an RGB map or "mesh" for a colour scale.
  #   - res (integer): Resolution of the ternary colour mesh.
  #   - hex1, hex2, hex3 (character): Colours for Crops, Mixed Vegetation, and Trees.
  #   - hexC (character): Centre colour for the ternary mesh.
  #   - blendPow (numeric): Blending parameter.
  
  # Effects:
  #   - Normalizes vegetation fractions and blends class colours for RGB maps.
  #   - Creates a triangular colour mesh for ternary plots.
  #   - Produces the selected plot.
  
  # Returns:
  #   - An RGB SpatRaster when type = "map".
  #   - A ggtern object when type = "mesh".
  
  type <- match.arg(type)
  
  if (type == "map") {
    
    # Sum layers to identify zero-sum pixels.
    sum_layers <- tern_values[[1]] + tern_values[[2]] + tern_values[[3]]
    zero_mask <- sum_layers == 0
    sum_layers[zero_mask] <- 1
    
    # Normalize ternary layers.
    A <- tern_values[[1]] / sum_layers
    B <- tern_values[[2]] / sum_layers
    C <- tern_values[[3]] / sum_layers
    
    # Convert hex colours to normalized RGB vectors.
    rgb1 <- col2rgb(hex1) / 255
    rgb2 <- col2rgb(hex2) / 255
    rgb3 <- col2rgb(hex3) / 255
    
    # Blend the three class colours directly.
    R <- A * rgb1[1] + B * rgb2[1] + C * rgb3[1]
    G <- A * rgb1[2] + B * rgb2[2] + C * rgb3[2]
    B <- A * rgb1[3] + B * rgb2[3] + C * rgb3[3]
    
    # Force zero-sum pixels to black.
    R_vals <- values(R); G_vals <- values(G); B_vals <- values(B)
    zero_vals <- values(zero_mask) == 1
    R_vals[zero_vals] <- 0; G_vals[zero_vals] <- 0; B_vals[zero_vals] <- 0
    R <- setValues(R, R_vals); G <- setValues(G, G_vals); B <- setValues(B, B_vals)
    
    # Stack into RGB SpatRaster and name layers.
    rgb_raster <- c(R, G, B) * 255
    names(rgb_raster) <- c("R", "G", "B")
    
    # Plot the RGB raster.
    plotRGB(rgb_raster, r = 1, g = 2, b = 3, scale = 255, smooth = FALSE, axes = FALSE)
    
    return(invisible(rgb_raster))
  }
  
  # Convert hex colours to normalized RGB vectors.
  rgb1 <- col2rgb(hex1) / 255; rgb2 <- col2rgb(hex2) / 255
  rgb3 <- col2rgb(hex3) / 255; rgbC <- col2rgb(hexC) / 255
  
  # Build mesh of small triangles covering the ternary plot.
  tri_list <- vector("list", length = res * res * 2)
  idx <- 1
  
  for (i in 0:(res - 1)) {
    for (j in 0:(res - 1 - i)) {
      A0 <- i / res; B0 <- j / res; C0 <- 1 - A0 - B0
      A1 <- (i + 1) / res; B1 <- j / res; C1 <- 1 - A1 - B1
      A2 <- i / res; B2 <- (j + 1) / res; C2 <- 1 - A2 - B2
      A3 <- (i + 1) / res; B3 <- (j + 1) / res; C3 <- 1 - A3 - B3
      
      tri_list[[idx]] <- tibble(tri = idx, A = c(A0, A1, A2), B = c(B0, B1, B2), C = c(C0, C1, C2))
      idx <- idx + 1
      
      if (i + j < res - 1) {
        tri_list[[idx]] <- tibble(tri = idx, A = c(A1, A3, A2), B = c(B1, B3, B2), C = c(C1, C3, C2))
        idx <- idx + 1
      }
    }
  }
  
  tri_list <- tri_list[seq_len(idx - 1)]
  mesh <- bind_rows(tri_list)
  
  # Calculate centroids and blend class colours.
  mesh_centroids <- mesh %>%
    group_by(tri) %>%
    summarize(A = mean(A), B = mean(B), C = mean(C), .groups = "drop") %>%
    mutate(
      R = A * rgb1[1] + B * rgb2[1] + C * rgb3[1],
      G = A * rgb1[2] + B * rgb2[2] + C * rgb3[2],
      Bc = A * rgb1[3] + B * rgb2[3] + C * rgb3[3],
      R = pmin(1, R), G = pmin(1, G), Bc = pmin(1, Bc),
      fill = rgb(R, G, Bc)
    )
  
  # Join fill colour back to all triangle vertices.
  mesh_full <- mesh %>% left_join(select(mesh_centroids, tri, fill), by = "tri")
  
  # Plot the ternary colour mesh.
  ggtern(mesh_full, aes(x = A, y = B, z = C, group = tri, fill = fill)) +
    geom_polygon(color = NA) + scale_fill_identity() +
    scale_T_continuous(breaks = 0.5) + scale_L_continuous(breaks = 0.5) +
    scale_R_continuous(breaks = 0.5) +
    labs(T = "Mixed Vegetation", L = "Crops", R = "Trees") +
    theme_showarrows() + theme_void() +
    theme(
      tern.axis.title.T = element_text(face = "bold", size = 24, color = hex2),
      tern.axis.title.L = element_text(face = "bold", size = 24, color = hex1),
      tern.axis.title.R = element_text(face = "bold", size = 24, color = hex3),
      tern.axis.text.T = element_text(angle = 30, size = 18),
      tern.axis.text.L = element_text(angle = 30, size = 18),
      tern.axis.text.R = element_text(angle = -60, size = 18)
    )
}