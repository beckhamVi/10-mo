library(spdep)
library(sf)
library(ggplot2)
library(viridis)

# ============================================================================
# 1. DATOS ESPACIALES DE EJEMPLO
# ============================================================================

set.seed(123)

# Crear una cuadrícula de 10x10 polígonos
n_rows <- 10
n_cols <- 10
cell_size <- 1

# Generar coordenadas de los centroides
coords <- expand.grid(x = 1:n_cols, y = 1:n_rows)

# Crear polígonos (cuadrados)
polys <- list()
for(i in 1:nrow(coords)) {
  x <- coords$x[i]
  y <- coords$y[i]
  
  poly_coords <- matrix(c(
    x - 0.5, y - 0.5,
    x + 0.5, y - 0.5,
    x + 0.5, y + 0.5,
    x - 0.5, y + 0.5,
    x - 0.5, y - 0.5
  ), ncol = 2, byrow = TRUE)
  
  polys[[i]] <- st_polygon(list(poly_coords))
}

# Crear objeto sf
spatial_data <- st_sf(
  id = 1:nrow(coords),
  geometry = st_sfc(polys),
  crs = NA
)

# Agregar variable de interés con autocorrelación espacial
# Crear clusters de valores altos y bajos
spatial_data$valor <- rnorm(nrow(spatial_data), mean = 50, sd = 10)

# Crear hotspot en la esquina superior derecha
hotspot_ids <- which(coords$x >= 7 & coords$y >= 7)
spatial_data$valor[hotspot_ids] <- spatial_data$valor[hotspot_ids] + 30

# Crear coldspot en la esquina inferior izquierda
coldspot_ids <- which(coords$x <= 4 & coords$y <= 4)
spatial_data$valor[coldspot_ids] <- spatial_data$valor[coldspot_ids] - 20

# ============================================================================
# 2. MATRICES DE PESOS ESPACIALES
# ============================================================================

cat("\n=== MATRICES DE PESOS ESPACIALES ===\n")

# 2.1 Matriz de vecindad por contigüidad (Queen)
nb_queen <- poly2nb(spatial_data, queen = TRUE)
cat("\nVecinos por contigüidad Queen (8 vecinos):\n")
print(summary(nb_queen))

# 2.2 Matriz de vecindad por contigüidad (Rook)
nb_rook <- poly2nb(spatial_data, queen = FALSE)
cat("\nVecinos por contigüidad Rook (4 vecinos):\n")
print(summary(nb_rook))

# 2.3 Convertir a matriz de pesos (row-standardized)
W_queen <- nb2listw(nb_queen, style = "W", zero.policy = TRUE)
W_rook <- nb2listw(nb_rook, style = "W", zero.policy = TRUE)

# 2.4 Matriz de pesos por distancia (k-vecinos más cercanos)
coords_centroids <- st_coordinates(st_centroid(spatial_data))
k <- 4
nb_knn <- knn2nb(knearneigh(coords_centroids, k = k))
W_knn <- nb2listw(nb_knn, style = "W")
cat("\nVecinos por k-NN (k=4):\n")
print(summary(nb_knn))

# ============================================================================
# 3. ÍNDICE DE MORAN (AUTOCORRELACIÓN ESPACIAL GLOBAL)
# ============================================================================

cat("\n\n=== ÍNDICE DE MORAN ===\n")

# Calcular el Índice de Moran I
moran_test <- moran.test(spatial_data$valor, W_queen)
cat("\nÍndice de Moran I:\n")
print(moran_test)

# Interpretación
cat("\nInterpretación del Índice de Moran:\n")
cat(sprintf("Moran's I = %.4f\n", moran_test$estimate[1]))
cat(sprintf("p-value = %.4f\n", moran_test$p.value))
if(moran_test$p.value < 0.05) {
  if(moran_test$estimate[1] > 0) {
    cat("Conclusión: Existe autocorrelación espacial POSITIVA significativa\n")
    cat("(valores similares tienden a agruparse espacialmente)\n")
  } else {
    cat("Conclusión: Existe autocorrelación espacial NEGATIVA significativa\n")
    cat("(valores diferentes tienden a estar juntos)\n")
  }
} else {
  cat("Conclusión: No hay evidencia de autocorrelación espacial\n")
}

# Test de permutación de Monte Carlo
moran_mc <- moran.mc(spatial_data$valor, W_queen, nsim = 999)
cat("\n\nTest de Monte Carlo (999 permutaciones):\n")
print(moran_mc)

# ============================================================================
# 4. ÍNDICE DE GEARY (AUTOCORRELACIÓN ESPACIAL GLOBAL)
# ============================================================================

cat("\n\n=== ÍNDICE DE GEARY ===\n")

# Calcular el Índice de Geary C
geary_test <- geary.test(spatial_data$valor, W_queen)
cat("\nÍndice de Geary C:\n")
print(geary_test)

# Interpretación
cat("\nInterpretación del Índice de Geary:\n")
cat(sprintf("Geary's C = %.4f\n", geary_test$estimate[1]))
cat(sprintf("p-value = %.4f\n", geary_test$p.value))
cat("\nNota: Geary's C varía de 0 a 2+\n")
cat("C < 1: autocorrelación positiva\n")
cat("C = 1: aleatoriedad espacial\n")
cat("C > 1: autocorrelación negativa\n")

# Test de permutación de Monte Carlo
geary_mc <- geary.mc(spatial_data$valor, W_queen, nsim = 999)
cat("\n\nTest de Monte Carlo (999 permutaciones):\n")
print(geary_mc)

# ============================================================================
# 5. ANÁLISIS DE HOTSPOTS (LISA - Local Indicators of Spatial Association)
# ============================================================================

cat("\n\n=== ANÁLISIS DE HOTSPOTS (LISA) ===\n")

# Calcular índices de Moran locales
local_moran <- localmoran(spatial_data$valor, W_queen)

# Agregar resultados al dataset
spatial_data$local_I <- local_moran[, 1]  # Índice local
spatial_data$local_pval <- local_moran[, 5]  # p-valor

# Estandarizar valores para clasificación
spatial_data$valor_std <- scale(spatial_data$valor)
spatial_data$lag_valor_std <- lag.listw(W_queen, spatial_data$valor_std)

# Clasificar en categorías LISA
spatial_data$lisa_cat <- "No significativo"
sig_level <- 0.05

# High-High (Hotspots)
spatial_data$lisa_cat[spatial_data$valor_std > 0 & 
                        spatial_data$lag_valor_std > 0 & 
                        spatial_data$local_pval < sig_level] <- "High-High (Hotspot)"

# Low-Low (Coldspots)
spatial_data$lisa_cat[spatial_data$valor_std < 0 & 
                        spatial_data$lag_valor_std < 0 & 
                        spatial_data$local_pval < sig_level] <- "Low-Low (Coldspot)"

# High-Low (Outliers)
spatial_data$lisa_cat[spatial_data$valor_std > 0 & 
                        spatial_data$lag_valor_std < 0 & 
                        spatial_data$local_pval < sig_level] <- "High-Low (Outlier)"

# Low-High (Outliers)
spatial_data$lisa_cat[spatial_data$valor_std < 0 & 
                        spatial_data$lag_valor_std > 0 & 
                        spatial_data$local_pval < sig_level] <- "Low-High (Outlier)"

# Resumen de clusters
cat("\nDistribución de clusters LISA:\n")
print(table(spatial_data$lisa_cat))

# ============================================================================
# 6. VISUALIZACIONES
# ============================================================================

# Gráfico 1: Variable original
p1 <- ggplot(spatial_data) +
  geom_sf(aes(fill = valor), color = "white", size = 0.3) +
  scale_fill_viridis(option = "plasma", name = "Valor") +
  labs(title = "Distribución de la Variable",
       subtitle = "Valores originales") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Gráfico 2: Índice de Moran Local
p2 <- ggplot(spatial_data) +
  geom_sf(aes(fill = local_I), color = "white", size = 0.3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Moran Local I") +
  labs(title = "Índice de Moran Local",
       subtitle = "Autocorrelación espacial local") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Gráfico 3: Significancia estadística
spatial_data$significativo <- ifelse(spatial_data$local_pval < 0.05, 
                                     "Significativo (p<0.05)", 
                                     "No significativo")
p3 <- ggplot(spatial_data) +
  geom_sf(aes(fill = significativo), color = "white", size = 0.3) +
  scale_fill_manual(values = c("Significativo (p<0.05)" = "red", 
                               "No significativo" = "gray90"),
                    name = "") +
  labs(title = "Significancia Estadística",
       subtitle = "Áreas con autocorrelación significativa") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Gráfico 4: Clusters LISA (Hotspots/Coldspots)
lisa_colors <- c("High-High (Hotspot)" = "#d7191c",
                 "Low-Low (Coldspot)" = "#2b83ba",
                 "High-Low (Outlier)" = "#fdae61",
                 "Low-High (Outlier)" = "#abd9e9",
                 "No significativo" = "gray90")

p4 <- ggplot(spatial_data) +
  geom_sf(aes(fill = lisa_cat), color = "white", size = 0.3) +
  scale_fill_manual(values = lisa_colors, name = "Categoría LISA") +
  labs(title = "Clusters LISA: Hotspots y Coldspots",
       subtitle = "High-High = Hotspots, Low-Low = Coldspots") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Mostrar gráficos
print(p1)
print(p2)
print(p3)
print(p4)

# Gráfico 5: Diagrama de dispersión de Moran
moran_plot <- moran.plot(spatial_data$valor, W_queen, 
                         labels = as.character(spatial_data$id),
                         xlab = "Valor estandarizado",
                         ylab = "Lag espacial del valor",
                         main = "Diagrama de Dispersión de Moran")

# ============================================================================
# 7. EXPORTAR RESULTADOS
# ============================================================================

# Crear resumen de resultados
resultados <- data.frame(
  id = spatial_data$id,
  valor = spatial_data$valor,
  moran_local = spatial_data$local_I,
  pvalor = spatial_data$local_pval,
  categoria_lisa = spatial_data$lisa_cat
)

cat("\n\nPrimeras filas del resultado:\n")
print(head(resultados, 10))

# Opcional: guardar resultados
# write.csv(resultados, "resultados_analisis_espacial.csv", row.names = FALSE)
# st_write(spatial_data, "datos_espaciales.shp")

cat("\n\n=== ANÁLISIS COMPLETADO ===\n")
cat("\nResumen de métricas globales:\n")
cat(sprintf("- Índice de Moran I: %.4f (p = %.4f)\n", 
            moran_test$estimate[1], moran_test$p.value))
cat(sprintf("- Índice de Geary C: %.4f (p = %.4f)\n", 
            geary_test$estimate[1], geary_test$p.value))
cat(sprintf("- Hotspots detectados: %d\n", 
            sum(spatial_data$lisa_cat == "High-High (Hotspot)")))
cat(sprintf("- Coldspots detectados: %d\n", 
            sum(spatial_data$lisa_cat == "Low-Low (Coldspot)")))

