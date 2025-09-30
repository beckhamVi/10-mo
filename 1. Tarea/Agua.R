library(haven)
library(dplyr)
library(sf)
library(leaflet)

# 1. Leer datos
caratula <- read_sav("E:\\10\\Estadistica espacial\\GitHub\\CARATULA.sav")
cap100a_04 <- read_sav("E:\\10\\Estadistica espacial\\GitHub\\01_CAP100A_04.sav")


# 2. Filtrar solo Puno
caratula_puno <- caratula %>% filter(NOMBREDD == "PUNO")
cap100a_puno <- cap100a_04 %>% filter(NOMBREDD == "PUNO")

# 3. Unir coordenadas
puno_merged <- cap100a_puno %>%
  left_join(caratula_puno %>% 
              select(ANIO, CCDD, CCPP, CCDI, NSEGM, ID_PROD, UA,
                     LATITUD, LONGITUD),
            by = c("ANIO", "CCDD", "CCPP", "CCDI", "NSEGM", "ID_PROD", "UA"))

# 4. Convertir a sf
puno_sf <- st_as_sf(puno_merged,
                    coords = c("LONGITUD", "LATITUD"),
                    crs = 4326)

# 5. Diccionarios de etiquetas
fuente_agua <- c("1"="Río", "2"="Manantial/puquio", "3"="Pozo",
                 "4"="Represa", "5"="Peq. reservorio", "6"="Otro")
sistema_riego <- c("1"="Exudación", "2"="Goteo", "3"="Microaspersión",
                   "4"="Aspersión", "5"="Multicompuertas",
                   "6"="Mangas", "7"="Gravedad", "8"="Otro")

# 6. Mapa A: Fuente de agua (P120)
puno_agua <- puno_sf %>% filter(!is.na(P120)) %>% slice_sample(n = 2000)

mapa_agua <- leaflet(puno_agua) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(
    radius = 3,
    color = ~case_when(
      P120 == 1 ~ "blue",
      P120 == 2 ~ "green",
      P120 == 3 ~ "orange",
      P120 == 4 ~ "purple",
      P120 == 5 ~ "brown",
      P120 == 6 ~ "red",
      TRUE ~ "gray"
    ),
    popup = ~paste0(
      "<b>Cultivo:</b> ", P115_NOM,
      "<br><b>Fuente agua:</b> ", fuente_agua[as.character(P120)]
    )
  ) %>%
  addLegend(
    position = "bottomright",
    colors = c("blue","green","orange","purple","brown","red"),
    labels = c("Río","Manantial","Pozo","Represa","Peq. reservorio","Otro"),
    title = "Fuente de agua"
  )

# 7. Mapa B: Sistema de riego (P121)
puno_riego <- puno_sf %>% filter(!is.na(P121)) %>% slice_sample(n = 2000)

mapa_riego <- leaflet(puno_riego) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(
    radius = 3,
    color = ~case_when(
      P121 == 1 ~ "darkblue",   # Exudación
      P121 == 2 ~ "green",      # Goteo
      P121 == 3 ~ "orange",     # Microaspersión
      P121 == 4 ~ "purple",     # Aspersión
      P121 == 5 ~ "brown",      # Multicompuertas
      P121 == 6 ~ "pink",       # Mangas
      P121 == 7 ~ "red",        # Gravedad
      P121 == 8 ~ "gray",       # Otro
      TRUE ~ "black"
    ),
    popup = ~paste0(
      "<b>Cultivo:</b> ", P115_NOM,
      "<br><b>Sistema riego:</b> ", sistema_riego[as.character(P121)]
    )
  ) %>%
  addLegend(
    position = "bottomright",
    colors = c("darkblue","green","orange","purple","brown","pink","red","gray"),
    labels = c("Exudación","Goteo","Microaspersión","Aspersión",
               "Multicompuertas","Mangas","Gravedad","Otro"),
    title = "Sistema de riego"
  )

# 8. Mostrar mapas (se muestran por separado en RStudio/Quarto)
mapa_agua
mapa_riego
