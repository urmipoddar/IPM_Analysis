## ************************************************************************** ##
## Author: Urmi Poddar

## Purpose:
## Interactive map of census plots
## ************************************************************************** ##
#Loading data and packages-----------------------------------------
library(tidyverse)
library(leaflet)

coords <- read_csv("Data/Plot_GPS_coords_corrected.csv")


# Data formatting------------------------------------------------------------
#  Converting degrees + decimal minutes to decimal degrees
coords <- coords %>%
  mutate(
    LAT = LATDEG + LATMIN / 60,
    LONG = -(LONGDEG + LONGMIN / 60))

#Checking the converted coordinates
coords %>%
  select(SITE, PLOTCODE, LATDEG, LATMIN, LONGDEG, LONGMIN, LAT, LONG)


#Creating interactive map-----------------------------------------
# Interactive OpenStreetMap
map <- leaflet(coords) %>%
  addProviderTiles(providers$OpenStreetMap) %>%
  addCircleMarkers(
    lng = ~LONG, 
    lat = ~LAT,
    radius = 6,
    stroke = TRUE,
    weight = 1,
    fillOpacity = 0.9,
    color = ~pal(SITE),
    label = ~PLOTCODE,
    popup = ~paste0(
      "<b>Site:</b> ", SITE,
      "<br><b>Plot:</b> ", PLOTCODE,
      "<br><b>Latitude:</b> ", round(LAT, 6),
      "<br><b>Longitude:</b> ", round(LONG, 6))) %>%
  addLegend(
    position = "bottomright",
    pal = pal,
    values = ~SITE,
    title = "Site",
    opacity = 1) %>%
  fitBounds(
    lng1 = min(coords$LONG, na.rm = TRUE),
    lat1 = min(coords$LAT, na.rm = TRUE),
    lng2 = max(coords$LONG, na.rm = TRUE),
    lat2 = max(coords$LAT, na.rm = TRUE))

map
