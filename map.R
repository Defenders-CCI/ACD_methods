centroids <- read.csv(file = 'C:/Users/mevans/Downloads/ACD_centroids.csv', header = TRUE)
centroids <- mutate(centroids, coords = str_extract(.geo, '\\[.*\\]'))
coords <- str_split_fixed(centroids$coords, ',', 2)
centroids$lon <- str_remove(coords[,1], '\\[')
centroids$lat <- str_remove(coords[,2], '\\]')


g <- list(
  scope = 'usa',
  projection = list(type = 'albers usa'),
  showland = TRUE,
  landcolor = toRGB("gray95"),
  subunitcolor = toRGB("gray85"),
  countrycolor = toRGB("gray85"),
  countrywidth = 0.5,
  subunitwidth = 0.5
)

map <- plot_geo(centroids, lat = ~lat, lon = ~lon,
                symbol = ~Habitat, marker = list(color = 'black', size = 10))%>%
  layout(geo = g)
map

bar <- plot_ly(group_by(centroids, Habitat)%>%summarize(count = n()),
               x = ~count, y= ~Habitat, type = 'bar', orientation = 'h',
               marker = list(color = 'grey', line = list(color = 'black', width = 2)))%>%
  layout(xaxis = list(
    titlefont = list(size = 14, color = 'black', family = 'Serif'),
    tickfont = list(size = 18, color = 'black', family = 'Serif'),
    showgrid = FALSE
    ),
    yaxis = list(
      titlefont = list(size = 14, color = 'black', family = 'Serif'),
      tickfont = list(size = 18, color = 'black', family = 'Serif'),
      showgrid = FALSE,
      zeroline = TRUE
    )
  )
bar

# add orca command line utility to R environmental path
Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "C:\\Users\\mevans\\AppData\\Local\\Programs\\orca", sep = .Platform$path.sep))

orca(map, file = 'map.png', format = tools::file_ext('png'), scale = 5)
orca(bar, file = 'bar.png', format = tools::file_ext('png'), scale = 5)
