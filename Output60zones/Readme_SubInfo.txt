SubInfo files contain tracer integrated concentrations over the water column in all MARS3D model cells that are part of the coastline. They are divided into 60 zones (named Traceur_1, ..Traceur_60).

File names (for example "subInfo2012_Apr_3w.txt") contains info on year, month and dispersal duration ("3w" = 3 weeks, "4w" = 4 weeks, "6w" = 6 weeks).

Column names are:
   Zone (1:60)
   Name (Traceur_1, ..Traceur_60): corresponds to receiver zones
   "i" , "j" are MARS3D model cells id
   "lat", "lon" are corresponding latitude and longitude
   Traceur_1, .., Traceur_60 columns corresponding to source zone