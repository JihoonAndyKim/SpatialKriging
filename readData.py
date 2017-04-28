import os
import numpy as np
import pandas as pd


lats = []
lons = []
windSpeeds = []

for filename in os.listdir(os.getcwd()+"/obs"):
    if(filename == ".DS_Store"):
        continue
    with open("obs/" + filename) as f:
        lines = f.readlines()
    if(len(lines) < 7):
        continue
    line = lines[6]
    currentline = line.split(",")

    if(len(currentline) <= 9):
        continue

    lon = currentline[2]
    lat = currentline[3]
    windSpeed = currentline[8]
    lons.append(lon)
    lats.append(lat)
    windSpeeds.append(windSpeed)

DF = pd.DataFrame({'lon' : lons, 'lat' : lats, 'wind.speed' : windSpeeds})
DF.to_csv(os.getcwd() + "obs.csv")
