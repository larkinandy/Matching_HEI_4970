![GitHub Logo](/Images/WindDirection.png)

**Summary** <br>
This folder contains scripts for downloading wind vectors, converting vectors into wind direction, and calculating hours downwind from high-traffic roads

- **[downloadWindVectors.py](https://github.com/larkinandy/Matching_HEI_4970/blob/main/wind%20metrics/scripts/downloadWindVectors.py)** - download wind vectors and convert vectors into wind direction (ignoring weed speed) <br>
- **[partitionWindByYear.py](https://github.com/larkinandy/Matching_HEI_4970/blob/main/wind%20metrics/scripts/partitionWindByYear.py)** - extract hourly wind direction for the birth cohort, stratified by year.  <br>
- **[createWindMatrix.py](https://github.com/larkinandy/Matching_HEI_4970/blob/main/wind%20metrics/scripts/partitionWindByYear.py)** - calculate hours upwind/downwind for each residence, stratified by year. <br>
- **[residentialWindParallel.py](https://github.com/larkinandy/Matching_HEI_4970/blob/main/wind%20metrics/scripts/residentialWindParallel.py)** - convert hours upwind/downwind in .csv file to a rose-wind shapefile. <br>
- **[roadExposureCategoriesByYear](https://github.com/larkinandy/Matching_HEI_4970/blob/main/wind%20metrics/scripts/roadExposureCategoriesByYear.ipynb)** - derive annual road network subsets based on traffic cutofffs
