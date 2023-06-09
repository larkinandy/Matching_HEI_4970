![GitHub Logo](/Images/Shielding.png)

**Summary** <br>
This folder contains scripts used to estimate the amount maternal residences were shielded from road exposures by building and trees

- **[preprocessBuildingFootprints.ipynb](https://github.com/larkinandy/Matching_HEI_4970/blob/main/building%20and%20tree%20shielding/scripts/preprocessBuildingFootprints.ipynb
)** - convert Microsoft Bing building footprints from geojson to shp format and restrict to within 5km of maternal residences <br>
- **[preprocesCoreLogic.ipynb](https://github.com/larkinandy/Matching_HEI_4970/blob/main/building%20and%20tree%20shielding/scripts/preprocessCoreLogic.ipynb)** - restrict Core Logic records to within 5km of maternal residences.
- **[preprocessParcels.ipynb](https://github.com/larkinandy/Matching_HEI_4970/blob/main/building%20and%20tree%20shielding/scripts/preprocessParcels.ipynb
)** - join Microsoft Bing, Core Logic, and parcel datasets to create year subsets <br>
- **[shieldingScript.py](https://github.com/larkinandy/Matching_HEI_4970/blob/main/building%20and%20tree%20shielding/scripts/shieldingScript.py)** - calculate tree and building shielding. <br>
- **[calcBufferAvgs.py](https://github.com/larkinandy/Matching_HEI_4970/blob/main/building%20and%20tree%20shielding/scripts/calcBufferAvgs.py
)** - calculate buffer averages for road, wind, and shielding exposure metrics <br>
