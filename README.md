# Ranked Surface Urban Heat Island detection for hotspots and cool spots  

**SUHI-HS** detection algorithm
![suhi](https://github.com/user-attachments/assets/32ad7df0-fab8-4ce1-8c92-6f6d10a2bbb4)


**DATA** folder structure:
- **LST** - Land Surface Temperature layer from Landsat imagery for 2023 with maximum 5% cloud cover, clipped to Oradea (Romania)
- **hotSpots** - layers with 10, 20 and 50 hotspots, with and without impreviousness
- **coolSpots** - layers with 10, 20 and 50 coolspots, with and without impreviousness
- **overall_Intensity** - layers with the overall intensities of the hotspots and coolspots considering at most 10, 20 and 50 hotspots, with and without impreviousness
- **persistence** - layers with the persistence of the hotspots and coolspots considering at most 10, 20 and 50 hotspots, with and without impreviousness
- **severity** - layers containing the combination of persistence and overall intensity classes for 10, 20 and 50 hotspots and cool spots with and without imperviousness

The **functions.R** file contains the functions used to generate the data. The parametrization and then naming convetion usd for the resulted files is decribed below.

**getHotSpot** function:

* parameters:
- **wd** - full path to the working directory where the LST (Landsat B10 band) is located
- **fn** - LST (Landsat B10 band) filename without extension
- **no** - – maximum number of HSs to be detected
- **minAcceptedValueTh** - percentile defining the minimum accepted value for a cell to be part of an HS
- **minThreshold** - percentile defining the minimum average value for an HS
* naming convention: **HS_<maxHotSpotNo>_<minAcceptedValuePercentile>_<minMeanValuePercentile>-<originalLSTfilename>** 
