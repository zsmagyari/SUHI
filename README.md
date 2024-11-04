## Ranked Surface Urban Heat Island detection for hotspots and cool spots


![suhi](https://github.com/user-attachments/assets/32ad7df0-fab8-4ce1-8c92-6f6d10a2bbb4)
<p align="center"> The pseudocode algorithm for Ranked Surface Urban Heat Island detection </p>

---

**DATA** folder structure:
- **LST** - Land Surface Temperature layer from Landsat imagery for 2023 with maximum 5% cloud cover, clipped to Oradea (Romania)
- **hotSpots** - layers with 10, 20 and 50 hotspots, with and without impreviousness
- **coolSpots** - layers with 10, 20 and 50 coolspots, with and without impreviousness
- **overall_Intensity** - layers with the overall intensities of the hotspots and coolspots considering at most 10, 20 and 50 hotspots, with and without impreviousness
- **persistence** - layers with the persistence of the hotspots and coolspots considering at most 10, 20 and 50 hotspots, with and without impreviousness
- **severity** - layers containing the combination of persistence and overall intensity classes for 10, 20 and 50 hotspots and cool spots with and without imperviousness

---

The **functions.R** file contains the functions used to generate the data. The parametrization and then naming convention used for the resulted files is decribed below.

**getHotSpots** function:
* parameters:
  - **wd** - full path to the working directory where the LST (Landsat B10 band) is located
  - **fn** - LST (Landsat B10 band) filename without extension
  - **no** - maximum number of HSs to be detected
  - **minAcceptedValueTh** - percentile defining the minimum accepted value for a cell to be part of an HS
  - **minThreshold** - percentile defining the minimum average value for an HS
* naming convention: **HS_**«_maxHotSpotNo_»**_**«_minAcceptedValuePercentile_»**_**«_minMeanValuePercentile_»**-**«_originalLSTfilename_» 

---

**hotspotPersistence** function:
* parameters:
  - **wd** - full path to working directory where a folder named spots should contain all the files resulting from the above-presented SUHI-HS detection algorithm
  - **type** - a string parameter representing a filter to differentiate between situations when imperviousness was and wasn’t considered for SUHI-HSs detection
  - **limit** - maximum number of SUHI-HS to be considered for the assessment
  - **partNo** - number of classes for reclassification
* naming convention: **T_nr**«_numberOfConsideredFiles_»**_limit**«_maxSpotNumber_»**_part**«_numberOfIntervals_»**_**«_typeParamValue_»

---

**hotspotOverallIntensity** function:
* parameters:
  - **wd** - full path to working directory where a folder named spots should contain all the files resulting from the above-presented SUHI-HS detection algorithm
  - **type** - a string parameter representing a filter to differentiate between situations when imperviousness was and wasn’t considered for SUHI-HSs detection
  - **limit** - maximum number of SUHI-HSs to be considered
  - **partNo** - number of intervals/classes
  - **bitNo** - number of bits to represent a hotspot rank
* naming convention: **I_nr**«_numberOfConsideredFiles_»**_limit**«_maxSpotNumber_»**_part**«_numberOfIntervals_»**_**«_typeParamValue_» 
