# Ranked Surface Urban Heat Island detection for hotspots and cool spots  

**DATA** folder structure:
- **LST** - Land Surface Temperature layer from Landsat imagery for 2023 with maximum 5% cloud cover, clipped to Oradea (Romania)
- **hotSpots** - layers with 10, 20 and 50 hotspots, with and without impreviousness
- **coolSpots** - layers with 10, 20 and 50 coolspots, with and without impreviousness
- **overall_Intensity** - layers with the overall intensities of the hotspots and coolspots considering at most 10, 20 and 50 hotspots, with and without impreviousness
- **persistence** - layers with the persistence of the hotspots and coolspots considering at most 10, 20 and 50 hotspots, with and without impreviousness
- **severity** - layers containing the combination of persistence and overall intensity classes for 10, 20 and 50 hotspots and cool spots with and without imperviousness

The **function.r** file contains the functions used to generate the data.
