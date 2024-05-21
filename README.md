# Urban Heat Island

DATA folder structure:
- **LST** - Land Surface Temperature layer from Landsat imagery for 2023 with maximum 5% cloud cover, clipped to Oradea (Romania)
- **hotSpots** - layers with 10, 20 and 50 hotspots, with and without impreviousness
- **coolSpots** - layers with 10, 20 and 50 hotspots, with and without impreviousness
- **overall_Intensity** - layers with the overall intensities of the hotspots and coolspots considering at most 10, 20 and 50 hotspots, with and without impreviousness
- **persistence** - layers with the persistence of the hotspots and coolspots considering at most 10, 20 and 50 hotspots, with and without impreviousness

The **function.r** file contains all the functions used to generate the data.
