## Figures replicability

The terrain figures shown in the article can be easily replicated using the provided application. Note that, in some cases, the lighting setup might differ slightly due to the reimplementation of the shading code for this release, or no shading at all used for some Figures. You can try enabling or disabling `Terrain shading` on the left panel.

The following figures are all based on the `Alps: Mont Blanc` example terrain, except for Figure 41 which used `Karakoram: Baltoro`. In the hydrology metrics, an asterisk (\*) means that the `Breach depressions` modifier must be applied to the terrain first. Figures 43, 47 and 48 overlay the rivers and ridges network on a uniformly shaded terrain, which can be set via Metrics / Local:elevation / Uniform. 

| Figure | Metric | Parameters |
|-------------|-------------|-------------|
|3 | Local:elevation / Slope (gradient norm)<br>Local:elevation / Aspect (CW from N) | - |
|4 left  | Local:elevation / Laplacian | - |
|4 right | Local:elevation / Fract Laplacian | s=0.5, n=30 |
|5 | Local:curvatures / Mean<br>Local:curvatures / Profile<br>Local:curvatures / Tangential<br>Local:curvatures / Contour | w=7 |
|8 | Local:relief / Topographic Position Index | w=33<br>w=129 |
|10 | Local:relief / Local Variance<br>Local:relief / Rugosity (area ratio)<br>Local:relief / Surface Roughness (normals spherical std)<br>Local:relief / Terrain Ruggedness Index | w=33 |
|11 | Local:elevation / Hillslope Asymmetry | w=55, 0º<br>w=55, 90º |
|12 left | Visibility / Total viewshed (out) | (EXTREMELY SLOW)<br> Rescale DEM or sample rays for approximations. |
|12 right | Visibility / Viewshed from location | Location: 324, 325<br>Observer: 1.70<br>DEM Rescale: 1 |
|14 | Visibility / Openness (positive)<br>Visibility / Openness (negative) | 8 dirs, 1000 m |
|16 | Visibility / Accessibility<br>Visibility / Sky-view factor | 2048 samples (SLOW) |
|18 left | Visibility / Sunlight | 45º |
|18 right| Visibility / DAHI | 202º |
|20 | Landforms / Curvatures | 31^2 window<br>Mask single landform: 1, 8, 9 and 14 |
|22 | Landforms / TPI landforms | 300 vs 2000 m<br>Mask single landform: 1, 4, 6 and 10 |
|23 | Orometry / Peakedness (percentile) | Radius: 1.5 km<br>Radius: 6 km |
|25 | Landforms / Black-White Top Hat | 7^2 window, t = 100m, 100m<br>7^2 window, t = 200m, 200m |
|28 | Landforms / Geomorphons | 1 km<br>Mask single landform: 5, 6, 8 and 9 |
|30 | Hydrology / Stream Area (log view) * | Flow algorithm Routing: D8<br>Flow algorithm Routing: MFD<br>(adjust colormap max to 6) |
|34 left | Hydrology / Wetness Index (TWI) * | Flow algorithm Routing: D8<br>(adjust colormap max to 16) |
|34 right | Hydrology / Depth to water * | To obtain the same output, modify in<br>Rivers net tab / Channel Initiation Threshold:<br>A \* s^0 > 200  |
|35 | Hydrology / Depth to water * | Flow algorithm Routing: MFD,<br>Left: A^1, S^1 (adjust colormap max to 1)<br>Right: A^0.4, S^1.3  |
|39 | Orometry / ORS | Radius: 1 km<br>Radius: 10 km |
|40 | Orometry / Jut<br>Orometry / Rut | Radius: 10 km |
|41 left | Orometry / Peak prototypicality | - |
|41 right| Orometry / Peakness similarity | 20 % |
|43 | Rivers net tab / Draw rivers * | Channel Initiation Threshold: A \* s^0 > 200<br>Segment width: Stream Area<br>Segment color: Uniform |
|47 | Ridges net tab / Draw ridges | drainage <= 1<br>min curvature <= -0.003<br>PPA width 5<br>Divide Tree, prom. 1m, detailed ridges |
|49 | Hydrology / Branch length max * | ^ 0.5 |
