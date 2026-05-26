# README
## R-scripts for spatial data preparation and analysis

Our spatial unit of analysis is the 1km-by-1km grid. Since the spatial data on the different dimensions of environmental quality are not readily available at that scale, we generally assign environmental conditions to grid cells via a spatial overlap approach: 1) intersecting grid cells with the spatial units at which the environmental data is obtained, 2) grouping the resulting fragments per grid cell, and 3) processing that information into grid-level measures of environmental quality.

Data preparation and analyses were carried out in R (R Core Team, 2024). We extensively used the following packages: tidyverse (Wickham et al., 2019), sf (Pebesma, 2018), sp (Bivand, Pebesma and Gómez-Rubio, 2013), spdep (Bivand, 2022), spatialreg (Bivand, Millo and Piras, 2021), parallel (R Core Team, 2024), sandwich (Zeileis, Köll and Graham, 2020), tmap (Tennekes, 2018), modelsummary (Arel-Bundock, 2022).

All code (R scripts) is available in a public GitHub repository (https://github.com/christian-koenig/Understanding_variation_EI). Below, we provide more detail on the individual data preparation and analysis steps and highlight the R script in which each step is carried out.

## Workflow and script ordering

The scripts are numbered to reflect the order in which they are run and the logical stages of the pipeline:

- **0**: initial and small preparatory steps (restricting the grid to FUA cores; preparing OSM residential areas; preparing city-level robustness inputs such as city hall geocoding and the gentrification proxy).
- **1**: preparation of the different environmental quality measures (noise, air pollution, green space) and projection onto the 1km-by-1km grid.
- **2**: combination of the grid-level environmental measures with the neighbourhood composition data and construction of the analytical variables.
- **3**: main analyses (Moran's I, descriptive statistics, SLX estimation, city-context regressions, robustness with alternative centrality measures).
- **4**: robustness check using alternative spatial regression models.

## Restricting the spatial grid to Functional Urban Area cores
**R-Script: "00a_prep_grid1km_fua_core.R"**

The first step restricts the nationwide 1km grid from the German Federal Agency for Cartography and Geodesy to grid cells belonging to Functional Urban Area (FUA) cores, with FUA polygons obtained from the European Environment Agency. The FUA core of Dresden erroneously contained surrounding districts, which we replaced with the city boundaries of Dresden obtained from the Federal Agency for Cartography and Geodesy. We first keep grid cells that overlap with FUA cores in at least one point, intersect them with the FUA core polygons to calculate each grid cell's surface share within an FUA core, and drop grid cells with less than 50 percent of their surface located within FUA cores.

## Preparation of residential areas from OpenStreetMap
**R-Script: "00b_prep_residential_areas.R"**

We load the grid restricted to FUA cores, filter 2018 OSM land use data to "residential" entries, intersect these with grid cells, and aggregate residential areas to the grid level. This grid-level residential area data is used in the construction of the noise exposure measure.

## Preparation of city-level centrality measures for robustness check (distance to city hall)
**R-script: "00c_robustness_centrality_prep.R"**

The main city-context analyses measure the residential centrality of poor residents and foreign minorities through the within-city correlation between the grid-level group share and grid population density (population density tends to peak near the city centre in German cities). To test the robustness of our findings to this operationalisation, this script constructs an alternative centrality measure based on each grid cell's Euclidean distance to the city hall. We reproject the analytical sample to ETRS89/LAEA Europe (EPSG:3035), reduce grid cells to centroids, and geocode a manually compiled list of city hall addresses via `tidygeocoder`, matching them to municipalities through the official AGS identifier. A `tmap` map is produced for visual inspection. Looping in parallel over sample cities, we compute the Euclidean distance (in kilometres) between each grid centroid and the corresponding city hall. We then aggregate to the city level by calculating the Pearson correlation between distance to city hall and the grid-level shares of SGBII recipients, foreign residents in the labour force, academics, and grid population. The resulting city-level dataset is consumed by "03_main_analyses.R".

## Preparation of city-level gentrification indicator
**R-script: "02b_city_gentrification_acad.R"**

This script constructs a city-level proxy for gentrification dynamics, used to examine whether environmental inequality varies with the geography of neighbourhood upgrading in central locations. The reasoning is that, in gentrifying cities, the highly educated (academics) increasingly concentrate in dense, central neighbourhoods, producing a positive within-city correlation between the grid-level share of academics and grid population density. We join the FUA-core grid to municipalities by largest areal overlap, merge in Federal Employment Agency data for 2013–2017, and harmonise the AGS identifier. Using the 2017 cross-section, we keep grid cells with non-zero population and at least one working-age resident, compute the grid-level share of academics as the ratio of academics to the working-age population (employed plus unemployed), and calculate, for each city, the within-city Pearson correlation between this share and grid population density. The resulting city-level dataset is used in "03_main_analyses.R" as a covariate in the scatterplot analyses linking environmental inequality to gentrification.

## Preparation of noise data and intersection with the spatial grid
**R-script: "01_prep_noise_grid1km.R"**

Noise data from the German Federal Environment Agency is obtained via the European Environment Agency's Central Data Repository. Road and aviation noise are stored separately by federal state and agglomeration; rail traffic noise is stored in a single nationwide directory. We process each source separately, simplifying geometries with a 10m tolerance (Douglas-Peucker algorithm) to save memory and speed up processing. We then intersect each source-specific noise dataset with the grid-level residential areas prepared previously, group and aggregate by grid cell to obtain the absolute residential area affected by noise from each source, and divide by total residential area in the grid to obtain relative shares. The three source-specific datasets are merged into one grid-level dataset with separate variables for road, rail, and aviation noise exposure. Finally, we drop grid cells with less than 50 percent overlap with cities for which noise data is available.

## Preparation of air pollution data and intersection with the spatial grid
**R-script: "01_prep_airpoll_grid1km.R"**

Air pollution data from the Federal Environment Agency is obtained at a 2km grid in pollutant-specific datasets. We combine data on nitrogen dioxide (NO2), sulfur dioxide (SO2), and fine particulate matter (PM25) and construct a toxicity-weighted air pollution measure as the weighted sum of pollutant-specific annual mean concentrations, using the toxicity weights from the Federal Environment Agency based on pollutant-specific health damage estimates (Matthey and Bünger, 2020). We project this measure from the 2km grid to the 1km grid by spatial intersection, aggregating to the 1km grid as the areally-weighted sum of toxicity-weighted pollution levels (equivalent to an areally-weighted mean since grid cells are of equal size).

## Preparation of green space data and intersection with the spatial grid
**R-script: "01_prep_greenspace_grid1km.R"**

Urban Atlas 2018 data from the European Environment Agency is obtained separately for the 96 German FUAs. We loop through the FUA-specific files and retain land use polygons with the codes "14100" (green urban areas), "14200" (sports and leisure facilities), and "31000" (forests). The combined green space data is intersected with the 1km grid, and the resulting fragments are aggregated by grid cell to obtain the absolute amount of green space per grid. Division by grid area yields the grid-level surface share covered by green space.

## Combination of neighbourhood composition data with environmental quality data, and further preparation
**R-Script: "02_merge_BA_EQ_grid.R"**

After assigning grid cells to municipalities and projecting each environmental quality dataset to the 1km grid, we merge these data with the neighbourhood composition data from the Federal Employment Agency using the grid cell ID. Some city-level context measures are merged at this point (dissimilarity indices for poor and foreign minorities; correlations between grid population density and poverty share / foreign minority share; city size), generated from the grid-level data by aggregation.

At the grid level, we standardise the environmental quality variables (air pollution, noise exposure, lack of green space) and construct our index of multiple environmental burdens as the standardised sum of these standardised components. The share of poor residents is calculated as the proportion of individuals under 65 receiving social assistance for the unemployed or low-income workers; the share of foreign minorities is the share of labour force residents without German citizenship. Both are standardised. Age group shares are computed as the absolute number in each group divided by total grid population.

In the final step, we attend to the spatial structure of the data. Of 11,990 grid cells, only one has zero neighbours (no adjacent grid cells, whether inhabited or not) and is dropped. For the remaining 11,989 cells, we calculate spatial lags of the environmental quality measures. Spatial lags of the neighbourhood composition variables can only be computed for the 9,707 cells with at least one inhabited neighbour. The calculation of spatial lags is carried out within "03_main_analyses.R" since some analyses require the spatially weighted neighbours lists (e.g. Moran's I). We recombine the two datasets and only then drop uninhabited grid cells, in order to retain inhabited grid cells surrounded by uninhabited ones (i.e., cells for which spatial lags can be computed for the independent but not the dependent variables).

## Main analyses: Moran's I, descriptives, SLX estimation, city-context regressions, and centrality-based robustness
**R-script: "03_main_analyses.R"**

This script implements the paper's main empirical analyses. We first generate two distinct sets of spatial lags using queen-contiguity weights matrices: spatial lags of the environmental quality variables are based on all urban grid cells (including uninhabited ones) to incorporate spillovers from non-residential surroundings; spatial lags of socioeconomic and demographic variables are restricted to inhabited grid cells to avoid non-finite values. Both sets are merged back, and key variables are z-standardised using population-weighted standard deviations. Grid cells are then spatially joined to FUA identifiers.

We compute global Moran's I for the main outcomes and treatments using the inhabited-grid weights, and visualise the resulting estimates in a bar chart. The descriptive section computes population-weighted summary statistics for grid- and city-level variables via `modelsummary::datasummary()` and saves them as Word documents.

The main SLX analyses estimate, for each combination of outcome (standardised shares of SGBII recipients and foreign minorities) and environmental quality treatment (standardised air pollution, green space, lack of green space, noise exposure, and the multiple-environmental-burdens index), a bivariate SLX model with treatment and treatment spatial lag as regressors. We estimate two specifications, pooled and with city fixed effects, both weighted by grid population and with city-clustered standard errors. The coefficient estimates are visualised in a forest plot.

To explore cross-city variation, we re-estimate the bivariate SLX models separately for each city subsample with heteroscedasticity-robust (HC1) standard errors. City-specific estimates are visualised as paired forest-plot/choropleth-map figures for each outcome-treatment combination, with NUTS1 boundaries as cartographic background. The combined figures for the multiple-environmental-burdens index form Figure 3 of the main text.

The city-context analyses then examine which city-level features explain variation in the strength of environmental inequality. We build a city-level dataset combining the city-specific SLX estimates (with inverse-variance weights), city-level dissimilarity indices, the centrality-density link variables, city-level mean environmental exposures, the city-level group share, log city population, and population density. For each outcome-treatment pair, we estimate weighted least-squares regressions of the city-specific environmental inequality estimate on these context variables in two specifications (parsimonious and full). Results are visualised in faceted coefficient plots. We also produce scatterplots linking the city-specific environmental inequality estimates (multiple-environmental-burdens index) to the gentrification proxy from "02b_city_gentrification_acad.R".

The final section replicates the city-context analyses with the alternative city-hall-distance-based centrality measure from "00c_robustness_centrality_prep.R", using the same estimation strategy and figure layout to enable direct visual comparison.

## Robustness check: alternative spatial regression models
**R-script: "04_spatial_regressions_robustness.R"**

This script assesses the sensitivity of our main SLX results to alternative spatial econometric specifications. The motivation is that SLX imposes a particular structure on spatial dependence (spillovers only via the regressors), whereas other commonly used spatial regression models allow spatial dependence via the outcome, the error term, or both. We build a queen-contiguity neighbours list with a 10m snapping tolerance, drop the single grid cell with no contiguous neighbours, construct a row-standardised spatial weights matrix, generate spatial lags, and standardise key variables by their population-weighted standard deviation as in the main analyses.

Looping in parallel over all combinations of outcome and environmental quality treatment, we estimate five spatial regression models per combination: an unweighted bivariate SLX (`lm`), a spatial autoregressive model (SAR, `lagsarlm` with `Durbin = FALSE`), a spatial Durbin model (SDM, `lagsarlm` with `Durbin = TRUE`), a spatial error model (SEM, `errorsarlm` with `Durbin = FALSE`), and a spatial Durbin error model (SDEM, `errorsarlm` with `Durbin = TRUE`). All models are estimated on the pooled FUA-core sample using the row-standardised weights matrix. Treatment coefficients from all five model types are visualised in a single faceted coefficient plot, allowing direct comparison with the main SLX estimates.
