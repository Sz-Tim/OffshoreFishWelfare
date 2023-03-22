# OffshoreFishWelfare

Data and code repository accompanying the article "Interactive effects of long-term wave exposure and other stressors on fish welfare on Atlantic salmon (Salmo salar) farms"

Szewczyk TM, Morro B, Díaz-Gil C, Gillibrand P, Hardwick JP, Davidson K, Aleynik D, Rey Planellas S



Scripts can be run in order to fit the models and generate the results summaries and figures.

The data represent a combination of in situ collections on Atlantic salmon farms and data extracted and aggregated from regional hydrodynamic models. 


## Description of the data and file structure

Data collected daily or weekly from the eight salmon farms between July 2018 and March 2020 were provided by Mowi Scotland Ltd. The precise range of dates available differed by farm (210-560 days). All fish samplings were non-lethal. Fish health data was collected daily or weekly by routine samplings performed by farm employees. Dead fish were counted, and 20 fish were sampled, anesthetized using Tricaine methane sulfonate (MS-222), and visually inspected for Amoebic Gill Disease (AGD) signs and the number of attached sea lice (Lepeophtheirus salmonis at each stage and Caligus spp. at the adult stage). 

In addition, the farms recorded the administration of anti-lice treatments (Alphamax, Florocol, Hydrolicer, Salmosan, SLICE, Thermolicer) and anti-AGD treatments (Freshwater, H2O2) (SEPA, https://www.sepa.org.uk/regulations/water/aquaculture/medicines-and-chemicals/). To summarise the type of stress imposed on the salmon by each treatment, anti-lice treatments were categorised as mechanical (Hydrolicer, Thermolicer) or medicinal (Alphamax, Florocol, Salmosan, SLICE). Mechanical treatments required handling of fish, while medicinal treatments were administered via the feed or as a bath. The number of weeks since each treatment type (anti-lice mechanical, anti-lice medicinal, anti-AGD) was calculated within each cage, with the mean duration used prior to the initial treatment. Farms also provided modelled fish weight, fish density, and total biomass, obtained using the software FishTalk (AKVA, https://www.akvagroup.com/fishtalk/). Due to interdependence and correlation, only biomass density (kg•m-3) was used. 

The datasets provided by each farm were subjected to quality control by checking for unlikely values, standardized, and aggregated to a weekly temporal resolution within each cage to match the weekly collection of sea lice counts and AGD prevalence. Total mortalities were summed within each week to calculate a weekly mortality rate, while biomass density was averaged. 

Farm wave exposure was assigned using the 50-year significant wave height at each location. Data from a 30-year hindcast wave model of western Scotland were extrapolated using extreme value theory to obtain the 50-year return value. The model was built using SWAN v. 41.31 (Hardwick et al., in review) on the unstructured mesh used in the WeStCOMS v.2 model domain ​(Aleynik et al., 2016, Davidson et al., 2021)​. The 50-year significant wave height is the significant wave height expected to be exceeded once in any 50-year span, where the significant wave height is the spectrally derived quantity which is close to the traditionally defined measure of the average of the largest 1/3 of measured waves (Hardwick et al., in review). SWAN is a 3rd generation spectral wave model that includes wave growth, propagation, and dissipation in coastal regions, as well as non-linear processes such as wave breaking. The 50-year significant wave height was highly correlated (r = 0.90) with wave fetch ​(Burrows, 2012)​, and was chosen as a more direct estimate of long-term wave exposure.  

Due to missing measurements and possible variation in collection methods and equipment calibration by the farms, we quantified environmental conditions at each location using the WeStCOMS model, a high-resolution hydrodynamic model on an unstructured mesh that models hourly conditions of a wide range of variables across depths within each mesh element ​(Aleynik et al., 2018, 2016)​. The WeStCOMS model was updated in 2019. Consequently, we used v.1 for dates prior to 2019 April 1 and v.2 for dates afterwards. As the mesh resolution differs between versions, values were averaged across elements within a 5 km radius to reduce localised differences. For each farm, we calculated the daily 90th quantile of temperature and water velocity at 6m depth, as well as wind velocity (Table 2). We then calculated weekly averages of these 90th quantiles for each variable to represent the environmental extremes experienced each week. 


**Column descriptions**
`data/pub_dataset.csv`:  
$ farm                 Site
$ cage                  Salmon cage identifier within farm
$ day                  Date of week start (yyyy-mm-dd)
$ waveHeight           50-yr significant wave height (m)
$ fetch                log10 wave fetch
$ days_since_start     Days since start of cycle
$ biomass_density      Weekly mean biomass density (kg/m3)
$ temperature          Weekly mean of daily 90th quantile temperatures
$ waterSpeed           Weekly mean of daily 90th quantile water speeds
$ windSpeed            Weekly mean of daily 90th quantile wind speeds
$ wksSinceAntiAGD      Weeks since most recent anti-AGD treatment
$ wksSinceMech         Weeks since most recent mechanical anti-lice treatment
$ wksSinceMedi         Weeks since most recent medicinal anti-lice treatment
$ lice_total           Total lice (mean of 20 fish, all stages of Lepeopheirus salmonis & Caligus spp.)
$ agd_score            Amoebic Gill Disease score (mean of 20 fish)
$ lep_fe               Lice count (mean of 20 fish, adult female L. salmonis)
$ mortality            Weekly mortality rate
$ lmortRate            Logit-transformed weekly mortality rate



## Sharing/Access information

Data was derived from the following sources:
  * Wave fetch: [https://figshare.com/articles/dataset/Wave_fetch_GIS_layers_for_the_UK_and_Ireland_at_200m_scale/12029682](https://figshare.com/articles/dataset/Wave_fetch_GIS_layers_for_the_UK_and_Ireland_at_200m_scale/12029682)  
  * Temperature, water speed, wind speed: [https://thredds.sams.ac.uk/catalog/SCOATS.html](https://thredds.sams.ac.uk/catalog/SCOATS.html)  


