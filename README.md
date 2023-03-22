# OffshoreFishWelfare

Data and code repository accompanying the article "Interactive effects of
long-term wave exposure and other stressors on fish welfare on Atlantic salmon
(Salmo salar) farms"

Szewczyk TM, Morro B, Díaz-Gil C, Gillibrand P, Hardwick JP, Davidson K, Aleynik
D, Rey Planellas S

Scripts can be run in order to fit the models and generate the results summaries
and figures.

The data represent a combination of in situ collections on Atlantic salmon farms
and data extracted and aggregated from regional hydrodynamic models.

**Column descriptions**  
`data/pub_dataset.csv`:

| Column             | Description                                                                            |
|--------------------|----------------------------------------------------------------------------------------|
| $ farm             | Farm                                                                                   |
| $ cage             | Salmon cage identifier within farm                                                     |
| $ day              | Date of week start (yyyy-mm-dd)                                                        |
| $ waveHeight       | 50-yr significant wave height (m)                                                      |
| $ fetch            | log10 wave fetch                                                                       |
| $ days_since_start | Days since start of cycle                                                              |
| $ biomass_density  | Weekly mean biomass density (kg/m^3)                                                   |
| $ temperature      | Weekly mean of daily 90th quantile (C)                                                 |
| $ waterSpeed       | Weekly mean of daily 90th quantiles (m/s)                                              |
| $ windSpeed        | Weekly mean of daily 90th quantiles (m/s)                                              |
| $ wksSinceAntiAGD  | Weeks since most recent ant-AGD treatment                                              |
| $ wksSinceMech     | Weeks since most recent mechanical anti-lice treatment                                 |
| $ wksSinceMedi     | Weeks since most recent medicinal anti-lice treatment                                  |
| $ lice_total       | Total lice (mean of 20 fish, all stages of *Lepeoptheirus salmonis* and *Caligus* spp. |
| $ agd_score        | Amoebic Gill Disease score (mean of 20 fish)                                           |
| $ lep_fe           | Lice count (mean of 20 fish, adult female *L. salmonis*)                               |
| $ mortality        | Weekly mortality rate                                                                  |
| $ lmortRate        | Logit-transformed weekly mortality rate                                                |

## Sharing/Access information

Data was derived from the following sources:

-   Wave fetch: Burrows 2012 MEPS via
    [figshare](https://figshare.com/articles/dataset/Wave_fetch_GIS_layers_for_the_UK_and_Ireland_at_200m_scale/12029682)  
-   Temperature, water speed, wind speed: WeStCOMS2 available on [SAMS thredds
    server](https://thredds.sams.ac.uk/catalog/SCOATS.html)
