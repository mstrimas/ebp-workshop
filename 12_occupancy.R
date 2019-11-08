## ----occupancy-data------------------------------------------------------
library(auk)
library(lubridate)
library(sf)
library(dggridR)
library(unmarked)
library(raster)
library(ebirdst)
library(MuMIn)
library(AICcmodavg)
library(fields)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

set.seed(1)

# ebird data
ebird <- read_csv("data/ebd_woothr_june_bcr27_zf.csv") %>% 
  mutate(year = year(observation_date),
         # occupancy modeling requires an integer response
         species_observed = as.integer(species_observed))

# modis land cover covariates
habitat <- read_csv("data/pland-elev_location-year.csv") %>% 
  mutate(year = as.integer(year))

# combine ebird and modis data
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))

# prediction surface
pred_surface <- read_csv("data/pland-elev_prediction-surface.csv")
# latest year of landcover data
max_lc_year <- pred_surface$year[1]
r <- raster("data/prediction-surface.tif")

# load gis data for making maps
map_proj <- st_crs(102003)
ne_land <- read_sf("data/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
bcr <- read_sf("data/gis-data.gpkg", "bcr") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("data/gis-data.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- read_sf("data/gis-data.gpkg", "ne_state_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()


## ----occupancy-prep-filter-----------------------------------------------
ebird_filtered <- filter(ebird_habitat, 
                         number_observers <= 5,
                         year == max(year))


## ----occupancy-prep-repeats----------------------------------------------
# subset for occupancy modeling
occ <- filter_repeat_visits(ebird_filtered, 
                            min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("latitude", "longitude", 
                                          "observer_id"))


## ----occupancy-prep-repeats-sol------------------------------------------
occ_days <- filter_repeat_visits(ebird_filtered, 
                                 min_obs = 2, max_obs = 10,
                                 n_days = 7,
                                 date_var = "observation_date",
                                 site_vars = c("latitude", "longitude", 
                                               "observer_id"))


## ----occupancy-prep-loss-sol---------------------------------------------
nrow(occ) / nrow(ebird_habitat)
n_distinct(occ$site)


## ----occupancy-prep-unmarked, class.source="livecode"--------------------
# format for unmarked, select occupancy and detection covariates
occ_wide <- format_unmarked_occu(occ, 
                                 site_id = "site", 
                                 response = "species_observed",
                                 site_covs = c("n_observations", 
                                               "latitude", "longitude", 
                                               # % deciduous forest
                                               "pland_04", 
                                               # % mixed forest
                                               "pland_05",
                                               # % cropland
                                               "pland_12",
                                               # % urban
                                               "pland_13"),
                                 obs_covs = c("time_observations_started", 
                                              "duration_minutes", 
                                              "effort_distance_km", 
                                              "number_observers", 
                                              "protocol_type",
                                              "pland_04", 
                                              "pland_05"))


## ----encounter-prep-sss, results = "hide"--------------------------------
# generate hexagonal grid with ~ 5 km betweeen cells
dggs <- dgconstruct(spacing = 5)
# get hexagonal cell id for each site
occ_wide_cell <- occ_wide %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)
# sample one checklist per grid cell
occ_ss <- occ_wide_cell %>% 
  group_by(cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  select(-cell)


## ----encounter-data-unmarked, class.source="livecode"--------------------
occ_um <- formatWide(occ_ss, type = "unmarkedFrameOccu")


## ----occupancy-model-fit, class.source="livecode"------------------------
# fit model
occ_model <- occu(~ time_observations_started + 
                    duration_minutes + 
                    effort_distance_km + 
                    number_observers + 
                    protocol_type +
                    pland_04 + pland_05
                  ~ pland_04 + pland_05 + pland_12 + pland_13, 
                  data = occ_um)
# look at the regression coefficients from the model
summary(occ_model)


## ----occupancy-model-assess, eval=FALSE, echo=1:2, class.source="dontrun"----
## occ_gof <- mb.gof.test(occ_model, nsim = 1, plot.hist = FALSE)
## print(occ_gof)
## saveRDS(occ_gof, "raw-data/woothr_occupancy-model_gof.rds")


## ----occupancy-model-assess-actual, echo=FALSE---------------------------
# read in saved gof test results 
occ_gof <- readRDS("raw-data/woothr_occupancy-model_gof.rds")
# print without chisq table
occ_gof$chisq.table <- NULL
print(occ_gof)


## ----occupancy-predict-predict, eval=FALSE, echo=1:10, class.source="livecode"----
## # make prediction for bcr 27
## occ_pred <- predict(occ_model,
##                     newdata = as.data.frame(pred_surface),
##                     type = "state")
## 
## # add to prediction surface
## pred_occ <- bind_cols(pred_surface,
##                       occ_prob = occ_pred$Predicted,
##                       occ_se = occ_pred$SE) %>%
##   select(latitude, longitude, occ_prob, occ_se)
## 
## saveRDS(pred_occ, "raw-data/woothr_occupancy-model_predictions.rds")


## ----occupancy-predict-predict-load, echo=FALSE--------------------------
pred_occ <- readRDS("raw-data/woothr_occupancy-model_predictions.rds")


## ----occupancy-predict-rasterize-----------------------------------------
r_pred <- pred_occ %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
r_pred <- r_pred[[c("occ_prob", "occ_se")]]


## ----occupancy-predict-map, fig.asp = 1.236------------------------------
# project predictions
r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")

par(mfrow = c(2, 1))
for (nm in names(r_pred)) {
  r_plot <- r_pred_proj[[nm]]
  
  par(mar = c(3.5, 0.25, 0.25, 0.25))
  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

  # occupancy probability or standard error
  if (nm == "occ_prob") {
    title <- "Wood Thrush Occupancy Probability"
    brks <- seq(0, 1, length.out = 21)
    lbl_brks <- seq(0, 1, length.out = 11) %>% 
      round(2)
  } else {
    title <- "Wood Thrush Occupancy Uncertainty (SE)"
    mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
    brks <- seq(0, mx, length.out = 21)
    lbl_brks <- seq(0, mx, length.out = 11) %>% 
      round(2)
  }
  pal <- abundance_palette(length(brks) - 1)
  plot(r_plot, 
       col = pal, breaks = brks, 
       maxpixels = ncell(r_plot),
       legend = FALSE, add = TRUE)
  
  # borders
  plot(bcr, border = "#000000", col = NA, lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  box()
  
  # legend
  par(new = TRUE, mar = c(0, 0, 0, 0))
  image.plot(zlim = range(brks), legend.only = TRUE, 
             breaks = brks, col = pal,
             smallplot = c(0.25, 0.75, 0.06, 0.09),
             horizontal = TRUE,
             axis.args = list(at = lbl_brks, labels = lbl_brks,
                              fg = "black", col.axis = "black",
                              cex.axis = 0.75, lwd.ticks = 0.5,
                              padj = -1.5),
             legend.args = list(text = title,
                                side = 3, col = "black",
                                cex = 1, line = 0))
}


## ----occupancy-model-select-dredge, class.source="livecode"--------------
# get list of all possible terms, then subset to those we want to keep
det_terms <- getAllTerms(occ_model) %>% 
  # retain the detection submodel covariates
  discard(str_detect, pattern = "psi")

# dredge all possibe combinations of the occupancy covariates
occ_dredge <- dredge(occ_model, fixed = det_terms)

# model comparison
select(occ_dredge, starts_with("psi(p"), df, AICc, delta, weight) %>% 
  mutate_all(~ round(., 3)) %>% 
  knitr::kable()


## ----occupancy-model-select-average, class.source="livecode"-------------
# select models with the most suport for model averaging
occ_dredge_95 <- get.models(occ_dredge, subset = cumsum(weight) <= 0.95)

# average models based on model weights 
occ_avg <- model.avg(occ_dredge_95, fit = TRUE)

# model coefficients
t(occ_avg$coefficients)


## ----occupancy-model-select-predict, eval=FALSE, class.source="dontrun"----
## occ_pred_avg <- predict(occ_avg,
##                         newdata = as.data.frame(pred_surface),
##                         type = "state")
## 
## # add to prediction surface
## pred_occ_avg <- bind_cols(pred_surface,
##                           occ_prob = occ_pred_avg$fit,
##                           occ_se = occ_pred_avg$se.fit) %>%
##   select(latitude, longitude, occ_prob, occ_se)


## ----occupancy-select-detection-define-----------------------------------
# define a null detection model
det_mod <- ~ time_observations_started + 
  duration_minutes + 
  effort_distance_km + 
  number_observers + 
  protocol_type ~
  pland_04 + pland_05 + pland_12 + pland_13

# define and fit candidate models
mods <- list(det_mod_null = det_mod, 
             det_mod_dec = update.formula(det_mod, ~ . + pland_04 ~ .),
             det_mod_mix = update.formula(det_mod, ~ . + pland_05 ~ .),
             global = update.formula(det_mod, 
                                     ~ . + pland_04 + pland_05 ~ .)) %>% 
  map(occu, data = occ_um)


## ----occupancy-select-detection-compare----------------------------------
mod_sel <- fitList(fits = mods) %>% 
  modSel()
mod_sel


## ----occupancy-select-detection-coef-------------------------------------
coef(occ_model) %>% 
  enframe() %>% 
  filter(str_detect(name, "pland_0"))


## ----encounter-purl, eval=FALSE, echo=FALSE------------------------------
## knitr::purl("12_occupancy.Rmd")

