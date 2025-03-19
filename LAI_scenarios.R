#Script for creating LAI scenarios to be sent to Zhejun and Qi for input into the UCM

### prepare work environment ##############
library(dplyr)
library(tidyr)
library(readr)
library(here)
library(ggplot2)
#library(ranger)
library(mgcv)


here::here()


#load in species coefficients
shading_coef <- read_csv("iTree_shading_coefficients_species.csv")

#load in tree polygon info
all_crowns <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/treeobjects_2021_nyc.csv") %>% 
          mutate(height_m = Height * 0.3048, #convert from ft to m
            radius_m = Radius * 0.3048, #convert from ft to m
            crown_area_m2 = SHAPE_Area * 0.09290304) %>% #convert from ft2 to m2 
          filter(height_m > 2 & height_m < 40) %>%  #current height record for a tree in NYS is 50 m
          filter(radius_m > 1 )

# all_crowns %>% 
#   sample_frac(.1) %>% 
# ggplot(aes( x = height_m, y = dbh)) + geom_point()

### quantify the relationships between DBH and crown size empirically with data from NYC Tree Points ===========================

#load in the tree points database merged with the canopy polygon database
tree_points_polygons <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/tree_pheno_pointextract_polyid_all_output.csv") %>% 
    select(species, Height, SHAPE_Area, dbh, Point_ID, Radius) %>% 
    distinct() %>% 
    drop_na() %>% 
  mutate(height_m = Height * 0.3048, #convert from ft to m
         radius_m = Radius * 0.3048, #convert from ft to m
         crown_area_m2 = SHAPE_Area * 0.09290304, #convert from ft2 to m2 
         dbh_cm = dbh * 2.54) #convert from inches to cm

## create regressions for the focal taxon
taxon_subset_raw <- 
  tree_points_polygons %>% 
  filter(species == "Acer platanoides") 
   
taxon_subset_train <- sample_frac(taxon_subset_raw, 0.7) %>% mutate(train = "train")
taxon_subset_test <- anti_join(taxon_subset_raw, taxon_subset_train, by = "Point_ID")%>% mutate(train = "test")
taxon_subset <- bind_rows(taxon_subset_train, taxon_subset_test)

#visualize relationships between dbh and other variables that are available for all polygons
taxon_subset %>% 
ggplot(aes(x = dbh, y = crown_area_m2, color = height_m)) + geom_point() + theme_bw() + geom_smooth()

# #GAM ends up being a bit better than linear models and RF but has some  weird shapes
# gam_train <- mgcv::gam(dbh_cm ~ s(crown_area_m2) + s(height_m) + s(radius_m), data = taxon_subset_train)
# summary(gam_train)
# #plot(gam_train)
# 
# gam_test <- predict(gam_train, taxon_subset_test)
# taxon_subset_test <- taxon_subset_test %>%
#   mutate(predicted_dbh_cm = gam_test)
# 
# taxon_subset_test %>%
#   ggplot(aes(x = predicted_dbh_cm, y = dbh_cm)) + geom_point() + geom_smooth(method = "lm")
# summary(lm(dbh_cm ~ predicted_dbh_cm, data = taxon_subset_test))
# 
# taxon_subset_test %>%
#   ggplot(aes(x = height_m, y = predicted_dbh_cm )) + geom_point() + geom_smooth(method = "lm")


## linear model ends up being okay
lm_train <- lm(dbh_cm ~ height_m * crown_area_m2 + radius_m, data = taxon_subset_train); summary(fit)
lm_test <- predict(lm_train, taxon_subset_test)
taxon_subset_test <- taxon_subset_test %>%
  mutate(predicted_dbh_cm = lm_test)

taxon_subset_test %>%
  ggplot(aes(x = predicted_dbh_cm, y = dbh_cm)) + geom_point() + geom_smooth(method = "lm")
summary(lm(dbh_cm ~ predicted_dbh_cm, data = taxon_subset_test))

taxon_subset_test %>%
  ggplot(aes(x = height_m, y = predicted_dbh_cm, color = dbh_cm)) + geom_point() + geom_smooth(method = "lm") + 
  scale_color_viridis_c() + theme_bw()


# # random forest actually ends up not fitting quite as well
# library(ranger)
# rf <- ranger(dbh_cm ~ crown_area_m2 + height_m + radius_m, data = taxon_subset_train)
# print(rf)
# rf_test <- predict(rf, taxon_subset_test)
# taxon_subset_test <- taxon_subset_test %>%
#   mutate(predicted_dbh_cm = rf_test$predictions)
# 
# taxon_subset_test %>%
#   ggplot(aes(x = predicted_dbh_cm, y = dbh_cm)) + geom_point() + geom_smooth(method = "lm")
# summary(lm(dbh_cm ~ predicted_dbh_cm, data = taxon_subset_test))
# 
# taxon_subset_test %>%
#   ggplot(aes(x = crown_area_m2 , y = predicted_dbh_cm, color = dbh_cm)) + geom_point() + geom_smooth(method = "lm") +
#   scale_color_viridis_c() + theme_bw()


### create predictions of dbh_cm as a function of crown area, height, etc., for all of the canopy crown objects ==========
dbh_pred_focal <- predict(lm_train, all_crowns) #takes a couple min to run for GAM

all_crowns_focal_sp <- all_crowns %>% 
  mutate(dbh_cm_pred = dbh_pred_focal)
  
all_crowns_focal_sp %>% 
  sample_frac(0.05) %>% 
  ggplot(aes(x = height_m, y = dbh_cm_pred)) + geom_point()


### functions to calculate LAI from crown height, diameter, dbh (estimated) and species using i-Tree methods ##############
#sources: 
#Understanding i-Tree Summary of Programs and Methods; https://www.fs.usda.gov/nrs/pubs/gtr/gtr_nrs200-2021.pdf 
#McPherson et al. 2018; https://research.fs.usda.gov/treesearch/55873
#Nowak 2005, Nowak 2002b, Nowak 2008

calc_LAI <- function(height_m, diam_m, species){

  #look up shading coefficient from appendix 3: this is a re-working of McPherson's 2018 equations of 
  #shading coefficients as a function of dbh and species
  # this appears to be an average across all climate zones
  #estimate DBH from 
  dbh_estimate
  
  species_param_coefficient <- shading_coef$Coefficient[shading_coef$Species == species]
  species_param <- 0.0617 * ln(dbh_estimate) + 0.615 + species_param_coefficient
  
  LAI <- -4.3309 + 0.2942 * height_m + 0.7312 * diam_m + 5.7217 *species_param + -0.0148 * C_param
}


#####

#tree monoculture scenario: Platanus 