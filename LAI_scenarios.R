#Script for creating LAI scenarios to be sent to Zhejun and Qi for input into the UCM

### prepare work environment ##############
library(dplyr)
library(tidyr)
library(readr)
library(here)
library(ggplot2)
#library(ranger)
#library(mgcv)
library(terra)
library(sf)


here::here()


#load in species coefficients
shading_coef <- read_csv("iTree_shading_coefficients_species.csv")

#load in tree object polygons for all of NYC
all_crowns <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/treeobjects_2021_nyc.csv") %>% 
          mutate(height_m = Height * 0.3048, #convert from ft to m
            radius_m = Radius * 0.3048, #convert from ft to m
            crown_area_m2 = SHAPE_Area * 0.09290304) %>% #convert from ft2 to m2 
          filter(height_m > 2 & height_m < 40) %>%  #current height record for a tree in NYS is 50 m
          filter(radius_m > 1 )

# all_crowns %>% 
#   sample_frac(.1) %>% 
# ggplot(aes( x = height_m, y = dbh)) + geom_point()

#load in LCZ for aggregating to the same scale as the UCM
bx <- rast("C:/Users/dsk273/Box/collaborations/trees heat and health/new_york_CONUS_LCZ_map_dlm_chm2021_bronx_mean.tif")
#plot(bx)

#load in tree canopy cover
##NEED TO GET THIS FROM DAVE


### new approach: sample from 2013 i-Tree census and plug into the canopy crown polygons =========================================

#load in the i-Tree Eco 2013 processed tree data that includes LAI estimates
it13 <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/2013_itree_plots/itree_2013_tree_measurements_processed.csv") %>% 
  janitor::clean_names() %>% 
  mutate(common_name = species_name,
         height_m = crown_height_ft / 3.281, 
         canopy_cover_m2 = canopy_cover_ft_2 / 10.764,
         dbh_cm = dbh_in * 2.54,
         basal_area_m2 = basal_area_ft_2 / 10.764,
         lai = leaf_area_index) 

#add scientific names
itree_speciesnames <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/2013_itree_plots/itree_SpeciesList.csv") %>% 
  janitor::clean_names() %>% 
  mutate(genus = stringr::word(scientific_name))
it13 <- left_join(it13, itree_speciesnames)

#which genera are most common
itree13_common_genera <- it13 %>%  
  group_by(genus) %>% 
  summarize(n_per_genus = n(),
            ba_sum_genus = sum(basal_area_m2),
            canopy_cover_genus = sum(canopy_cover_m2)) #%>% arrange(-n_per_species) %>%  top_n(15)
it13 <- left_join(it13, itree13_common_genera)

#weighted mean of LAI per canopy cover
test <- 
  it13 %>% ungroup() %>% 
  dplyr::group_by(genus) %>% 
  summarize(  LAI_wm = weighted.mean(lai, canopy_cover_m2),
              n = n()) %>% 
  arrange(-n) %>% 
  filter(n > 10)


#some data vis
it13 %>% filter(n_per_genus > 9) %>% 
  ggplot(aes(x = lai)) + geom_histogram() + facet_wrap(~genus)
#test <- 
  it13 %>% group_by(genus) %>% summarize(mean_lai = mean(lai), n_per_genus = mean(n_per_genus)) %>% arrange(-mean_lai) %>% print(n = 100)


### set focal taxon 
focal_genus <- "Acer"
focal_taxon_lai <- it13 %>%  filter(genus == focal_genus) %>% select(lai) %>% unlist(.)




### add LAI to crown objects 
all_crowns_sim <- all_crowns %>% 
  mutate(sim_lai = sample(focal_taxon_lai, size = nrow(all_crowns), replace = TRUE))

### create weighted mean of LAI per LCZ pixel 
mean(focal_taxon_lai)

### create dataframe of unique values for Zhejun 














### quantify the relationships between DBH and crown size empirically with data from itree 2013 census ===========================

#load in the i-Tree Eco 2013 tree data for checking dbh vs height vs canopy width
itree13 <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/2013_itree_plots/itree_2013_tree_measurements.csv") %>% 
  janitor::clean_names() %>% 
  mutate(height_m = total_height_ft * 0.3048, #convert from ft to m
         radius_m = ((crown_width_n_s_ft + crown_width_e_w_ft )/2) * 0.3048) #convert from ft to m

#function to convert dbh in inches to basal area in m2
dbh_in_to_ba_m2 <- function(x){
  dbh_cm <- 2.54 * x
  ba_m2 <- 0.00007854 * dbh_cm^2 
}

dbh_cm <- itree13 %>% 
  select(contains("dbh") & contains("in")) %>% 
  dbh_in_to_ba_m2(.) %>% 
  mutate(ba_m2_total = rowSums(., na.rm = TRUE),
         dbh_cm = sqrt(ba_m2_total/0.00007854)) %>% 
  select(dbh_cm)

itree13 <- bind_cols(itree13, dbh_cm)

itree13_common_sp <- itree13 %>%  
  group_by(species) %>% 
  summarize(n = n()) %>% 
  arrange(-n) %>% 
  top_n(15)

itree13 %>% 
  filter(species %in% itree13_common_sp$species) %>% 
  ggplot(aes(x = dbh_cm, y = height_m)) + geom_point() + geom_abline(intercept = 0, slope = 1) + geom_smooth() + theme_bw() + 
  facet_wrap(~species)

itree13 %>% 
  filter(species %in% itree13_common_sp$species) %>% 
  ggplot(aes(x = dbh_cm)) + geom_histogram() + theme_bw() + 
  facet_wrap(~species)


### quantify the relationships between DBH and crown size empirically with data from NYC Tree Points ===========================
#load in the tree points database merged with the canopy polygon database
tree_points_polygons <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/tree_pheno_pointextract_polyid_all_output.csv") %>% 
    select(species, Height, SHAPE_Area, dbh, Point_ID, Radius) %>% 
    distinct() %>% 
    drop_na() %>% 
  mutate(height_m = Height * 0.3048, #convert from ft to m
         radius_m = Radius * 0.3048, #convert from ft to m
         crown_area_m2 = SHAPE_Area * 0.09290304, #convert from ft2 to m2 
         dbh_cm = dbh * 2.54) %>%  #convert from inches to cm
  filter(dbh_cm > 1 & dbh_cm < 200) %>% #remove some obvious mistakes
  filter(height_m < 50) %>% 
  filter(height_m/dbh_cm < 1) #this is to remove off some obviously wrong individuals

#ggplot(tree_points_polygons, aes(x = dbh_cm, y = radius_m)) + geom_jitter(alpha = 0.01) + theme_bw() + scale_x_continuous(limits = c(0,25))
#summary(tree_points_polygons)


## create regressions for the focal taxon
taxon_subset_raw <- 
  tree_points_polygons %>% 
  filter(species == "Acer rubrum") 
   
taxon_subset_train <- sample_frac(taxon_subset_raw, 0.7) %>% mutate(train = "train")
taxon_subset_test <- anti_join(taxon_subset_raw, taxon_subset_train, by = "Point_ID")%>% mutate(train = "test")
taxon_subset <- bind_rows(taxon_subset_train, taxon_subset_test)



#visualize relationships between dbh and other variables that are available for all polygons
taxon_subset %>% 
ggplot(aes(x = dbh_cm, y = height_m, color = radius_m)) + geom_point() + theme_bw() + geom_smooth()

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