#This script is for 

library(rnpn)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(readr)
library(prism)
library(sf)
library(purrr)

install.packages("sf")

### select a few top anemophilous taxa from NPN #####################################
#anemophilous angiosperms of potential allergenic concern
#npn_species() %>%  filter(functional_type == "Deciduous broadleaf") %>% arrange(genus)
acer_species_list <- c(777,1843,59,778,1,2,1591,60,779,780,3,781,61,1199)
alnus_species_list <- c(62,63,319)
betula_species_list <- c(97, 1439, 98, 430, 1850, 1339, 1851, 99, 1805)
carya_species_list <- c(1176, 1177, 824, 67, 68)
celtis_species_list <- c(829, 1342, 1605)
corylus_species_list <- c(71, 72, 2241)
fraxinus_species_list <- c(74,872,873,75,1350)
ginkgo_species_list <- c(884)
gleditsia_species_list <- c(885)
juglans_species_list <- c(1353, 80, 2143, 2260, 2259)
morus_species_list <- c(2007, 2266)
platanus_species_list <- c(2271, 970)
populus_species_list <- c(1361,320,976,977,1188,27,1481)
quercus_species_list <- c(705,100,1365,757,1870,987,1690,1484,988,316,297,1485,1190,765,1486,
                          301,704,101,1691,1212,989,1366,102,1756,1213,1755,1487,1159,305)
salix_species_list <- c(1007, 717, 1875, 1494, 20166, 293, 322, 1493, 1006, 1163, 77, 1371, 1009, 2278, 1008, 1372, 1010, 1876)
tilia_species_list <- c(93, 1775, 1776, 1777)
ulmus_species_list <- c(1192,1048,1049,1215,1216)
zelkova_species_list <- c(1374)


list_all_focal_taxa <- c(acer_species_list, alnus_species_list, betula_species_list, carya_species_list, celtis_species_list, corylus_species_list,
                         fraxinus_species_list, ginkgo_species_list, gleditsia_species_list, juglans_species_list, morus_species_list,
                         platanus_species_list, populus_species_list, quercus_species_list, salix_species_list, tilia_species_list,
                         ulmus_species_list, zelkova_species_list
                         )

###download and process data #####################################
# to assess which phenophases to use, I downloaded all of the relevant phenophase classes for several taxa
# pheno_class_ids = c(1,2,3,4,6,7,8)
# and then compared which ones were available to the phenophase definitions (e.g., p. 36 of the Phenophase Primer)
# https://usanpn.org/files/education/Phenophase_Primer_for_Flowering_Plants_Dec2024.pdf 
# npn_direct %>% 
#   group_by(phenophase_id, phenophase_description) %>% 
#   summarize(n = n()) %>% 
#   arrange(-n)

npn_direct <- npn_download_status_data(
  request_source = 'Daniel Katz, Cornell University',
  species_ids = list_all_focal_taxa, # list_all_focal_taxa,
  years = c(as.character(2012:2025)), #years to include
  phenophase_ids = c( 467, # increasing leaf size phenophase
                          # For deciduous tree and shrub species that have a single large flush of leaves at the beginning of the growing season, 
                          # we observe “Increasing leaf size”. This phenophase starts with the first unfolded leaf and ends when the majority of 
                          # leaves from that flush have expanded to their full mature size. It represents "green-up"—the period of time it takes the 
                          # plant to produce a full canopy of leaf tissue and build its photosynthetic capacity. 
                      501, 502), #angiosperms: 501 == "Open flowers", 502 == "Pollen release (flowers)" 
   additional_fields = c("Observed_Status_Conflict_Flag", "partner_group", "observedby_person_id",  "Observation_Comments",
                         "AGDD", "gdd", "Daylength","Greenup_0", "MidGreenup_0", "Maturity_0") #
  #additional_fields documentation : https://docs.google.com/document/d/1yNjupricKOAXn6tY1sI7-EwkcfwdGUZ7lxYv7fcPjO8/edit?tab=t.0
 # six_leaf_layer = TRUE, six_bloom_layer = TRUE # agdd_layer = 32 
)
#npn_direct_raw <- npn_direct
#write_csv(npn_direct, "C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/npn_flow_leaf_download250507.csv")
#npn_direct <- read_csv("C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/npn_flow_leaf_download250507.csv")

#clean up and organize data
npn_direct <- subset(npn_direct, observed_status_conflict_flag == "-9999")
npn_direct <- filter(npn_direct, !(partner_group %in% c("CSU Chico NSCI 102", "SCMF Naturalists", 
                                                        "Sycamore Canyon Wilderness Park - Riverside",
                                                        "Marist College", "Sam Hughes Neighborhood Association",
                                                        "UNCO BIO 111", "Maricopa Cooperative Extension",
                                                        "Pima County Extension", "Lasell College",
                                                        "UofL campus", "Ursinus College", "U of A Campus Arboretum",
                                                        "RMC Campus Phenology", "SUNY Geneseo", "AZ Project WET")))
npn_direct <- npn_direct %>% 
  mutate(year_obs = year(observation_date)) %>% 
  mutate(intensity_p = case_when(
                        phenophase_status == 0  ~ 0,
                        intensity_value == "Less than 25%" ~ 0.125,
                        intensity_value == "25-49%" ~ (0.25+0.49)/2,
                        intensity_value == "50-74%" ~ (0.5+0.74)/2,
                        intensity_value == "75-94%" ~ (0.75+0.94)/2,
                        intensity_value == "95% or more" ~ 0.97,
                        phenophase_status == 1  ~ NA)) %>%  #
  arrange(individual_id, observation_date)
#test <- npn_direct %>% filter(individual_id == 289346 & phenophase_id == 467)

#add in how many observations each observer made
obs_per_observer <- npn_direct %>% 
  group_by(observedby_person_id) %>% 
  summarize(n_obs_per_person = n()) %>% 
  arrange(-n_obs_per_person) #summary() #ggplot(aes(x = n_obs_per_person)) + geom_histogram()
  
npn_direct <- left_join(npn_direct, obs_per_observer)

# assess an individual observer
# test <- npn_direct %>% 
#   filter(observedby_person_id == "2002") %>% 
#   filter(individual_id == 4318) %>% 
#   filter(phenophase_id == 467) %>% 
#   ggplot(aes(x = ymd(observation_date), y = intensity_p, color = phenophase_description)) + geom_point()


#take the average date of leaf expansion for each individual tree in each year
indiv_leafout <- npn_direct %>% filter(phenophase_id == 467) %>% 
  filter(phenophase_status == 1) %>% 
  filter(day_of_year < 173) %>%  #solstice
  group_by(genus, species, latitude, longitude, individual_id, year_obs) %>% 
  slice(which.min(abs(intensity_p - 0.5))) %>% 
  filter(intensity_p > 0.2 & intensity_p < 0.8)

indiv_leafout_summary <- npn_direct %>% filter(phenophase_id == 467) %>% 
  filter(phenophase_status == 1) %>% 
  filter(day_of_year < 173) %>%  #solstice
  filter(intensity_p > 0.2 & intensity_p < 0.8) %>% 
  group_by(genus, species, individual_id, year_obs) %>% 
  summarize(leaf_mean = mean(day_of_year),
            leaf_early = min(day_of_year),
            leaf_late = max(day_of_year),
            leaf_nobs = n()) %>% 
  mutate(leafout_duration = leaf_late - leaf_early)
  
indiv_leafout <- left_join(indiv_leafout, indiv_leafout_summary) %>% 
  filter(leafout_duration < 15)

#take the average date of open flowers for each individual in each year
indiv_flow <- npn_direct %>% filter(phenophase_id == 501) %>% 
  filter(phenophase_status == 1) %>% 
  filter(day_of_year < 173) %>%  #solstice
  group_by(genus, species, latitude, longitude, individual_id, year_obs) %>% 
  slice(which.min(abs(intensity_p - 0.5))) %>% 
  filter(intensity_p > 0.2) #"what percentage of all fresh flowers on the plant are open?"

indiv_flow_summary <- npn_direct %>% filter(phenophase_id == 501) %>% 
  filter(phenophase_status == 1) %>% 
  filter(day_of_year < 173) %>%  #solstice
  group_by(genus, species, individual_id, year_obs) %>% 
  filter(intensity_p > 0.2) %>% 
  summarize(flow_mean = mean(day_of_year),
            flow_early = min(day_of_year),
            flow_late = max(day_of_year),
            flow_nobs = n()) %>% 
  mutate(flow_duration = flow_late - flow_early)

indiv_flow <- left_join(indiv_flow, indiv_flow_summary) %>% 
  filter(flow_duration < 15)

indiv_flow_join <- indiv_flow %>% 
  select(year_obs, individual_id, latitude, longitude, genus, species, common_name, flow_mean, flow_early, flow_late, flow_nobs, flow_duration, flow_intensity_p = intensity_p)

#join leaf expansion and flowering and retain observations that have both
lf <- left_join(indiv_leafout, indiv_flow_join) %>% 
  filter(!is.na(flow_mean)) 


#take the average date of pollen release for each individual in each year
indiv_pol_release <- npn_direct %>% filter(phenophase_id == 502) %>%
  filter(phenophase_status == 1) %>%
  filter(intensity_value %in% c("Little","Some", "Lots")) %>% 
  filter(day_of_year < 173) %>%  #solstice
  group_by(genus, species, latitude, longitude, individual_id, year_obs) %>%
  summarize(pol_mean = mean(day_of_year),
            pol_early = min(day_of_year),
            pol_late = max(day_of_year),
            pol_dif = pol_late - pol_early,
            pol_nobs = n()) %>% 
  filter(pol_dif < 15)  #exclude observations that implausibly claim a tree was releasing pollen for over 2 wks 
        # indiv_pol_release %>% 
        #   group_by(pol_dif) %>% 
        #   summarize(nobs = n())


###extract temperature data for each location in each year ####################
prism_set_dl_dir("C:/Users/dsk273/Documents/prism")
#get_prism_monthlys(years = 2024:2025, type = "tmean") #updating my local dataset
#get_prism_dailys(minDate = ymd("2025-05-01"), maxDate = ymd("2025-05-05"), type = "tmean") #updating my local dataset


#prepare data to extract
lf_extract <- lf %>% ungroup() %>% 
  mutate(year_obs = year(observation_date),
         too_recent = case_when(year_obs == 2025 & flow_mean > 120 ~ "too recent",
                                year_obs == 2025 & leaf_mean > 120 ~ "too recent",
                                TRUE ~ "not too recent")) %>% 
  filter(too_recent == "not too recent") %>% #removing observations from this month
  dplyr::select(individual_id, longitude, latitude, year_obs, flow_mean, leaf_mean) %>% 
  distinct() 


### function to extract temperature for days around flowering/leaf out
extract_temp_daily <- function(individual_id, year_obs, flow_mean, leaf_mean, longitude, latitude){ 
  individual_id_focal <- individual_id
  longitude_focal <- longitude
  latitude_focal <- latitude
  yr_focal <- year_obs
  #year_obs <- lf_sf$year_obs[1]; flow_mean <- lf_sf$flow_mean; leaf_mean <- lf_sf$leaf_mean[1]; individual_id_focal <- lf_sf$individual_id[1]
  date_min_doy <- max(flow_mean) - 30 #can substitute in flowering or leaf or the max of both here
  date_max_doy <- max(flow_mean)
  date_min <- ymd(paste0(yr_focal, "/1/1")) + date_min_doy
  date_max <- ymd(paste0(yr_focal, "/1/1")) + date_max_doy
  
  #print(date_min); print(date_max)
  tmean_rast_d <- prism_archive_subset(temp_period = "daily", type = "tmean", minDate = date_min, maxDate = date_max)
  tmean_rast2_d <- pd_stack(tmean_rast_d)
  r_mean <- raster::calc(tmean_rast2_d, mean) #raster::plot(r_mean)
  
  lf_sub <- filter(lf_extract, year_obs == yr_focal, individual_id == individual_id_focal) 
  #print(lf_sub)
  lf_sub_sf <- lf_sub %>% st_as_sf(coords = c( "longitude", "latitude"), crs = 4326) 
  
  tmean_data <- unlist(raster::extract(x = r_mean,
                                       y = lf_sub_sf)) %>% as.data.frame()

  lf_sub2 <- lf_sub %>% mutate(tmean = as.numeric(unlist(tmean_data)))

  print(tmean_data)
  return(lf_sub2)
}

#apply function to calculate temp for each observation individually
#lf_test <- lf_extract[20:21,]
spring_temp <- purrr::pmap_dfr(.l = lf_extract, .f = extract_temp_daily)

npn_active_flow <- left_join(npn_active_flow, spring_temp)



### extract monthly spring temperature in each year, one year at a time
tmean_data_all_years_list <- list()
yr_list <- 2012:2025
for(i in 1:length(yr_list)){
#create a raster of mean spring temperature for a particular year
  tmean_rast_d <- prism_archive_subset(temp_period = "monthly", type = "tmean", mon = c(1,4), years = yr_list[i])
  tmean_rast2_d <- pd_stack(tmean_rast_d)
  tmean_focal_year <- raster::calc(tmean_rast2_d, mean) #raster::plot(r_mean)

#convert the points to extract to sf
  lf_yr <- lf_extract %>% filter(year_obs == yr_list[i])
  lf_yr_sf <- lf_yr %>% st_as_sf(coords = c( "longitude", "latitude"), crs = 4326) 
  
#extract temperature for all points in that year
  tmean_data_focal_year <- mutate(lf_yr, tmean_spring = unlist(
      raster::extract(x = tmean_focal_year, y = lf_yr_sf)))
  
#combine results into a list
  tmean_data_all_years_list[[i]] <- tmean_data_focal_year
}

#merge all the years together
tmean_data_all_years_spring <- bind_rows(tmean_data_all_years_list)
  
#join back into the parent dataset
lf <- left_join(lf, tmean_data_all_years_spring)


### joining various datasets for data visualization
lfp <- left_join(lf, indiv_pol_release)


#write_csv(lfp, "C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/npn_lfp_250508.csv")
#lfp <- read_csv("C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/npn_lfp_250508.csv")



#### data visualization ########################################################
# test <- lf %>% mutate(lf_dif = leaf_mean - flow_mean) %>% 
#   arrange(lf_dif) %>% 
#   filter(longitude > -90 & longitude < -60) %>% 
#   filter(tmean_spring > 3 & tmean_spring < 9) 

# test2 <- filter(npn_direct, individual_id == 132692) %>% 
#   filter(phenophase_id == 467)

lf %>% 
  #filter(latitude > 40.7 - 3 & latitude < 40.7 + 3) %>% 
  #filter(longitude > -74 - 3 & longitude < -74 + 3) %>% 
  filter(longitude > -90 & longitude < -60) %>% 
  filter(n_obs_per_person > 200) %>% 
  ggplot(aes(x = leaf_mean, y = flow_mean, color = species)) + geom_point(alpha = 0.2) +
  facet_wrap(~genus) + geom_abline(slope = 1, intercept = 0, lty = 2) + ggthemes::theme_few() +
  geom_smooth(method = lm, se = FALSE)


lf %>% 
 # filter(latitude > 40.7 - 3 & latitude < 40.7 + 3) %>% 
 # filter(longitude > -74 - 3 & longitude < -74 + 3) %>% 
  filter(longitude > -90 & longitude < -60) %>% 
  ggplot(aes(x = flow_early - leaf_early, fill = species)) + geom_histogram() +  facet_wrap(~genus, scales = "free") + theme_bw() + geom_vline(xintercept = 0, col = "red", lty = 2)


lf %>% 
#  filter(latitude > 40.7 - 3 & latitude < 40.7 + 3) %>% 
 # filter(longitude > -74 - 3 & longitude < -74 + 3) %>% 
  filter(longitude > -90) %>% 
  ggplot(aes(x = longitude, y = latitude, col = flow_early - leaf_early)) + geom_point(size = 3) + theme_bw() + facet_wrap(~genus) +
  scale_color_viridis_c()

lf %>% 
  filter(longitude > -90) %>% 
  ggplot(aes(x = tmean_spring , y = leaf_early, color = species)) + geom_point() + theme_bw() + facet_wrap(~genus) +
  scale_color_viridis_d() + geom_smooth(method = "lm") 


lf %>% 
  filter(longitude > -90) %>% 
  ggplot(aes(x = tmean_spring , y = flow_early )) + geom_point() + theme_bw() + facet_wrap(~genus, scales = "free") +
  scale_color_viridis_c() + geom_smooth(method = "lm")

lf %>% 
  filter(longitude > -90) %>% 
  filter(genus == "Quercus") %>% 
  filter(tmean_spring > 3 & tmean_spring < 9) %>% 
  ggplot(aes(x = leaf_early , y = flow_early, col = tmean_spring )) + geom_point() + theme_bw() + facet_wrap(~species, scales = "free") +
  scale_color_viridis_c() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, lty= 2)

lf %>% 
  filter(longitude > -90) %>% 
  filter(genus == "Acer") %>% 
   ggplot(aes(x = leaf_early, y = flow_early )) + geom_point() + theme_bw() + facet_wrap(~species) +
  scale_color_viridis_c() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, lty= 2)




lf %>% 
  filter(longitude > -90) %>% 
  ggplot(aes(x = longitude, y = latitude, col = tmean_spring)) + geom_point(size = 2, alpha = 0.5) + theme_bw() + 
  scale_color_viridis_c()


lf %>% 
  filter(longitude > -90) %>% 
  filter(genus == "Quercus") %>% 
  filter(tmean_spring > 3 & tmean_spring < 9) %>% 
  ggplot(aes(x = leaf_early- flow_early)) + geom_histogram() + theme_bw() + facet_wrap(~species, scales = "free") +
  scale_color_viridis_c() 



#visualize pollen release
lfp %>% 
  group_by(genus) %>%
  filter(!is.na(pol_mean)) %>% 
  summarize(n = n())
  
lfp %>% 
#  filter(tmean_spring > 3 & tmean_spring < 9) %>% 
  ggplot(aes(x = flow_mean - pol_mean)) + geom_histogram() + theme_bw() + facet_wrap(~genus, scales = "free") +
  scale_color_viridis_c() 
  
  summary(lf$n_obs_per_person)
  
fitted_models <- lf %>% 
  filter(n_obs_per_person > 1000) %>% 
  mutate(gen_sp = paste(genus, species)) %>% 
  nest(data = -gen_sp) %>% 
  mutate(mod1 = map(data, ~lm(flow_mean ~ leaf_mean, data = .)), tidied = map(mod1, tidy)) %>% unnest(tidied)

glance <- fitted_models %>% mutate(tidy = map(mod1, broom::tidy),
                                glance = map(mod1, broom::glance),
                                augment = map(mod1, broom::augment),
                                rsq = glance %>% map_dbl('r.squared'),
                                slope = tidy %>% map_dbl(function(x) x$estimate[2]))
glance %>% 
  filter(gen_sp %in% c("Quercus rubra", "Acer platanoides", "Populus deltoides"))


lf %>% 
  filter(genus == "Acer" & species == "rubrum") %>% 
  filter(longitude > -90) %>% 
  filter(greenup_0 != -9999) %>% 
  ggplot(aes(x = greenup_0, y = flow_mean)) + geom_point()

fitted_models
fitted_models$model
library(broom)
fitted_models %>% tidy(model)
rowwise(fitted_models) %>% tidy(model)
