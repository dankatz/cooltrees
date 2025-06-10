#NYC tree flowering based on start of season from planetscope and NPN relationships
# this is a companion script for NPN_flower_vs_leaf.R

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(readr)

#loading in the start of season data from Dave's PlanetScope analysis
nyc_sos <- read_csv ("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv")
#names(nyc_sos)

#get a summary of median start of spring for each species of tree for use in the NPN analysis
  nyc_sos_summary <- 
    nyc_sos %>% 
    filter(R2 > 0.8 & dbh > 3.9) %>% 
    filter(Year == 2024) %>% 
    group_by(genus, species) %>% 
    summarize(med_sos = median(SOS_50),
              p25_sos = quantile(SOS_50, 0.25),
              p75_sos = quantile(SOS_50, 0.75),
              ntrees = n())
  write_csv(nyc_sos_summary, "C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/tree_pheno_sp_summary_sos50_2024.csv")



### predicting the proportion of flowering on each particular day using NPN models ##################

## wind-pollinated species whose flowering will be predicted in nyc
  # species that are available from planetscope
  nyc_sos_species_list <- nyc_sos_summary %>% 
    filter(ntrees > 10) #exclude species with few good individuals

## species that are available from npn
#table_lf_all <- read_csv
 npn_species_available <- table_lf_all %>% filter(nobs > 10) %>% filter(!is.na(lf_dif_pred_nyc_2024_mean )) %>% 
   mutate(species = paste(genus, species),
          nobs_npn = nobs) %>%  select(species, nobs_npn) %>% arrange(species)#to match the nyc_sos format

#species represented in both nyc_sos and npn data
 nyc_npn_species_list <- left_join(nyc_sos_species_list, npn_species_available) %>% 
   filter(nobs_npn > 10)


#start species loop
for(focal_sp_i in 1:1){
  
  #inputs for loop
  focal_genus <- npn_species_to_analyze$genus[focal_sp_i] #focal_genus <-"Quercus"
  focal_species <- npn_species_to_analyze$species[focal_sp_i] #focal_species <-"rubra"

  ### load in season start date for each individual tree that we have sos estimated for
  nyc_sos_focal_sp_raw <- filter(nyc_sos, species == paste(focal_genus, focal_species))
  
  ### filtering out records where sos is questionable
  nyc_sos_focal_sp <- nyc_sos_focal_sp_raw %>% 
    filter(R2 > 0.8 & dbh > 3.9) %>% 
    filter(n_SOS_50_7day > 0) %>% #removing trees that didn't have at least one image within a week of sos50
    mutate(leaf_mean = SOS_50,  #to match the naming of the NPN model
           latitude = Lat)
  
  ### add in temperature data
    #create empty list to hold temperature data
    nyc_t_data_all_years_list <- list()
      
    yr_list <- 2017:2024 #range of NYC PlanetScope data
      for(focal_year_i in 1:length(yr_list)){ #start year loop 
        
        focal_year <- yr_list[focal_year_i] #focal_year <- 2024 
        
        
        #create a raster of mean spring temperature for a particular year
        tmean_rast_d <- pd_stack(prism_archive_subset(temp_period = "monthly", type = "tmean", mon = c(1,2,3,4), years = focal_year))
        
        #convert the points to extract to sf
        nyc_tree_yr <- nyc_sos_focal_sp %>% filter(Year == focal_year)
        nyc_tree_yr_sf <- nyc_tree_yr %>% st_as_sf(coords = c( "Lon", "Lat"), crs = 4326) 
        
        #extract temperature for all points in that year
        t_mean_nyc_months_raw <- as.data.frame(raster::extract(x = tmean_rast_d, y = nyc_tree_yr_sf)) %>% 
          rename(t_month_1 = 1, t_month_2 = 2, t_month_3 = 3, t_month_4 = 4)
        
        #add temperature data back to tree data
        t_mean_nyc_months <- bind_cols(nyc_tree_yr, t_mean_nyc_months_raw)

        #combine results into a list
        nyc_t_data_all_years_list[[focal_year_i]] <- t_mean_nyc_months
      } #end year loop
      
      #merge all the years together
      nyc_t_data_all_years <- bind_rows(nyc_t_data_all_years_list)
      #ggplot(nyc_t_data_all_years, aes(x = Year, y = t_month_4)) + geom_point() #visual check of temp data, looks fine now
    
      #join temperature back into the parent dataset
      nyc_sos_focal_sp <- left_join(nyc_sos_focal_sp, nyc_t_data_all_years)
      
   ### add in elevation in m
      dem_nyc <- terra::rast("C:/Users/dsk273/Box/Katz lab/NYC/elevation_nyc_SRTM_in_m.tif")
      #plot(dem_nyc) #downloaded from STRM manually
      
      #convert an sf version
      nyc_sos_focal_sp_sf <- nyc_sos_focal_sp %>% st_as_sf(coords = c( "Lon", "Lat"), crs = 4326) 
      
      #extract elevation and add back to parent dataset
      nyc_sos_focal_sp <- nyc_sos_focal_sp %>% 
        mutate(elevation_in_meters = terra::extract(dem_nyc, nyc_sos_focal_sp_sf)[,2])
      
    

  ### load in focal species model from NPN data
     npn_fit <- read_rds(paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/NPN_flower_leaves/", "NPN_model_",
                                  focal_genus, "_", focal_species, ".rds"))
      
   ### predicted probability of peak flowering on each individual day
      #create predictions for prediction intervals spanning from 0% to 99% of prediction interval
      nyc_sos_focal_sp_sub <- nyc_sos_focal_sp %>% #slice(nyc_sos_focal_sp, 1:200) %>% 
        #filter(Point_ID == 2542070) %>% 
        select(Point_ID, species, 
               Year, leaf_mean,
               latitude, t_month_1, t_month_2, t_month_3, t_month_4, elevation_in_meters)
      
      flow_dif_preds_nyc_list_j <- list()
      pred_levels <- seq(0, 0.99, by = 0.01)
      
      for(j in 1:100){
        flow_dif_preds_nyc_j <- predict.lm(object = fit, newdata = nyc_sos_focal_sp_sub, interval = "prediction", level = pred_levels[j],
                                           se.fit = TRUE)
        flow_dif_preds_nyc_list_j[[j]] <-
          nyc_sos_focal_sp_sub %>% 
          mutate(pred_level = pred_levels[j],
                 fit = flow_dif_preds_nyc_j$fit[,1],
                 upr = flow_dif_preds_nyc_j$fit[,3])
      }
      #combine list of prediction levels
      flow_dif_preds_nyc_v1 <- bind_rows(flow_dif_preds_nyc_list_j)
      
      #extract the proportion of the prediction distribution on each day
      flow_dif_preds_nyc_v2 <- flow_dif_preds_nyc_v1 %>% 
        mutate( fit_r = round(fit),
                difup = upr - fit, 
                difup_char = as.character(trunc(difup)),
                difup_num = as.numeric(difup_char)) %>% 
        group_by(difup_char, Point_ID, Year) %>% 
        slice_max( difup) %>% 
        group_by(Point_ID, Year) %>% 
        arrange(Point_ID, Year, pred_level) %>% 
        mutate(level_up = pred_level/2,
               prop_in_day_raw = level_up - lag(level_up, n = 1),
               prop_in_day = case_when(is.na(prop_in_day_raw) ~ level_up, .default = prop_in_day_raw)) %>% 
        ungroup()
      
      flow_dif_preds_nyc_v3 <- flow_dif_preds_nyc_v2 %>% 
        group_by(Point_ID, Year) %>% 
        mutate(difup_num = (difup_num * -1) - 1) %>% 
        #   filter(difup_num != 0) %>% 
        bind_rows(., flow_dif_preds_nyc_v2) %>% 
        arrange(Point_ID, Year, difup_num) %>% 
        mutate(dif_lf_day = round(fit) + difup_num,
               doy = leaf_mean + dif_lf_day) %>% 
        select(Point_ID, species, Year, doy, leaf_mean, dif_lf_day, prop_in_day) 
      
      # filter(flow_dif_preds_nyc_v3) %>% 
      #   ggplot(aes(x = dif_lf_day, y = prop_in_day, color = as.factor(Point_ID))) + geom_line() + geom_point()+ theme_bw() + facet_wrap(~Year)
      # 
      # test <- flow_dif_preds_nyc_v3 %>% 
      #   #filter(Point_ID == 2542070)
      #   group_by(Point_ID, Year) %>% 
      #   summarize(total_prop = sum(prop_in_day))
      # #create a long dataframe with the prediction on each day
      
      # flow_dif_preds_nyc_v3 %>% #filter(Point_ID < 200000) %>% 
      #   #ggplot(aes(x = doy, y = Point_ID)) + geom_point() + facet_wrap(~Year)
      # ggplot(aes(x = doy, y = prop_in_day)) + geom_point() + facet_wrap(~Year)
      
     flow_dif_preds_nyc_v3 %>% 
        group_by(doy, Year) %>% 
        summarize(sum_day = sum(prop_in_day)) %>% 
        ggplot(aes(x = doy, y = sum_day)) + geom_point() + facet_wrap(~Year) + theme_bw() + coord_cartesian(xlim = c(50, 200))
      
  
  ### aggregating daily flowering to pixel level
  
  
}






 
#### old stuff ########################

 
 
 
 
 
 
 
 
 
 #################################

#range of dif_leaf_flow to predict over
#NEEDS TO BE EXPANDED TO INCLUDE PREDICTIONS
pred_range_dif_leaf_flow <- round(min(ds4$dif_leaf_flow)) : round(max(ds4$dif_leaf_flow))


#create flowering predictions for all trees of a sp based on predicted dif_leaf_flow
tree_all <- list() #should probably turn this into a function instead of a loop
for(i in 1:nrow(ds5)){
  tree_i_preds <- data.frame( 
    pred_range_dif_leaf_flow = pred_range_dif_leaf_flow,
    rel_flow_day = dnorm(pred_range_dif_leaf_flow, mean = flow_dif_preds$fit[i,1], sd = flow_dif_preds$se.fit[i]))
  
  tree_i_vars <- ds5 %>%  slice(i) %>% 
    select(genus, species, individual_id, year_obs, leaf_mean) %>% 
    slice(rep(1, nrow(tree_i_preds))) 
  
  tree_i <- bind_cols(tree_i_vars, tree_i_preds)
  tree_all <- bind_rows(tree_all, tree_i)
}


#checking to see how it went across all trees
test <- tree_all %>% 
  mutate(pred_day = round(leaf_mean) + pred_range_dif_leaf_flow) %>% 
  group_by(year_obs, pred_day) %>% 
  summarize(day_total = sum(rel_flow_day))

#%>% ggplot(aes(x = pred_day, y = day_total)) + geom_point() + facet_wrap(~year_obs)




#some messing around as I think about how to effectively predict the probability of peak flowering being on each day

test <- data.frame(level = seq(0, 0.99, by = 0.01), fit = rep(NA, 100), lwr = rep(NA, 100), upr = rep(NA,100))

for(i in 1:100){
  flow_dif_preds_nyc <- predict.lm(object = fit, newdata = nyc_example_tree_pred_data, interval = "prediction", level = test$level[i],
                                   se.fit = TRUE)
  
  test$fit[i] <- flow_dif_preds_nyc$fit[,1] 
  test$lwr[i] <- flow_dif_preds_nyc$fit[,2]
  test$upr[i] <- flow_dif_preds_nyc$fit[,3]
}

#function to find the rows that correspond to the cutoff for each day
day_cutoff_fun <- function(x){
  x %>% which.
  }





test2 <- test %>% mutate( fit_r = round(fit),
                 full = upr - lwr,
                 difup = upr - fit,
                 difdown = fit - lwr) 

test2 %>% mutate(diff = abs(7 - difup)) %>%
  slice_min(diff, n = 1) %>% 
  mutate(up_level = level * 0.5)
  


# Continuous vector
v <- c(0.3, 1.6, 2.1, 2.8, 3.5, 4.9, 5.2, 5.9)

# Target integers
target_ints <- floor(1):ceiling(max(test2$difup))

# Find closest value in `v` to each integer
closest_vals <- map_dbl(target_ints, function(x) v[which.min(abs(v - x))])

# Combine results
tibble(
  target = target_ints,
  closest_value = closest_vals
)






#%>%   ggplot(aes(x = level, y = difup)) + geom_point()

plot(test2, test)
test2<-  seq(0, 0.98, by = 0.01)

test3 <- data.frame(test2, test)




test3 <- test2 %>% 
  mutate(difup_char = as.character(trunc(difup))) %>% 
  group_by(difup_char) %>% 
  slice_max( difup) %>% 
  ungroup() %>% 
  mutate(level_up = level/2,
         prop_in_day_raw = level_up - lag(level_up, n = 1),
         prop_in_day = case_when(is.na(prop_in_day_raw) ~ level_up, .default = prop_in_day_raw))

