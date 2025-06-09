#NYC tree flowering based on start of season from planetscope and NPN relationships
# this is a companion script for NPN_flower_vs_leaf.R

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(readr)


nyc_sos <- read_csv ("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv")
names(nyc_sos)

#get a summary of median start of spring for each species of tree
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





### creating predictions for the proportion of flowering on each particular day using results from prediction ##################

#list of species whose flowering will be predicted in nyc

#start species loop
for(focal_sp_i in 1:1){
  
  #inputs for loop
  focal_genus <- npn_species_to_analyze$genus[focal_sp_i] #focal_genus <-"Quercus"
  focal_species <- npn_species_to_analyze$species[focal_sp_i] #focal_species <-"rubra"

  ### load in season start date for each individual tree that we have sos estimated for
  nyc_sos_focal_sp_raw <- filter(nyc_sos, species == paste(focal_genus, focal_species))
  
  ### filtering out records where sos is questionable
  nyc_sos_focal_sp <- nyc_sos_focal_sp_raw %>% filter(R2 > 0.8 & dbh > 3.9) 
  #NEED TO ADD IN ANOTHER CHECK FOR TIME SERIES THAT HAVE HIGH R2 BUT FEW IMAGES IN SPRING
  
  ### add in temperature and other required data

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
        nyc_t_data_all_years_list[[i]] <- t_mean_nyc_months
      } #end year loop
      
      #merge all the years together
    nyc_t_data_all_years <- bind_rows(nyc_t_data_all_years_list)
      
ggplot(nyc_t_data_all_years, aes(x = Year, y = t_month_4)) + geom_point()
    
    
    
    
    #join back into the parent dataset
    indiv_leafout_t <- left_join(indiv_leafout, t_data_all_years_months)
    
  
  ### load in focal species model from NPN data

  ### create predictions of flowering at each level 
  
  ### predicted probability of peak flowering on each individual day
  
  ### creating prob peak flowering/day dataframe for each tree in long format
  
  ### aggregating daily flowering to pixel level
  
  
}













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

test %>% mutate( full = upr - lwr,
                 difup = upr - fit,
                 difdown = fit - lwr) %>% 
  ggplot(aes(x = level, y = full)) + geom_point()

plot(test2, test)
test2<-  seq(0, 0.98, by = 0.01)

test3 <- data.frame(test2, test)

