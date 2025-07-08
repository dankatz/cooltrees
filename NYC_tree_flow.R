#NYC tree flowering based on start of season from planetscope and NPN relationships
# this is a companion script for NPN_flower_vs_leaf.R

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(readr)
library(sf)
library(prism)
library(ggspatial)
library(terra)
library(tidyterra)
library(gifski)
library(stringr) 
library(purrr)
prism_set_dl_dir("C:/Users/dsk273/Documents/prism")


### PlanetScope-derived start of spring ###########################################################
# loading in the start of season data from Dave's PlanetScope analysis
nyc_sos <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv")
#names(nyc_sos)

# #stats for paper
#   length(nyc_sos$Poly_ID)
#   length(unique(nyc_sos$Poly_ID))
#   length(unique(nyc_sos$genusspeci))
#   nyc_sos %>%     
#     filter(R2 > 0.7 & dbh > 3.9) %>% 
#     filter(n_SOS_50_7day > 0) %>% 
#     filter(tpconditio == "Good" | tpconditio == "Excellent") %>% 
#     filter(tpstructur == "Full") %>% 
#     summarize(mean_r2 = mean(R2),
#               sd_r2 = sd(R2),
#               rmse_mean = mean(RMSE),
#               rmse_sd = sd(RMSE))

## QA/QC for SOS dates
nyc_sos_percentile_sp_yr <- nyc_sos %>%     
  filter(R2 > 0.7 & dbh > 3.9) %>% 
  filter(n_SOS_50_7day > 0) %>% 
  filter(tpconditio == "Good" | tpconditio == "Excellent") %>% 
  filter(tpstructur == "Full") %>% 
  group_by(genus, species, Year) %>%
  summarize(# nobs = n(),
             SOS_min = min(SOS_50),
            SOS_p2.5 = quantile(SOS_50, 0.025),
            # SOS_p25 = quantile(SOS_50, 0.25),
            # SOS_p50 = quantile(SOS_50, 0.50),
            # SOS_p75 = quantile(SOS_50, 0.75),
            SOS_p97.5 = quantile(SOS_50, 0.975),
            SOS_max = max(SOS_50)) 

nyc_sos %>% filter(species == "Acer rubrum") %>% 
  filter(tpconditio == "Good" | tpconditio == "Excellent") %>%
  filter(tpstructur == "Full") %>%
    filter(R2 > 0.7 & dbh > 3.9) %>% 
  filter(n_SOS_50_7day > 0) %>% 
  ggplot(aes(x = SOS_50)) + geom_histogram() + facet_wrap(~Year) + theme_bw()


### Table X: a summary of median start of spring  for each species of tree for use in the NPN analysis #sort(unique(nyc_sos$genus))
genera_to_include <- c("Acer", "Alnus", "Betula", "Carpinus", "Carya", "Celtis", "Corylus", "Fagus", 
                       "Fraxinus", "Ginkgo", "Gleditsia", "Juglans", "Liquidambar", "Morus", "Ostrya",
                       "Platanus", "Populus", "Quercus", "Salix", "Tilia", "Ulmus", "Zelkova")
high_nobs_spp <- c("Acer platanoides", "Acer rubrum", "Acer saccharinum", "Ginkgo biloba", "Gleditsia triacanthos",
                   "Liquidambar styraciflua", "Platanus x acerifolia", "Quercus palustris", "Quercus rubra", "Tilia americana",
                   "Tilia cordata", "Tilia tomentosa", "Ulmus americana", "Zelkova serrata")
  nyc_sos_summary <- 
    nyc_sos %>% 
    filter(R2 > 0.7 & dbh > 3.9) %>% 
    filter(n_SOS_50_7day > 0) %>% 
    filter(tpconditio == "Good" | tpconditio == "Excellent") %>% 
    filter(tpstructur == "Full") %>% 
    #filter(Year == 2024) %>% 
    group_by(genus, species, Year) %>% 
    summarize(med_sos = median(SOS_50),
              p25_sos = quantile(SOS_50, 0.25),
              p75_sos = quantile(SOS_50, 0.75),
              ntrees = n()) %>% 
    mutate(SOS_results = paste0(med_sos, " (", p25_sos, " - ", p75_sos, "): ", ntrees)) %>% 
    ungroup() %>% 
    filter(genus %in% genera_to_include) %>% 
    dplyr::select(species, SOS_results, Year) %>% 
    filter(species %in% high_nobs_spp) %>% 
    pivot_wider(names_from = Year, values_from = SOS_results)
    # filter(ntrees > 10000) %>% 
    # filter(! species %in% c("Acer", "Fraxinus", "Morus", "Quercus"))
  #write_csv(nyc_sos_summary, "C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/tree_pheno_sp_summary_sos50_2024.csv")


### Fig 2: map of Quru SOS50 in nyc, one panel per year ############################
  fig2_focal_sp_data <- 
    nyc_sos %>% 
    filter(R2 > 0.7 & dbh > 3.9) %>% 
    filter(n_SOS_50_7day > 0) %>% 
    filter(tpconditio == "Good" | tpconditio == "Excellent") %>% 
    filter(tpstructur == "Full") %>% 
    filter(species == "Quercus rubra")
  
  sos50_preds_nyc_sf <- 
    fig2_focal_sp_data %>% 
    # filter(Year == focal_year_i) %>% 
    dplyr::select(Point_ID, species, Lon, Lat, Year, SOS_50_date ) %>% 
    mutate(doy = yday(SOS_50_date),
           date_for_map = as.Date(doy, origin = "2024-01-01")) %>% 
    st_as_sf(coords = c( "Lon", "Lat"), crs = 4326)
  #hist(sos50_preds_nyc_sf$doy)
  
  #load in nyc boundary polygon
  nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
    st_union() #combine the different boroughs
  nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= crs(nyc_boundary))
  nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)
  

  # create map
  fig2_map_sos50 <- 
    ggplot() + ggthemes::theme_few() +   
    theme(panel.background = element_rect(fill='gray94', color = NA)) + 
    #geom_spatraster(data = flow_peak_rast) + 
    #scale_fill_distiller(palette = "Spectral", na.value = "gray94", name = "peak flowering (day)") +
    geom_sf(data = nyc_boundary_invert, fill = "gray94", color = "gray94") +  #
    geom_sf(data = nyc_boundary, fill = "white", color = "black") +
    geom_sf(data = sos50_preds_nyc_sf, aes(color = date_for_map), size = 1, alpha = 0.3) + 
    scale_color_gradientn(colors = c("blue", "royalblue1","plum1","mediumvioletred", "red3"),
                          limits = c(mdy("04-14-2024"), mdy("05-28-2024")),
                          name = "day of year", trans = "date",
                          na.value = NA) +
    # scale_color_viridis_c(limits = c(mdy("04-14-2024"), mdy("05-28-2024")),
    #                       name = "day of year", trans = "date",
    #                       na.value = NA, option = "turbo") +
    # scale_color_gradient2(low = "blue", mid = "plum1", high = "red2", midpoint = 19849, 
    #                       limits = c(mdy("04-14-2024"), mdy("05-31-2024")),
    #                       name = "day of year", trans = "date",
    #                       na.value = NA) +
    annotation_scale(location = "br") +  # annotation_north_arrow(location = "tl") 
    facet_wrap(~Year, nrow = 4)
  
  map_sos50_title <- paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/fig_2_Quercus_rubra_sos50.png")
  ggsave(filename = map_sos50_title, plot = fig2_map_sos50, units = "px", width = 3000, height = 3800)
  
### predicting the proportion of flowering on each particular day using NPN models ##################

## wind-pollinated species whose flowering will be predicted in nyc
  # species that are available from planetscope
  nyc_sos_species_list <- 
     nyc_sos %>% 
    filter(R2 > 0.7 & dbh > 3.9) %>% 
    filter(n_SOS_50_7day > 0) %>% 
    filter(tpconditio == "Good" | tpconditio == "Excellent") %>% 
    filter(tpstructur == "Full") %>% 
    #filter(Year == 2024) %>% 
    group_by(genus, species) %>% 
    summarize(ntrees_nyc = n()) %>% 
    filter(ntrees_nyc > 10) #exclude species with few good individuals

## species that are available from npn
table_lf_all <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/NPN_flower_leaves/table_lf_all_250708_weights_without_D.csv")
 npn_species_available <- table_lf_all %>% filter(nobs > 10) %>% filter(!is.na(lf_dif_pred_nyc_2024_mean )) %>% 
   mutate(species = paste(genus, species),
          nobs_npn = nobs) %>%  
   dplyr::select(species, nobs_npn) %>% arrange(species)#to match the nyc_sos format

#species represented in both nyc_sos and npn data
 nyc_npn_species_list <- left_join(nyc_sos_species_list, npn_species_available) %>% 
   mutate(species = gsub("^\\S+ ", "", species)) %>% 
   arrange(-ntrees_nyc) %>% 
   filter( species != "montana") %>% 
   filter(nobs_npn > 10) #print(nyc_npn_species_list, n = 30)

#species with problems: Populus tremuloides #Quercus phellos, Populus grandidentata
 
#start species loop
for(focal_sp_i in 24:26){
  
  #inputs for loop
  focal_genus <- nyc_npn_species_list$genus[focal_sp_i] #focal_genus <-"Quercus"
  focal_species <- nyc_npn_species_list$species[focal_sp_i] #focal_species <-"alba"
  print(paste(focal_genus, focal_species, "started at", Sys.time()))
  
  ### load in season start date for each individual tree that we have sos estimated for
  nyc_sos_focal_sp_raw <- filter(nyc_sos, species == paste(focal_genus, focal_species))
  
  ### filtering out records where sos is questionable
  nyc_sos_focal_sp <- nyc_sos_focal_sp_raw %>% 
    filter(R2 > 0.7 & dbh > 3.9) %>%  #dbh is in inches so 3.9 inch = 10 cm
    filter(n_SOS_50_7day > 1) %>% #removing trees that didn't have at least one image within a week of sos50
    filter(tpconditio == "Good" | tpconditio == "Excellent") %>% 
    filter(tpstructur == "Full") %>% 
    mutate(leaf_mean = SOS_50,  #to match the naming of the NPN model
           latitude = Lat) %>% 
    left_join(., nyc_sos_percentile_sp_yr) %>%  #filtering to the central 95% of the distribution to remove outliers
    filter(SOS_50 > SOS_p2.5,
           SOS_50 < SOS_p97.5)
  
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
        dplyr::select(Point_ID, species, 
               Year, leaf_mean,
               latitude, t_month_1, t_month_2, t_month_3, t_month_4, elevation_in_meters) %>% 
        rowwise() %>%  mutate(t_month_1_4 = mean(c(t_month_1, t_month_2, t_month_3, t_month_4))) %>% ungroup() #create all spring temp var
      
      
      
      flow_dif_preds_nyc_list_j <- list()
      pred_levels <- seq(0, 0.99, by = 0.01)
      
      for(j in 1:100){
        flow_dif_preds_nyc_j <- predict.lm(object = npn_fit, newdata = nyc_sos_focal_sp_sub, interval = "prediction", level = pred_levels[j],
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
        mutate(difup_num = (difup_num * -1) - 1) %>% #minus one to shift the negative dates so day 0 isn't duplicated
        bind_rows(., flow_dif_preds_nyc_v2) %>% 
        arrange(Point_ID, Year, difup_num) %>% 
        mutate(dif_lf_day = round(fit) + difup_num + 1, #shifting the whole curve back again 
               doy = leaf_mean + dif_lf_day ) %>%  
        dplyr::select(Point_ID, species, Year, doy, leaf_mean, fit, dif_lf_day, prop_in_day) %>% 
        mutate(obs_date = as.Date(doy - 1, paste0(Year, "-01-01")))
      
      #normalizing by annual totals
      flow_dif_preds_nyc_v3_annual_sum <-
        flow_dif_preds_nyc_v3 %>% 
        group_by(Year) %>% 
        summarize(sum_year = sum(prop_in_day)) 

      flow_dif_preds_nyc_v3 <-
        left_join(flow_dif_preds_nyc_v3, flow_dif_preds_nyc_v3_annual_sum) %>% 
        mutate(prop_in_day_n = prop_in_day/sum_year)
      
      
      
    ### saving the day of each percentile of peak flowering 
      #function to extract doy at a certain percentile 
      pred_level_extraction <- function(focal_percentile){
      pred_level_extraction_x <- flow_dif_preds_nyc_v3 %>% 
        group_by(Point_ID, Year) %>% 
        mutate( sum_flowpred = cumsum(prop_in_day)) %>% 
        mutate(dif_05 = abs(sum_flowpred - focal_percentile)) %>% 
        slice_min(dif_05, with_ties = FALSE) %>% 
        ungroup() %>% 
        dplyr::select(doy) %>% 
        rename_with(.cols = 1, ~paste0("focal_perc_", focal_percentile))
      return(pred_level_extraction_x)
      }
      
      #apply function to the doy when x percentage of predicted flowering interval has occurred
        pred_cols_raw <- purrr::map_dfc( list(0.05, 0.25, 0.5, 0.75, 0.95), pred_level_extraction)
      
      #add metadata to predictions      
        nyc_sos_focal_sp_sub_fpred <-
          flow_dif_preds_nyc_v3 %>% 
          dplyr::select(Point_ID, species, Year, leaf_mean, fit) %>%
          distinct() %>% 
          bind_cols(., pred_cols_raw)
      
      #save csv file for the species: summary
        preds_file_name <- paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/flowering_pred_spp/",
                                  "F_preds_nyc_", focal_genus, "_", focal_species, ".csv")
        write_csv(nyc_sos_focal_sp_sub_fpred, preds_file_name)
  
      #save csv file for the species: full length
        preds_full_file_name <- paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/flowering_pred_spp/long_version/",
                                  "F_preds_nyc_long_", focal_genus, "_", focal_species, ".csv")
        write_csv(flow_dif_preds_nyc_v3, preds_full_file_name)
        
        
### summary figure for predicted flowering distribution in each year #####################################      
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
      
  plot_interannual_timing <- 
        flow_dif_preds_nyc_v3 %>% 
        group_by(doy, Year) %>% 
        summarize(sum_day = sum(prop_in_day_n)) %>% 
        ggplot(aes(x = (doy + mdy("01-01-2024")), y = sum_day, color = as.factor(Year))) + geom_line() + theme_bw() + 
       scale_x_date(name = "date") + ylab("predicted peak flowering (proportion)") + scale_color_viridis_d(name = "Year") +
        ggtitle((paste(focal_genus, focal_species)))
    plot_interannual_timing
    ggsave(plot_interannual_timing, filename = paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/flowering_pred_maps/",
                  focal_genus,"/", focal_species, "/plot_flowering_interannual_timing.png"), create.dir = TRUE)
      
     
### maps of flowering of individual trees for each species on each day of each year ##########################
#   #load in nyc boundary polygon
#    nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>%
#                   st_union() #combine the different boroughs
#    nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= crs(nyc_boundary))
#    nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)
# 
# 
#     #creating a version of predictions that has zeros so I can sum up rasters and distinguish zeros and NAs
#       flow_dif_preds_nyc_v3_all_zeros <- flow_dif_preds_nyc_v3 %>% ungroup() %>%
#         expand(Point_ID, species, Year, doy)
#       # mutate(prop_in_day = 0,
#       #        prop_in_day_n = 0)
# 
#       flow_dif_preds_nyc_v3_with_0s <-
#         flow_dif_preds_nyc_v3 %>% ungroup() %>%
#         right_join(., flow_dif_preds_nyc_v3_all_zeros) %>%
#         arrange(Point_ID, Year, doy) %>%
#         mutate(prop_in_day = case_when(is.na(prop_in_day) ~ 0, .default = prop_in_day),
#                prop_in_day_n = case_when(is.na(prop_in_day_n) ~ 0, .default = prop_in_day_n))
# 
# #   #convert dataframe to sf
#     flow_dif_preds_nyc_v3_sf <-  nyc_sos_focal_sp %>%
#      dplyr::select(Point_ID, species, Lon, Lat) %>%
#      distinct() %>% #remove same tree from multiple years
#      left_join(flow_dif_preds_nyc_v3_with_0s, .) %>%
#      st_as_sf(coords = c( "Lon", "Lat"), crs = 4326) %>%
#        mutate(perc_in_day = prop_in_day * 100)
# 
# #create a map of peak flowering for each day
#   #create a directory to save individual frames to
#       output_directory <- paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/flowering_pred_maps/",
#                                 focal_genus,"/", focal_species, "/")
# 
#   #what days should be included?
#       doy_min_sp <- as.numeric(quantile(flow_dif_preds_nyc_v3$doy, 0.025) )
#       doy_max_sp <- as.numeric(quantile(flow_dif_preds_nyc_v3$doy, 0.975) ) #min(nyc_sos_focal_sp$SOS_p97.5)
# 
# 
# 
#   #start loop for creating frames
#   for(focal_year_i in 2017:2024){
# 
#     #create a directory if it doesn't exist yet
#     output_directory_year <- paste0(output_directory, "year_", focal_year_i)
#     if(!dir.exists(output_directory_year)){dir.create(path = output_directory_year, recursive = TRUE)}
# 
#      #start the day loop
#     for(focal_doy_i in doy_min_sp:doy_max_sp){
# 
#            flow_dif_preds_nyc_v3_sf_frame <- flow_dif_preds_nyc_v3_sf %>%
#            filter(Year == focal_year_i) %>%
#            filter(doy == focal_doy_i)
# 
#            frame_focal_date <- as.Date(focal_doy_i - 1, paste0(focal_year_i, "-01-01"))
# 
#            # aggregate daily flowering points to raster
#            flow_dif_preds_nyc_v3_terra <- terra::vect(flow_dif_preds_nyc_v3_sf_frame)
#            nyc_grid <- rast(ext(flow_dif_preds_nyc_v3_terra), resolution = 0.03, crs = crs(flow_dif_preds_nyc_v3_terra))
#            flow_rast <- terra::rasterize(flow_dif_preds_nyc_v3_terra, nyc_grid, field = "perc_in_day", fun = "mean", background = NA)
# 
#            #create map
#             map_frame <- ggplot() + ggthemes::theme_few() +
#               theme(panel.background = element_rect(fill='gray94', color = NA)) +
# 
#                 geom_spatraster(data = flow_rast) +
#                  scale_fill_gradientn(colors = c("white","khaki2","gold", "red", "red3"), limits=c(0, 10),
#                                       na.value = NA, name = "peak flowering (%)") +
# 
#                 geom_sf(data = nyc_boundary_invert, fill = "gray94", color = "gray94") +  #
#                 geom_sf(data = nyc_boundary, fill = NA, color = "black") +
# 
#                 geom_sf(data = flow_dif_preds_nyc_v3_sf_frame, aes(color = perc_in_day), size = 0.5, alpha = 0.9) +
#                     scale_color_gradientn(colors = c("white","khaki2","gold", "red", "red3"), limits=c(0, 10),
#                                           na.value = NA, name = "peak flowering (%)", guide="none") +
# 
#                  annotation_scale(location = "br") +  # annotation_north_arrow(location = "tl")
#                  ggtitle(frame_focal_date)
# 
#             map_frame_title <- paste0(output_directory_year,
#                                      "/map_frame_", focal_genus,"_", focal_species, "_", focal_year_i, "_", focal_doy_i,".png")
#             ggsave(filename = map_frame_title, plot = map_frame, units = "px", width = 3000, height = 2000)
# 
#          } ## end loop for creating map of that day
# 
#     ### create gifs from the individual frames
#     #list of png files in directory
#     gif_file_list <- list.files(path = output_directory_year, full.names = TRUE)
#     gifski(gif_file_list, gif_file = paste0(output_directory, "animation_", focal_year_i,".gif"), delay = 0.5)
# 
#     } ## end loop for year for animations
    
# ### map of peak flowering time for each species in a year #############################################
#     #start loop for each year
#     for(focal_year_i in 2017:2024){
#       
#     #load in nyc boundary polygon
#     nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
#       st_union() #combine the different boroughs
#     nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= crs(nyc_boundary))
#     nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)
#     
#     
#       flow_dif_preds_nyc_peak_date <-
#         flow_dif_preds_nyc_v3 %>% 
#         filter(Year == focal_year_i) %>% 
#         group_by(Point_ID, Year) %>% 
#         slice_max(prop_in_day_n) %>% 
#         summarize(obs_date2 = mean(obs_date),
#                   doy2 = mean(doy))
#       
#       flow_dif_preds_nyc_peak_date_sf <- 
#         nyc_sos_focal_sp %>% 
#         filter(Year == focal_year_i) %>% 
#         dplyr::select(Point_ID, species, Lon, Lat) %>% 
#         distinct() %>% #remove same tree from multiple years
#         left_join(flow_dif_preds_nyc_peak_date, .) %>% 
#         st_as_sf(coords = c( "Lon", "Lat"), crs = 4326)
#       
#       
#     # aggregate daily flowering points to raster
#       flow_dif_preds_nyc_v3_terra <- terra::vect(flow_dif_preds_nyc_peak_date_sf)
#       nyc_grid <- rast(ext(flow_dif_preds_nyc_v3_terra), resolution = 0.03, crs = crs(flow_dif_preds_nyc_v3_terra))
#       flow_peak_rast <- terra::rasterize(flow_dif_preds_nyc_v3_terra, nyc_grid, field = "doy2", fun = "mean", background = NA)
#       
#     # create map
#       map_max <- 
#         ggplot() + ggthemes::theme_few() +   
#         theme(panel.background = element_rect(fill='gray94', color = NA)) + 
#           geom_spatraster(data = flow_peak_rast) + 
#              scale_fill_distiller(palette = "Spectral", na.value = "gray94", name = "peak flowering (day)") +
#           geom_sf(data = nyc_boundary_invert, fill = "gray94", color = "gray94") +  #
#           geom_sf(data = nyc_boundary, fill = NA, color = "black") +
#         geom_sf(data = flow_dif_preds_nyc_peak_date_sf, color = "gray", size = 0.5, alpha = 0.25) + 
#           annotation_scale(location = "br")   # annotation_north_arrow(location = "tl") 
#        
#   
#       map_max_title <- paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/flowering_pred_maps/",
#                                 "/map_max_", focal_genus,"_", focal_species, "_", focal_year_i, ".png")
#       ggsave(filename = map_max_title, plot = map_max, units = "px", width = 3000, height = 2000)
#         
#     } ## end loop for year for peak flowering maps

      print(paste(focal_genus, focal_species, "finished at", Sys.time()))
      #gc()
      
} ### end species loop




### cross-species data visualization ##############################################

## load in csv files of prediction intervals for each tree that could be predicted
 nyc_fpred <-
   list.files( path = "C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/flowering_pred_spp/", full = TRUE) %>% 
    discard(., grepl("long_version", .)) %>% #remove subfolder
   map_dfr(read_csv)
 
 ### add back to the full sos data 
 nyc_sos_fpred <- nyc_sos %>% 
   filter(R2 > 0.7 & dbh > 3.9) %>%  #dbh is in inches so 3.9 inch = 10 cm
   filter(n_SOS_50_7day > 1) %>% #removing trees that didn't have at least one image within a week of sos50
   filter(tpconditio == "Good" | tpconditio == "Excellent") %>% 
   filter(tpstructur == "Full") %>% 
   dplyr::select(Point_ID, species, Year, Lon, Lat, Height, Radius, dbh, tpstructur, tpconditio, riskrating, 
          SOS_50, SOS_20, SOS_80, n_SOS_50_7day, n_SOS_50_14day,
          R2, RMSE) %>% 
      left_join(., nyc_fpred)
   # #filtering to the central 95% of the distribution to remove outliers
   # filter(SOS_50 > SOS_p2.5,
   #        SOS_50 < SOS_p97.5)
   
 

### table 3: flower prediction interval in nyc ####################################
 table3 <- nyc_sos_fpred %>% 
   group_by(species, Year) %>% 
   summarize(nobs = n(),
             flowmax_mean = round(mean(focal_perc_0.5, na.rm = TRUE), 0),
             flowmax_25 = round(mean(focal_perc_0.25 , na.rm = TRUE), 0),
             flowmax_75 = round(mean(focal_perc_0.75 , na.rm = TRUE), 0),
             leafmax_mean = round(mean(leaf_mean, na.rm = TRUE), 0),
             leafmax_mean = round(mean(leaf_mean, na.rm = TRUE), 0)
             ) %>% 
   pivot_wider(id_cols = species, names_from = c(Year), values_from = c(flowmax_mean, flowmax_25, flowmax_75, leafmax_mean, nobs)) %>% 
   rowwise() %>% 
   mutate(y2017 = paste0(flowmax_mean_2017, " (", flowmax_25_2017, " - ", flowmax_75_2017, "): ", nobs_2017),
          y2018 = paste0(flowmax_mean_2018, " (", flowmax_25_2018, " - ", flowmax_75_2018, "): ", nobs_2018),
          y2019 = paste0(flowmax_mean_2019, " (", flowmax_25_2019, " - ", flowmax_75_2019, "): ", nobs_2019),
          y2020 = paste0(flowmax_mean_2020, " (", flowmax_25_2020, " - ", flowmax_75_2020, "): ", nobs_2020),
          y2021 = paste0(flowmax_mean_2021, " (", flowmax_25_2021, " - ", flowmax_75_2021, "): ", nobs_2021),
          y2022 = paste0(flowmax_mean_2022, " (", flowmax_25_2022, " - ", flowmax_75_2022, "): ", nobs_2022),
          y2023 = paste0(flowmax_mean_2023, " (", flowmax_25_2023, " - ", flowmax_75_2023, "): ", nobs_2023),
          y2024 = paste0(flowmax_mean_2024, " (", flowmax_25_2024, " - ", flowmax_75_2024, "): ", nobs_2024)
          ) %>% 
   filter(!is.na(flowmax_mean_2023) ) %>% 
   dplyr::select(species, contains("y2"))
 
 
 

### Fig 3: map of Quru Fmax in nyc, one panel per year ############################
 fig3_focal_sp_data <- nyc_sos_fpred %>% 
   filter(species == "Quercus palustris")
   
   #load in nyc boundary polygon
   nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
      st_union() #combine the different boroughs
   nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= crs(nyc_boundary))
   nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)
 
 flow_preds_nyc_peak_date_sf <- 
   fig3_focal_sp_data %>% 
     # filter(Year == focal_year_i) %>% 
   dplyr::select(Point_ID, species, Lon, Lat, Year, focal_perc_0.5) %>% 
   mutate(date_for_map = as.Date(focal_perc_0.5, origin = "2024-01-01")) %>% 
   st_as_sf(coords = c( "Lon", "Lat"), crs = 4326)
 
 ## anomaly within each year
 flow_preds_nyc_year_means <- 
   flow_preds_nyc_peak_date_sf %>% 
   group_by(species, Year) %>% 
   summarize(mean_fmax_year = mean(focal_perc_0.5, na.rm = TRUE)) %>% 
   st_drop_geometry()
 
 flow_preds_nyc_peak_date_sf <- left_join(flow_preds_nyc_peak_date_sf, flow_preds_nyc_year_means) %>% 
   mutate(dif_within_year = focal_perc_0.5 - mean_fmax_year )
 
 # flow_dif_preds_nyc_peak_date <-
 #   flow_dif_preds_nyc_v3 %>% 
 #   filter(Year == focal_year_i) %>% 
 #   group_by(Point_ID, Year) %>% 
 #   slice_max(prop_in_day_n) %>% 
 #   summarize(obs_date2 = mean(obs_date),
 #             doy2 = mean(doy))
 # aggregate daily flowering points to raster
 # flow_dif_preds_nyc_v3_terra <- terra::vect(flow_dif_preds_nyc_peak_date_sf)
 # nyc_grid <- rast(ext(flow_dif_preds_nyc_v3_terra), resolution = 0.03, crs = crs(flow_dif_preds_nyc_v3_terra))
 # flow_peak_rast <- terra::rasterize(flow_dif_preds_nyc_v3_terra, nyc_grid, field = "doy2", fun = "mean", background = NA)
 
 # create map
 fig3_map_flowmax <- 
   ggplot() + ggthemes::theme_few() +   
   theme(panel.background = element_rect(fill='gray94', color = NA)) + 
   #geom_spatraster(data = flow_peak_rast) + 
   #scale_fill_distiller(palette = "Spectral", na.value = "gray94", name = "peak flowering (day)") +
   geom_sf(data = nyc_boundary_invert, fill = "gray94", color = "gray94") +  #
   geom_sf(data = nyc_boundary, fill = "white", color = "black") +
   geom_sf(data = flow_preds_nyc_peak_date_sf, aes(color = date_for_map), size = 1, alpha = 0.5) + 
   scale_color_gradient2(low = "springgreen4", mid = "yellow2", high = "orangered", midpoint = 19845, name = "day of year", trans = "date",
                         na.value = NA) +
   annotation_scale(location = "br") +  # annotation_north_arrow(location = "tl") 
   facet_wrap(~Year, nrow = 4)
     
 map_max_title <- paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/fig_3_Quercus_rubra_flow_max.png")
 ggsave(filename = map_max_title, plot = fig3_map_flowmax, units = "px", width = 4000, height = 4500)
 
 
 # #version with anomalies per year
 # fig3_map_flowmax_anomo <- 
 #   ggplot() + ggthemes::theme_few() +   
 #   theme(panel.background = element_rect(fill='gray94', color = NA)) + 
 #   #geom_spatraster(data = flow_peak_rast) + 
 #   #scale_fill_distiller(palette = "Spectral", na.value = "gray94", name = "peak flowering (day)") +
 #   geom_sf(data = nyc_boundary_invert, fill = "gray94", color = "gray94") +  #
 #   geom_sf(data = nyc_boundary, fill = "white", color = "black") +
 #   geom_sf(data = flow_preds_nyc_peak_date_sf, aes(color = dif_within_year), size = 1, alpha = 0.5) + 
 #   scale_color_gradient2(low = "springgreen4", mid = "yellow2", high = "orangered", midpoint = 0, name = "difference (days)", 
 #                         na.value = NA) +
 #   annotation_scale(location = "br") +  # annotation_north_arrow(location = "tl") 
 #   facet_wrap(~Year, nrow = 4)
 # 
 # map_max_title <- paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/fig_3_Quercus_rubra_flow_max_anomol.png")
 # ggsave(filename = map_max_title, plot = fig3_map_flowmax_anomo, units = "px", width = 3000, height = 5000)
 

 ### Fig 4: map of Quercus species Fmax in nyc, one panel per year ############################
 focal_year_i <- 2023
 
 #load in nyc boundary polygon
 nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
   st_union() #combine the different boroughs
 nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= crs(nyc_boundary))
 nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)

 oak_spp <- c("Quercus alba", "Quercus falcata", "Quercus macrocarpa", "Quercus palustris", "Quercus rubra",  "Quercus velutina")
 flow_preds_nyc_peak_date_sf <- 
   nyc_sos_fpred %>% 
   filter(species %in% oak_spp) %>%  #make sure to leave a space
   #filter(species %in% str_subset(unique(nyc_sos_fpred$species), "Quercus ")) %>%  #make sure to leave a space
   filter(Year == focal_year_i) %>% 
   dplyr::select(Point_ID, species, Lon, Lat, Year, focal_perc_0.5) %>% 
   mutate(date_for_map = as.Date(focal_perc_0.5, origin = "2024-01-01")) %>% 
   st_as_sf(coords = c( "Lon", "Lat"), crs = 4326)
 

 # create map
 fig4_map_flowmax <- 
   ggplot() + ggthemes::theme_few() +   
   theme(panel.background = element_rect(fill='gray94', color = NA), strip.text = element_text(face = "italic")) + 
   #geom_spatraster(data = flow_peak_rast) + 
   #scale_fill_distiller(palette = "Spectral", na.value = "gray94", name = "peak flowering (day)") +
   geom_sf(data = nyc_boundary_invert, fill = "gray94", color = "gray94") +  #
   geom_sf(data = nyc_boundary, fill = "white", color = "black") +
   geom_sf(data = flow_preds_nyc_peak_date_sf, aes(color = date_for_map), size = 1.3, alpha = 0.9) + 
   scale_color_gradient2(low = "springgreen4", mid = "yellow2", high = "orangered",
                         midpoint = 19840, name = "day of year", trans = "date", na.value = NA) +
   # scale_color_gradient2(low = "blue", mid = "gold", high = "red", midpoint = 19840, name = "day of year", trans = "date",
   #                       na.value = NA) +
   annotation_scale(location = "br") +  # annotation_north_arrow(location = "tl") 
   facet_wrap(~species, nrow = 3)
 
 map_max_title <- paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/fig_4_Quercus_spp_flow_max.png")
 ggsave(filename = map_max_title, plot = fig4_map_flowmax, units = "px", width = 4000, height = 4600)
 
 
 
 ### Fig x: map of various Acer species Fmax in nyc, shown in one year ############################
 focal_year_i <- 2022
 
 #load in nyc boundary polygon
 nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
   st_union() #combine the different boroughs
 nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= crs(nyc_boundary))
 nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)
 
 acer_spp <- c("Acer rubrum", "Acer platanoides", "Acer saccharum", "Acer saccharinum", "Acer negundo")
 flow_preds_nyc_peak_date_sf <- 
   nyc_sos_fpred %>% 
   filter(species %in% acer_spp) %>%  #make sure to leave a space
   #filter(species %in% str_subset(unique(nyc_sos_fpred$species), "Acer ")) %>%  #make sure to leave a space
   filter(Year == focal_year_i) %>% 
   dplyr::select(Point_ID, species, Lon, Lat, Year, focal_perc_0.5) %>% 
   mutate(date_for_map = as.Date(focal_perc_0.5, origin = "2024-01-01")) %>% 
   st_as_sf(coords = c( "Lon", "Lat"), crs = 4326)
 
 
 # create map
 figx_map_flowmax <- 
   ggplot() + ggthemes::theme_few() +   
   theme(panel.background = element_rect(fill='gray94', color = NA), strip.text = element_text(face = "italic")) + 
   #geom_spatraster(data = flow_peak_rast) + 
   #scale_fill_distiller(palette = "Spectral", na.value = "gray94", name = "peak flowering (day)") +
   geom_sf(data = nyc_boundary_invert, fill = "gray94", color = "gray94") +  #
   geom_sf(data = nyc_boundary, fill = "white", color = "black") +
   geom_sf(data = flow_preds_nyc_peak_date_sf, aes(color = date_for_map), size = 1.2, alpha = 0.9) + 
   scale_color_gradient2(low = "springgreen4", mid = "yellow2", high = "orangered",
                         midpoint = 19820, name = "day of year", trans = "date", na.value = NA) +
   # scale_color_gradient2(low = "blue", mid = "gold", high = "red", midpoint = 19840, name = "day of year", trans = "date",
   #                       na.value = NA) +
   annotation_scale(location = "br") +  # annotation_north_arrow(location = "tl") 
   facet_wrap(~species, nrow = 3)
 
 map_max_title <- paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/fig_x_Acer_spp_flow_max.png")
 ggsave(filename = map_max_title, plot = figx_map_flowmax, units = "px", width = 4000, height = 4600)
 

# # ### Fig X: Distribution of flowering probability for each oak species in each year #######################################
#  oak_spp <- c("Quercus alba", "Quercus macrocarpa", "Quercus palustris", "Quercus rubra",  "Quercus velutina")
# 
#   nyc_sos_fpred_quercus <-
#    nyc_sos_fpred %>%
#    filter(species %in% oak_spp) %>%
#    group_by(species, focal_perc_0.5, Year) %>%
#     summarize(sum_day = n()) %>%
#    arrange(species, Year, focal_perc_0.5) %>%
#    group_by(species, Year) %>%
#    mutate(cumulative_number_peaked = cumsum(sum_day))
# 
#    nyc_sos_fpred_quercus_cumsum <-
#      nyc_sos_fpred_quercus %>%
#      slice_max(cumulative_number_peaked) %>%
#      mutate(cum_total = cumulative_number_peaked) %>%
#      dplyr::select(species, Year, cum_total)
# 
#    nyc_sos_fpred_quercus <- left_join( nyc_sos_fpred_quercus,  nyc_sos_fpred_quercus_cumsum) %>%
#      mutate(cum_prop = cumulative_number_peaked/cum_total)
# 
#    ggplot(nyc_sos_fpred_quercus, aes(x = (focal_perc_0.5 + mdy("01-01-2024")), y = cum_prop, color = species)) +
#      geom_line() + geom_point(size = 0.6)+ ggthemes::theme_few() +
#      scale_x_date(name = "date", limits = c( mdy("04-01-2024"),  mdy("06-01-2024"))) + ylab("predicted peak flowering (proportion)") +
#      scale_color_viridis_d(name = "Year") + facet_wrap(~Year)

### Fig X: comparing distribution of flowering times for some species ####################################
 oak_spp <- c("Quercus alba", "Quercus macrocarpa", "Quercus palustris", "Quercus rubra",  "Quercus velutina")
 maple_spp <- c("Acer rubrum", "Acer saccharum", "Acer saccharinum", "Acer negundo", "Acer platanoides")
 birch_spp <- c("Betula lenta", "Betula nigra", "Betula papyrifera")
 hickory_spp <- c("Carya glabra", "Carya ovata")
 
 
 #load in the long format predictions (amount per individual per day)
   nyc_fpred_long <-
     list.files( path = "C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/flowering_pred_spp/long_version/", full = TRUE) %>% 
     map_dfr(read_csv)
   
   ### add back to the full sos data 
   nyc_sos_fpred_long <- nyc_sos %>% 
     filter(R2 > 0.7 & dbh > 3.9) %>%  #dbh is in inches so 3.9 inch = 10 cm
     filter(n_SOS_50_7day > 0) %>% #removing trees that didn't have at least one image within a week of sos50
     filter(tpconditio == "Good" | tpconditio == "Excellent") %>% 
     filter(tpstructur == "Full") %>% 
     dplyr::select(Point_ID, species, Year, Lon, Lat, Height, Radius, dbh, tpstructur, tpconditio, riskrating, 
                   SOS_50, SOS_20, SOS_80, n_SOS_50_7day, n_SOS_50_14day,
                   R2, RMSE) %>% 
     left_join(nyc_fpred_long, .)
   # #filtering to the central 95% of the distribution to remove outliers
   # filter(SOS_50 > SOS_p2.5,
   #        SOS_50 < SOS_p97.5)
   
   
   fig_x_plot_interannual_timing <- 
     nyc_sos_fpred_long %>% 
     filter(species %in% maple_spp) %>% 
     #filter(species != "Quercus phellos" & species != "Quercus falcata") %>% 
     group_by(species, doy, Year) %>% 
     summarize(sum_day = sum(prop_in_day_n)) %>% 
     ggplot(aes(x = (doy + mdy("01-01-2024")), y = sum_day, color = as.factor(species))) +  
     geom_vline(aes(xintercept = mdy("04-01-2024")), color = "gray") + 
     geom_line(lwd = 1.1) + 
     ggthemes::theme_few() + 
     scale_x_date(name = "date", limits = c(mdy("03-01-2024"), mdy("05-15-2024"))) + ylab("predicted flowering interval (proportion)") + 
     scale_color_viridis_d(name = "species") +
     facet_wrap(~Year, ncol = 2) + guides(color = guide_legend(theme = theme(legend.text = element_text(face = "italic"))))
   
   fig_x_plot_interannual_timing
   # ggsave(plot_interannual_timing, filename = paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/flowering_pred_maps/",
   #                                                   focal_genus,"/", focal_species, "/plot_flowering_interannual_timing.png"), create.dir = TRUE)
   
   
   
 #### SI section: comparison to airborne pollen concentrations ###############################################
 #Lincoln Center from Guy Robinson
 #https://docs.google.com/spreadsheets/d/1m2Op3h86NFPQUpPDugqACCJatEhR_-ea_C_0nLNdSOE/edit#gid=0

 ## processing Guy's dataset
     #manually remove the first blank 2 rows and initial couple of columns that have an unusual format
     pol_raw <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/airborne pollen/Lincoln_Center_semi_processed_google_sheets_download_250620.csv", 
                         name_repair = "minimal") %>% 
       t() %>% 
       as_tibble() %>% 
       `colnames<-`(c("date_obs", "p_tot", "nothing", "perc_by_taxa")) %>% 
       dplyr::select(-nothing) 
    
     pol_entry_list <- list()
     
     for(i in 1:17){ #most number of taxa recorded on a single day: 17
     pol_entry_i <- 
       pol_raw %>% 
       mutate(p_total = parse_number(p_tot),
              entry1 = str_split_i(string = perc_by_taxa, pattern = "%", i = i),
              taxon = tolower(str_trim(str_replace_all(entry1, "[^a-zA-Z ]", ""))),
              pol_perc = parse_number(entry1),
              pol_perc = case_when(entry1 == "No Pollen" ~ 0, .default = pol_perc),
              pol_pcm = p_total * (pol_perc/100) ) %>% 
       dplyr::select(date_obs, taxon, pol_pcm) %>% 
       filter(!is.na(taxon)) %>% 
       filter(taxon != "")
     pol_entry_list[[i]] <- pol_entry_i
     }
     
     pol <- bind_rows(pol_entry_list) %>% 
       arrange(date_obs, taxon) %>% 
       mutate(date_obs = mdy(date_obs),
              year_obs = year(date_obs),
              doy = yday(date_obs))
     
     #adding in rows for missing dates
     date_obs <- data.frame(date_obs = seq(ymd('2022-01-01'),ymd('2025-06-15'), by = '1 day'))
     pol_NAs <-  left_join(date_obs, pol) 

    #amount of pollen by each taxon
     pol_NAs %>% group_by(taxon) %>% summarize(total = sum(pol_pcm, na.rm = TRUE)) %>% arrange(-total) %>% print(n = 20)
     
 #create figure for a genus for each year #  sort(unique(nyc_sos_fpred$species))
     oak_spp <- c("Quercus rubra", "Quercus palustris", "Quercus alba", "Quercus macrocarpa", "Quercus velutina") #"Quercus falcata", "Quercus phellos"
     maple_spp <- c("Acer","Acer rubrum", "Acer saccharum", "Acer saccharinum", "Acer negundo", "Acer platanoides")
     birch_spp <- c("Betula","Betula lenta", "Betula nigra", "Betula papyrifera")
     
     
     dist_param <- 10
     common_name <- "maple"
     species_list <- maple_spp
       
       #which individiuals are within a certain distance of Lincoln Center?
       LC <- st_sfc(st_point(c(-73.9834, 40.7724)), crs = 4326) #lincolnc center
       nyc_sos_fpred_LC <- 
         nyc_sos_fpred %>% 
         filter(species %in% species_list) %>% 
         st_as_sf(coords = c( "Lon", "Lat"), crs = 4326) %>% 
         mutate(LC_dis = as.numeric(st_distance(., LC)/1000)) %>% 
         st_drop_geometry() %>% 
         filter(n_SOS_50_7day > 0) %>% 
         filter(tpconditio == "Good" | tpconditio == "Excellent") %>% 
         filter(tpstructur == "Full") 
         #ggplot(nyc_sos_fpred_LC) +  geom_sf(aes(color = LC_dis))
       
      #50 percent peak 
        nyc_sos_fpred_quercus_50 <- 
         nyc_sos_fpred_LC %>% 
         filter(species %in% species_list) %>% 
         filter(LC_dis < dist_param) %>% 
         group_by(Year, focal_perc_0.5) %>% 
         summarize(total_per_day_at_peak = n()) %>% 
         mutate(cumulative_number_peaked = cumsum(total_per_day_at_peak)) 
       nyc_sos_fpred_quercus_50_cumsum <- nyc_sos_fpred_quercus_50 %>% slice_max(cumulative_number_peaked) %>% 
         mutate(cum_total = cumulative_number_peaked) %>% 
         dplyr::select(Year, cum_total)
       nyc_sos_fpred_quercus_50 <- left_join( nyc_sos_fpred_quercus_50,  nyc_sos_fpred_quercus_50_cumsum) %>% 
         mutate(cum_prop = cumulative_number_peaked/cum_total) %>% 
         rename(doy = focal_perc_0.5) %>% 
         mutate(year_obs = Year) %>% 
         filter(year_obs > 2021)
       
      #25 percent peak 
       nyc_sos_fpred_quercus_25 <- 
         nyc_sos_fpred_LC %>% 
         filter(species %in% species_list) %>% 
         filter(LC_dis < dist_param) %>% 
         group_by(Year, focal_perc_0.25) %>% 
         summarize(total_per_day_at_peak = n()) %>% 
         mutate(cumulative_number_peaked = cumsum(total_per_day_at_peak)) 
       nyc_sos_fpred_quercus_25_cumsum <- nyc_sos_fpred_quercus_25 %>% slice_max(cumulative_number_peaked) %>% 
         mutate(cum_total = cumulative_number_peaked) %>% 
         dplyr::select(Year, cum_total)
       nyc_sos_fpred_quercus_25 <- left_join( nyc_sos_fpred_quercus_25,  nyc_sos_fpred_quercus_25_cumsum) %>% 
         mutate(cum_prop = cumulative_number_peaked/cum_total) %>% 
         rename(doy = focal_perc_0.25) %>% 
         mutate(year_obs = Year) %>% 
         filter(year_obs > 2021)
       
      #75 percent peak 
       nyc_sos_fpred_quercus_75 <- 
         nyc_sos_fpred_LC %>% 
         filter(species %in% species_list) %>% 
         filter(LC_dis < dist_param) %>% 
         group_by(Year, focal_perc_0.75) %>% 
         summarize(total_per_day_at_peak = n()) %>% 
         mutate(cumulative_number_peaked = cumsum(total_per_day_at_peak)) 
       nyc_sos_fpred_quercus_75_cumsum <- nyc_sos_fpred_quercus_75 %>% slice_max(cumulative_number_peaked) %>% 
         mutate(cum_total = cumulative_number_peaked) %>% 
         dplyr::select(Year, cum_total)
       nyc_sos_fpred_quercus_75 <- left_join( nyc_sos_fpred_quercus_75,  nyc_sos_fpred_quercus_75_cumsum) %>% 
         mutate(cum_prop = cumulative_number_peaked/cum_total) %>% 
         rename(doy = focal_perc_0.75) %>% 
         mutate(year_obs = Year) %>% 
         filter(year_obs > 2021)
       
       #05 percent peak 
       nyc_sos_fpred_quercus_05 <- 
         nyc_sos_fpred_LC %>% 
         filter(species %in% species_list) %>% 
         filter(LC_dis < dist_param) %>% 
         group_by(Year, focal_perc_0.05) %>% 
         summarize(total_per_day_at_peak = n()) %>% 
         mutate(cumulative_number_peaked = cumsum(total_per_day_at_peak)) 
       nyc_sos_fpred_quercus_05_cumsum <- nyc_sos_fpred_quercus_05 %>% slice_max(cumulative_number_peaked) %>% 
         mutate(cum_total = cumulative_number_peaked) %>% 
         dplyr::select(Year, cum_total)
       nyc_sos_fpred_quercus_05 <- left_join( nyc_sos_fpred_quercus_05,  nyc_sos_fpred_quercus_05_cumsum) %>% 
         mutate(cum_prop = cumulative_number_peaked/cum_total) %>% 
         rename(doy = focal_perc_0.05) %>% 
         mutate(year_obs = Year) %>% 
         filter(year_obs > 2021)
       
       #95 percent peak 
       nyc_sos_fpred_quercus_95 <- 
         nyc_sos_fpred_LC %>% 
         filter(species %in% species_list) %>% 
         filter(LC_dis < dist_param) %>% 
         group_by(Year, focal_perc_0.95) %>% 
         summarize(total_per_day_at_peak = n()) %>% 
         mutate(cumulative_number_peaked = cumsum(total_per_day_at_peak)) 
       nyc_sos_fpred_quercus_95_cumsum <- nyc_sos_fpred_quercus_95 %>% slice_max(cumulative_number_peaked) %>% 
         mutate(cum_total = cumulative_number_peaked) %>% 
         dplyr::select(Year, cum_total)
       nyc_sos_fpred_quercus_95 <- left_join( nyc_sos_fpred_quercus_95,  nyc_sos_fpred_quercus_95_cumsum) %>% 
         mutate(cum_prop = cumulative_number_peaked/cum_total) %>% 
         rename(doy = focal_perc_0.95) %>% 
         mutate(year_obs = Year) %>% 
         filter(year_obs > 2021)
       
       
       pol_NAs_qusp_max <- pol_NAs %>% 
         filter(taxon == common_name | is.na(taxon)) %>% 
         group_by(year_obs) %>% 
         summarize(pol_yr_pcm_max = max(pol_pcm, na.rm = TRUE))
        
       pol_NAs %>% 
         filter(taxon == common_name | is.na(taxon)) %>% 
         left_join(., pol_NAs_qusp_max) %>% 
         mutate(pol_pcm_prop = pol_pcm / pol_yr_pcm_max) %>% 
          filter(year_obs > 2019 & year_obs < 2025) %>% 
         filter(doy < 200) %>% 
         filter(doy >75) %>% 
       ggplot() + 
         geom_point(aes(x = doy, y = pol_pcm_prop)) + 
         geom_line(aes(x = doy, y = pol_pcm_prop), alpha = 0.5) +
          facet_wrap(~year_obs, scales = "free_y", ncol = 1) +
         geom_line(data = nyc_sos_fpred_quercus_50, aes(x = doy, y = cum_prop ), color = "red")+
         geom_line(data = nyc_sos_fpred_quercus_25, aes(x = doy, y = cum_prop ), color = "orange")+
         geom_line(data = nyc_sos_fpred_quercus_75, aes(x = doy, y = cum_prop ), color = "orange")+
         geom_line(data = nyc_sos_fpred_quercus_05, aes(x = doy, y = cum_prop ), color = "yellow")+
         geom_line(data = nyc_sos_fpred_quercus_95, aes(x = doy, y = cum_prop ), color = "yellow")+
         scale_x_continuous(limits = c(100, 150))+
         ylab("airborne pollen (relative to annual maximum)") +
         scale_y_continuous(sec.axis = sec_axis(~., name="predicted flowering (cumulative proportion of trees per year)") ) +
         ggthemes::theme_few()  +
         theme(
           axis.line.y.right = element_line(color = "red"),  # Change secondary axis line color
           axis.text.y.right = element_text(color = "red"),   # Change secondary axis text color
           axis.title.y.right = element_text(color = "red")
         ) 

 
       
       
       
       
# ### create figure for Quercus for each year ######
#        #days of peak flowering
#        #oak_spp <- c("Quercus rubra", "Quercus pallustris", "Quercus alba", "Quercus macrocarpa", "Quercus velutina")
#       
#       #50 percent peak 
#         nyc_sos_fpred_quercus_50 <- 
#          nyc_sos_fpred %>% 
#          filter(species %in% oak_spp) %>% 
#          group_by(Year, focal_perc_0.5) %>% 
#          summarize(total_per_day_at_peak = n()) %>% 
#          mutate(cumulative_number_peaked = cumsum(total_per_day_at_peak)) 
#        nyc_sos_fpred_quercus_50_cumsum <- nyc_sos_fpred_quercus_50 %>% slice_max(cumulative_number_peaked) %>% 
#          mutate(cum_total = cumulative_number_peaked) %>% 
#          dplyr::select(Year, cum_total)
#        nyc_sos_fpred_quercus_50 <- left_join( nyc_sos_fpred_quercus_50,  nyc_sos_fpred_quercus_50_cumsum) %>% 
#          mutate(cum_prop = cumulative_number_peaked/cum_total) %>% 
#          rename(doy = focal_perc_0.5) %>% 
#          mutate(year_obs = Year) %>% 
#          filter(year_obs > 2021)
#        
#       #25 percent peak 
#        nyc_sos_fpred_quercus_25 <- 
#          nyc_sos_fpred %>% 
#          filter(species %in% oak_spp) %>% 
#          group_by(Year, focal_perc_0.25) %>% 
#          summarize(total_per_day_at_peak = n()) %>% 
#          mutate(cumulative_number_peaked = cumsum(total_per_day_at_peak)) 
#        nyc_sos_fpred_quercus_25_cumsum <- nyc_sos_fpred_quercus_25 %>% slice_max(cumulative_number_peaked) %>% 
#          mutate(cum_total = cumulative_number_peaked) %>% 
#          dplyr::select(Year, cum_total)
#        nyc_sos_fpred_quercus_25 <- left_join( nyc_sos_fpred_quercus_25,  nyc_sos_fpred_quercus_cumsum) %>% 
#          mutate(cum_prop = cumulative_number_peaked/cum_total) %>% 
#          rename(doy = focal_perc_0.25) %>% 
#          mutate(year_obs = Year) %>% 
#          filter(year_obs > 2021)
#        
#       #75 percent peak 
#        nyc_sos_fpred_quercus_75 <- 
#          nyc_sos_fpred %>% 
#          filter(species %in% oak_spp) %>% 
#          group_by(Year, focal_perc_0.75) %>% 
#          summarize(total_per_day_at_peak = n()) %>% 
#          mutate(cumulative_number_peaked = cumsum(total_per_day_at_peak)) 
#        nyc_sos_fpred_quercus_75_cumsum <- nyc_sos_fpred_quercus_75 %>% slice_max(cumulative_number_peaked) %>% 
#          mutate(cum_total = cumulative_number_peaked) %>% 
#          dplyr::select(Year, cum_total)
#        nyc_sos_fpred_quercus_75 <- left_join( nyc_sos_fpred_quercus_75,  nyc_sos_fpred_quercus_cumsum) %>% 
#          mutate(cum_prop = cumulative_number_peaked/cum_total) %>% 
#          rename(doy = focal_perc_0.75) %>% 
#          mutate(year_obs = Year) %>% 
#          filter(year_obs > 2021)
#        
#        #05 percent peak 
#        nyc_sos_fpred_quercus_05 <- 
#          nyc_sos_fpred %>% 
#          filter(species %in% oak_spp) %>% 
#          group_by(Year, focal_perc_0.05) %>% 
#          summarize(total_per_day_at_peak = n()) %>% 
#          mutate(cumulative_number_peaked = cumsum(total_per_day_at_peak)) 
#        nyc_sos_fpred_quercus_05_cumsum <- nyc_sos_fpred_quercus_05 %>% slice_max(cumulative_number_peaked) %>% 
#          mutate(cum_total = cumulative_number_peaked) %>% 
#          dplyr::select(Year, cum_total)
#        nyc_sos_fpred_quercus_05 <- left_join( nyc_sos_fpred_quercus_05,  nyc_sos_fpred_quercus_cumsum) %>% 
#          mutate(cum_prop = cumulative_number_peaked/cum_total) %>% 
#          rename(doy = focal_perc_0.05) %>% 
#          mutate(year_obs = Year) %>% 
#          filter(year_obs > 2021)
#        
#        #95 percent peak 
#        nyc_sos_fpred_quercus_95 <- 
#          nyc_sos_fpred %>% 
#          filter(species %in% oak_spp) %>% 
#          group_by(Year, focal_perc_0.95) %>% 
#          summarize(total_per_day_at_peak = n()) %>% 
#          mutate(cumulative_number_peaked = cumsum(total_per_day_at_peak)) 
#        nyc_sos_fpred_quercus_95_cumsum <- nyc_sos_fpred_quercus_95 %>% slice_max(cumulative_number_peaked) %>% 
#          mutate(cum_total = cumulative_number_peaked) %>% 
#          dplyr::select(Year, cum_total)
#        nyc_sos_fpred_quercus_95 <- left_join( nyc_sos_fpred_quercus_95,  nyc_sos_fpred_quercus_cumsum) %>% 
#          mutate(cum_prop = cumulative_number_peaked/cum_total) %>% 
#          rename(doy = focal_perc_0.95) %>% 
#          mutate(year_obs = Year) %>% 
#          filter(year_obs > 2021)
#        
#        
#        pol_NAs_qusp_max <- pol_NAs %>% 
#          filter(taxon == "oak" | is.na(taxon)) %>% 
#          group_by(year_obs) %>% 
#          summarize(pol_yr_pcm_max = max(pol_pcm, na.rm = TRUE))
#         
#        pol_NAs %>% 
#          filter(taxon == "oak" | is.na(taxon)) %>% 
#          left_join(., pol_NAs_qusp_max) %>% 
#          mutate(pol_pcm_prop = pol_pcm / pol_yr_pcm_max) %>% 
#           filter(year_obs > 2019 & year_obs != 2025) %>% 
#          filter(doy < 200) %>% 
#          filter(doy >75) %>% 
#        ggplot() + 
#          geom_point(aes(x = doy, y = pol_pcm_prop)) + 
#          geom_line(aes(x = doy, y = pol_pcm_prop), alpha = 0.5) +
#           facet_wrap(~year_obs, scales = "free_y", ncol = 1) +
#          geom_line(data = nyc_sos_fpred_quercus_50, aes(x = doy, y = cum_prop ), color = "red")+
#          geom_line(data = nyc_sos_fpred_quercus_25, aes(x = doy, y = cum_prop ), color = "orange")+
#          geom_line(data = nyc_sos_fpred_quercus_75, aes(x = doy, y = cum_prop ), color = "orange")+
#          geom_line(data = nyc_sos_fpred_quercus_05, aes(x = doy, y = cum_prop ), color = "yellow")+
#          geom_line(data = nyc_sos_fpred_quercus_95, aes(x = doy, y = cum_prop ), color = "yellow")+
#          scale_x_continuous(limits = c(100, 150))+
#          ylab("airborne pollen (relative to annual maximum)") +
#          scale_y_continuous(sec.axis = sec_axis(~., name="predicted flowering (cumulative proportion of trees per year)") ) +
#          ggthemes::theme_few()  +
#          theme(
#            axis.line.y.right = element_line(color = "red"),  # Change secondary axis line color
#            axis.text.y.right = element_text(color = "red"),   # Change secondary axis text color
#            axis.title.y.right = element_text(color = "red")
#          ) 
# 
#  #length(unique(nyc_slength(unique(nyc_slength(unique(nyc_sos$Point_ID))
# 
#        
       
### figure SI X: comparing flowering times for different oak species vs airborne pollen #################
       dist_param <- 5
       
       LC <- st_sfc(st_point(c(-73.9834, 40.7724)), crs = 4326) #lincolnc center
       nyc_sos_fpred_LC <- 
         nyc_sos_fpred %>% 
         filter(species %in% oak_spp) %>% 
         st_as_sf(coords = c( "Lon", "Lat"), crs = 4326) %>% 
         mutate(LC_dis = as.numeric(st_distance(., LC)/1000)) %>% 
         st_drop_geometry() %>% 
         filter(LC_dis < dist_param)
       
       
  #50 percent peak 
  nyc_sos_fpred_quercus_50 <- 
    nyc_sos_fpred_LC %>% 
    group_by(species, Year, focal_perc_0.5) %>% 
    summarize(total_per_day_at_peak = n()) %>% 
    mutate(cumulative_number_peaked = cumsum(total_per_day_at_peak)) 
  nyc_sos_fpred_quercus_50_cumsum <- nyc_sos_fpred_quercus_50 %>% slice_max(cumulative_number_peaked) %>% 
    mutate(cum_total = cumulative_number_peaked) %>% 
    dplyr::select(Year, cum_total)
  nyc_sos_fpred_quercus_50 <- left_join( nyc_sos_fpred_quercus_50,  nyc_sos_fpred_quercus_50_cumsum) %>% 
    mutate(cum_prop = cumulative_number_peaked/cum_total) %>% 
    rename(doy = focal_perc_0.5) %>% 
    mutate(year_obs = Year) %>% 
    filter(year_obs > 2021)
  
  pol_NAs_qusp_max <- pol_NAs %>% 
    filter(taxon == "oak" | is.na(taxon)) %>% 
    group_by(year_obs) %>% 
    summarize(pol_yr_pcm_max = max(pol_pcm, na.rm = TRUE))
  
  pol_NAs %>% 
    filter(taxon == "oak" | is.na(taxon)) %>% 
    left_join(., pol_NAs_qusp_max) %>% 
    mutate(pol_pcm_prop = pol_pcm / pol_yr_pcm_max) %>% 
    filter(year_obs > 2019 & year_obs != 2025) %>% 
    filter(doy < 200) %>% 
    filter(doy >75) %>% 
    ggplot() + 
    geom_point(aes(x = doy, y = pol_pcm_prop)) + 
    geom_line(aes(x = doy, y = pol_pcm_prop), alpha = 0.5) +
    facet_wrap(~year_obs, scales = "free_y", ncol = 1) +
    geom_line(data = nyc_sos_fpred_quercus_50, aes(x = doy, y = cum_prop, color = species)) +
    # geom_line(data = nyc_sos_fpred_quercus_25, aes(x = doy, y = cum_prop, color = species), lty = 2) +
    # geom_line(data = nyc_sos_fpred_quercus_75, aes(x = doy, y = cum_prop, color = species), lty = 2) +
    # geom_line(data = nyc_sos_fpred_quercus_05, aes(x = doy, y = cum_prop, color = species), lty = 3) +
    # geom_line(data = nyc_sos_fpred_quercus_95, aes(x = doy, y = cum_prop, color = species), lty = 4) +
    scale_x_continuous(limits = c(100, 150))+
    ylab("airborne pollen (relative to annual maximum)") +
    scale_y_continuous(sec.axis = sec_axis(~., name="predicted flowering (cumulative proportion of trees per year)") ) +
    ggthemes::theme_few()  +
    theme(
      axis.line.y.right = element_line(color = "red"),  # Change secondary axis line color
      axis.text.y.right = element_text(color = "red"),   # Change secondary axis text color
      axis.title.y.right = element_text(color = "red")
    ) 
       

  
  
### figure SI X: comparing flowering times for different oak species vs airborne pollen #################
  focal_genus_common <- "oak"
  #focal_genus_spp <- str_subset(unique(nyc_sos_fpred$species), "Ginkgo ") #make sure to leave a space
  focal_genus_spp <- oak_spp  
  dist_param <- 100
  
  LC <- st_sfc(st_point(c(-73.9834, 40.7724)), crs = 4326) #lincolnc center
  nyc_sos_fpred_LC <- 
    nyc_sos_fpred %>% 
    filter(species %in% focal_genus_spp) %>% 
    st_as_sf(coords = c( "Lon", "Lat"), crs = 4326) %>% 
    mutate(LC_dis = as.numeric(st_distance(., LC)/1000)) %>% 
    st_drop_geometry()
  
  #50 percent peak 
  nyc_sos_fpred_quercus_50 <- 
    nyc_sos_fpred_LC %>% 
    filter(species %in% focal_genus_spp) %>% 
    filter(LC_dis < dist_param) %>% 
    filter(n_SOS_50_7day > 0) %>% 
    filter(tpconditio == "Good" | tpconditio == "Excellent") %>% 
    filter(tpstructur == "Full") %>% 
    group_by(species, Year, focal_perc_0.5) %>% 
    summarize(total_per_day_at_peak = n()) %>% 
    mutate(cumulative_number_peaked = cumsum(total_per_day_at_peak)) 
  nyc_sos_fpred_quercus_50_cumsum <- nyc_sos_fpred_quercus_50 %>% slice_max(cumulative_number_peaked) %>% 
    mutate(cum_total = cumulative_number_peaked) %>% 
    dplyr::select(Year, cum_total)
  nyc_sos_fpred_quercus_50 <- left_join( nyc_sos_fpred_quercus_50,  nyc_sos_fpred_quercus_50_cumsum) %>% 
    mutate(cum_prop = cumulative_number_peaked/cum_total) %>% 
    rename(doy = focal_perc_0.5) %>% 
    mutate(year_obs = Year) %>% 
    filter(year_obs > 2016)
  
  pol_NAs_qusp_max <- pol_NAs %>% 
    filter(taxon == focal_genus_common | is.na(taxon)) %>% 
    group_by(year_obs) %>% 
    summarize(pol_yr_pcm_max = max(pol_pcm, na.rm = TRUE))
  
  pol_NAs %>% 
    filter(taxon == focal_genus_common | is.na(taxon)) %>% 
    left_join(., pol_NAs_qusp_max) %>% 
    mutate(pol_pcm_prop = pol_pcm / pol_yr_pcm_max) %>% 
    filter(year_obs > 2016 ) %>% 
    filter(doy < 200) %>% 
    filter(doy >75) %>% 
    ggplot() + 
    geom_point(aes(x = doy, y = pol_pcm_prop)) + 
    geom_line(aes(x = doy, y = pol_pcm_prop), alpha = 0.5) +
    facet_wrap(~year_obs, scales = "free_y", ncol = 1) +
    geom_line(data = nyc_sos_fpred_quercus_50, aes(x = doy, y = cum_prop, color = species)) +
    # geom_line(data = nyc_sos_fpred_quercus_25, aes(x = doy, y = cum_prop, color = species), lty = 2) +
    # geom_line(data = nyc_sos_fpred_quercus_75, aes(x = doy, y = cum_prop, color = species), lty = 2) +
    # geom_line(data = nyc_sos_fpred_quercus_05, aes(x = doy, y = cum_prop, color = species), lty = 3) +
    # geom_line(data = nyc_sos_fpred_quercus_95, aes(x = doy, y = cum_prop, color = species), lty = 4) +
    scale_x_continuous(limits = c(100, 150))+
    ylab("airborne pollen (relative to annual maximum)") +
    scale_y_continuous(sec.axis = sec_axis(~., name="predicted flowering (cumulative proportion of trees per year)") ) +
    ggthemes::theme_few()  +
    theme(
      axis.line.y.right = element_line(color = "darkred"),  # Change secondary axis line color
      axis.text.y.right = element_text(color = "darkred"),   # Change secondary axis text color
      axis.title.y.right = element_text(color = "darkred")
    ) 
  