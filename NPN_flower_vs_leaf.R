# This script is for calculating the relationships between leaf out and flowering from NPN data
# and from additional data collected by DK 

# set up work environment
  library(rnpn)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(readr)
  library(prism)
  library(sf)
  library(purrr)
  library(broom)
  library(robustbase)


### select top anemophilous taxa from NPN #####################################
  #anemophilous angiosperms of potential allergenic concern that are present in the Eastern US
  #npn_species() %>%  filter(functional_type == "Deciduous broadleaf") %>% arrange(genus) -> test
  acer_species_list <- c(1843,59,778,2,1591,60,779,780,3,781,61,1199)
  alnus_species_list <- c(62,63,319)
  betula_species_list <- c(97, 1439, 98, 1850, 1339, 1851, 99, 1805)
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
  quercus_species_list <- c(705,100,1365,757,1870,987,1690,1484,988,316,297,1485,1190,1486,
                            301,101,1691,1212,989,1366,102,1756,1213,1755,1487,1159,305)
  salix_species_list <- c(1007, 717, 1875, 1494, 20166, 293, 322, 1493, 1006, 1163, 77, 1371, 1009, 2278, 1008, 1372, 1010, 1876)
  tilia_species_list <- c(93, 1775, 1776, 1777)
  ulmus_species_list <- c(1192,1048,1049,1215,1216)
  zelkova_species_list <- c(1374)
  
  list_all_focal_taxa <- c(acer_species_list, alnus_species_list, betula_species_list, carya_species_list, celtis_species_list, corylus_species_list,
                           fraxinus_species_list, ginkgo_species_list, gleditsia_species_list, juglans_species_list, morus_species_list,
                           platanus_species_list, populus_species_list, quercus_species_list, salix_species_list, tilia_species_list,
                           ulmus_species_list, zelkova_species_list
                           )

### download and process NPN data #####################################
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
     additional_fields = c("Observed_Status_Conflict_Flag", "partner_group", "observedby_person_id",  "Observation_Comments",#
                           "AGDD", "gdd", "Daylength") #"Greenup_0", "MidGreenup_0", "Maturity_0") #
    #additional_fields documentation : https://docs.google.com/document/d/1yNjupricKOAXn6tY1sI7-EwkcfwdGUZ7lxYv7fcPjO8/edit?tab=t.0
   # six_leaf_layer = TRUE, six_bloom_layer = TRUE # agdd_layer = 32 
  )
  
  #get missing elevation data and add it back in
    elevation_data <- filter(npn_direct, elevation_in_meters == -9999) %>% 
    dplyr::select(longitude, latitude, observation_id) %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  elev_epqs <- elevatr::get_elev_point(elevation_data, prj = 4326, src = "epqs")
  elev_epqs_nosf <- sf::st_drop_geometry(elev_epqs) %>% dplyr::select(-elev_units)
  
  npn_direct <- left_join(npn_direct, elev_epqs_nosf) %>% 
    mutate(elevation_in_meters  = case_when(elevation_in_meters == -9999 ~ elevation, .default = elevation_in_meters ))
  #write_csv(npn_direct, "C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/npn_flow_leaf_download250701.csv")
  #npn_direct <- read_csv("C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/npn_flow_leaf_download250701.csv")
  

  
  ## clean up and organize data
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

    #length(unique(npn_direct$species_id))
    # assess an individual observer
    # test <- npn_direct %>% 
    #   filter(observedby_person_id == "2002") %>% 
    #   filter(individual_id == 4318) %>% 
    #   filter(phenophase_id == 467) %>% 
    #   ggplot(aes(x = ymd(observation_date), y = intensity_p, color = phenophase_description)) + geom_point()
  
  ### calculate the average date of leaf expansion phenophase for each individual tree in each year
  indiv_leafout <- npn_direct %>% filter(phenophase_id == 467) %>% 
    filter(phenophase_status == 1) %>% 
    filter(day_of_year < 173) %>%  #solstice
    group_by(genus, species, latitude, longitude, individual_id, year_obs) %>% 
    slice(which.min(abs(intensity_p - 0.5))) %>% #select the observation with the highest intensity
    filter(intensity_p > 0.24 & intensity_p < 0.75) # select observations with high proportion of leaves that are expanding
  
  indiv_leafout_summary <- npn_direct %>% filter(phenophase_id == 467) %>% 
    filter(phenophase_status == 1) %>% 
    filter(day_of_year < 173) %>%  #solstice
    filter(intensity_p > 0.24 & intensity_p < 0.75) %>% 
    group_by(genus, species, individual_id, year_obs) %>% 
    summarize(leaf_mean = mean(day_of_year),
              leaf_early = min(day_of_year),
              leaf_late = max(day_of_year),
              leaf_nobs = n()) %>% 
    mutate(leafout_duration = leaf_late - leaf_early)
    
  indiv_leafout <- left_join(indiv_leafout, indiv_leafout_summary) %>% 
    filter(leafout_duration < 15)

  ### calculate the average date of open flowers for each individual in each year
    indiv_flow <- npn_direct %>% filter(phenophase_id == 501) %>% 
      filter(phenophase_status == 1) %>% 
      filter(day_of_year < 173) %>%  #solstice
      group_by(genus, species, latitude, longitude, individual_id, year_obs) %>% 
      slice(which.min(abs(intensity_p - 1))) %>% 
      filter(intensity_p > 0.24) #"what percentage of all fresh flowers on the plant are open?"
    
    indiv_flow_summary <- npn_direct %>% filter(phenophase_id == 501) %>% 
      filter(phenophase_status == 1) %>% 
      filter(day_of_year < 173) %>%  #solstice
      group_by(genus, species, individual_id, year_obs) %>% 
      filter(intensity_p > 0.24) %>% 
      summarize(flow_mean = mean(day_of_year),
                flow_early = min(day_of_year),
                flow_late = max(day_of_year),
                flow_nobs = n()) %>% 
      mutate(flow_duration = flow_late - flow_early)
    
    indiv_flow <- left_join(indiv_flow, indiv_flow_summary) %>% 
      filter(flow_duration < 15)
    
    indiv_flow_join <- indiv_flow %>% 
      dplyr::select(year_obs, individual_id, latitude, longitude, genus, species, common_name, flow_mean, flow_early, flow_late, flow_nobs, flow_duration, flow_intensity_p = intensity_p)

 ## join leaf expansion and flowering and retain observations that have both
  lf <- left_join(indiv_leafout, indiv_flow_join) %>% 
        filter(!is.na(flow_mean)) 
  # length(unique(lf$observedby_person_id))
  # length(unique(lf$individual_id))
  
  ### calculate the average date of pollen release for each individual in each year
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

### process non-NPN data from DK and other sources #############################
   # #data from Detroit in 2017
   #  #this is from the script: "C:\Users\dsk273\Box\MIpostdoc\trees\Phenology and daily variation in pollen release\parsing_pheno_measurements.R"
   #  #documentation is available: "phenology and relative release SOP 170626.docx" in the same folder
   #   phenotree <- readr::read_csv("C:/Users/dsk273/Box/MIpostdoc/trees/Phenology and daily variation in pollen release/phenotree250418.csv") %>%
   #     st_as_sf(coords = c( "POINT_X", "POINT_Y"), crs = 3857) %>% #turns out the field collected data were in WGS84 Web Mercator
   #     st_transform(., 4326) %>%
   #     mutate(longitude = sf::st_coordinates(.)[,1],
   #            latitude = sf::st_coordinates(.)[,2]) %>%
   #     mutate(year_obs = 2017,
   #            genus = stringr::word(Species, 1),
   #            species = stringr::word(Species, 2),
   #            species = case_when(genus == "Platanus" & species == "x" ~ "acerifolia", .default = species),
   #            elevation_in_meters = 200) %>% # actual elevation of sites was all within ~25 m of this
   #    st_drop_geometry() %>%
   #    filter(!(pheno_ID == "76d7b4e0-e525-439c-a03d-9be53b3968b3" & created_date == "4/21/2017 19:37")) %>%
   #    filter(!(pheno_ID == "8040fab7-2224-4c6e-8b47-356df9b13605" & created_date == "4/28/2017 18:01"))
   # 
   #   #phenotree_sf$lon; phenotree_sf$lat; plot(phenotree_sf[,1])  #checks
   #      #test <- phenotree %>% filter(Species == "Platanus x acerfolia") %>% arrange( tree, date)
   # 
   #  # #add in additional observations of Platanus acerifolia from DK from several years based on photos
   #  #   phenotree_additional_obs <- read_csv("C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/manual Plac observations from DK 250613.csv") %>%
   #  #     mutate( date = mdy(date),
   #  #             julia = yday(date),
   #  #             year_obs = year(date))
   #    #
   #    # phenotree <- bind_rows(phenotree, phenotree_additional_obs)
   # 
   #   ## extract monthly spring temperature in each year, one year at a time
   #       t_data_all_years_list <- list()
   #       t_data_all_years_months <- phenotree
   #       yr_list <- unique(phenotree$year_obs)
   #       for(j in 1:4){ #start month loop
   #         for(i in 1:length(yr_list)){ #start year loop
   #           #create a raster of mean spring temperature for a particular year
   #           tmean_rast_d <- prism_archive_subset(temp_period = "monthly", type = "tmean", mon = j, years = yr_list[i])
   #           tmean_rast2_d <- pd_stack(tmean_rast_d)
   #           tmean_focal_year <- raster::calc(tmean_rast2_d, mean) #raster::plot(r_mean)
   # 
   #           #convert the points to extract to sf
   #           l_yr <- phenotree %>% filter(year_obs == yr_list[i])
   #           l_yr_sf <- phenotree %>% st_as_sf(coords = c( "longitude", "latitude"), crs = 4326)
   # 
   #           #extract temperature for all points in that year
   #            t_data_focal_year <- mutate(l_yr, "{paste0('t_month_', j)}" := unlist( #dynamically renaming variable
   #             raster::extract(x = tmean_focal_year, y = l_yr_sf)))
   # 
   #           #combine results into a list
   #           t_data_all_years_list[[i]] <- t_data_focal_year
   #         } #end year loop
   # 
   #         #merge all the years together
   #         t_data_all_years_month <- bind_rows(t_data_all_years_list)
   # 
   #         #join the different months
   #         t_data_all_years_months <- left_join(t_data_all_years_months, t_data_all_years_month)
   #       } #end month loop
   # 
   # 
   #   #join back into the parent dataset
   #       phenotree <- left_join(phenotree, t_data_all_years_months) #names(phenotree)
   # 
   #   #calculate the day of peak leaf expansion ("Stage D") for each tree
   #     #"Stage D: leaves partly expanded" vs. "Stage C: new shoots emerging" vs. "Stage E: leaves fully expanded"
   #    d_leafmax <- phenotree %>%
   #      group_by(genus, species, longitude, latitude, elevation_in_meters, year_obs, tree, t_month_1, t_month_2, t_month_3, t_month_4) %>%
   #      slice_max(leaf_d) %>%
   #      summarize(date_mean = mean(date), #taking the middle date if there are multiple
   #                leaf_mean = mean(julia),
   #                intensity = mean(leaf_d)/100)%>%
   #      ungroup() %>%
   #      filter(intensity > 0.25) #remove the one tree where peak leaf expansion wasn't observed
   # 
   #  # calculate the day of peak mature flowers ("stage Y") for each tree
   #  # "Stage y: flowers fully developed" vs. "Stage x: immature inflorescence visible" vs. "Stage z: flowers senescent"
   #    d_flowmax <- phenotree %>%
   #      group_by(genus, species, Species, longitude, latitude, elevation_in_meters, year_obs, tree, t_month_1, t_month_2, t_month_3, t_month_4) %>%
   #      slice_max(flower_y) %>%
   #      summarize(date_mean = mean(date),
   #                flow_mean = mean(julia),
   #                intensity = mean(flower_y)/100) %>%
   #      filter(!(Species == "Platanus x acerfolia" & intensity < 0.01)) %>%
   #      filter(!(Species != "Platanus x acerfolia" & intensity < 0.25))%>%
   #      ungroup()
   #  
   #    #write_csv(d_leafmax, "C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/detroit_leafmax_250630.csv")
   #    #write_csv(d_flowmax, "C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/detroit_flowmax_250630.csv")
   #    #d_leafmax <- read_csv("C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/detroit_leafmax_250630.csv")
   #    #d_flowmax <- read_csv("C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/detroit_flowmax_250630.csv")
      
### extract temperature data for each location in each year ####################
  prism_set_dl_dir("C:/Users/dsk273/Documents/prism")
  # get_prism_monthlys(years = 2025:2025, type = "tmean") #updating my local dataset
  # prism_archive_clean(type = "tmean", temp_period = "monthly", years = 2012:2025) #needed to clean out some provisional data
  # get_prism_dailys(minDate = ymd("2025-05-01"), maxDate = ymd("2025-06-05"), type = "tmean") #updating my local dataset


  ## prepare data to extract
    lf_extract <- lf %>% ungroup() %>% 
      mutate(year_obs = year(observation_date),
             too_recent = case_when(year_obs == 2025 & flow_mean > 150 ~ "too recent", #currently set to include through may 2025
                                    year_obs == 2025 & leaf_mean > 150 ~ "too recent",
                                    TRUE ~ "not too recent")) %>% 
      filter(too_recent == "not too recent") %>% #removing observations from this month
      dplyr::select(individual_id, longitude, latitude, year_obs, flow_mean, leaf_mean) %>% 
      distinct() 


    ### function to extract temperature for days around flowering/leaf out
        # extract_temp_daily <- function(individual_id, year_obs, flow_mean, leaf_mean, longitude, latitude){ 
        #   individual_id_focal <- individual_id
        #   longitude_focal <- longitude
        #   latitude_focal <- latitude
        #   yr_focal <- year_obs
        #   #year_obs <- lf_sf$year_obs[1]; flow_mean <- lf_sf$flow_mean; leaf_mean <- lf_sf$leaf_mean[1]; individual_id_focal <- lf_sf$individual_id[1]
        #   date_min_doy <- max(flow_mean) - 30 #can substitute in flowering or leaf or the max of both here
        #   date_max_doy <- max(flow_mean)
        #   date_min <- ymd(paste0(yr_focal, "/1/1")) + date_min_doy
        #   date_max <- ymd(paste0(yr_focal, "/1/1")) + date_max_doy
        #   
        #   #print(date_min); print(date_max)
        #   tmean_rast_d <- prism_archive_subset(temp_period = "daily", type = "tmean", minDate = date_min, maxDate = date_max)
        #   tmean_rast2_d <- pd_stack(tmean_rast_d)
        #   r_mean <- raster::calc(tmean_rast2_d, mean) #raster::plot(r_mean)
        #   
        #   lf_sub <- filter(lf_extract, year_obs == yr_focal, individual_id == individual_id_focal) 
        #   #print(lf_sub)
        #   lf_sub_sf <- lf_sub %>% st_as_sf(coords = c( "longitude", "latitude"), crs = 4326) 
        #   
        #   tmean_data <- unlist(raster::extract(x = r_mean,
        #                                        y = lf_sub_sf)) %>% as.data.frame()
        # 
        #   lf_sub2 <- lf_sub %>% mutate(tmean = as.numeric(unlist(tmean_data)))
        # 
        #   print(tmean_data)
        #   return(lf_sub2)
        # }
        # 
        # #apply function to calculate temp for each observation individually
        # #lf_test <- lf_extract[20:21,]
        # spring_temp <- purrr::pmap_dfr(.l = lf_extract, .f = extract_temp_daily)
        # 
        # npn_active_flow <- left_join(npn_active_flow, spring_temp)



  ### extract monthly spring temperature in each year, one year at a time
    tmean_data_all_years_list <- list()
    yr_list <- 2012:2025
    for(i in 1:length(yr_list)){
      #create a raster of mean spring temperature for a particular year
        tmean_rast_d <- prism_archive_subset(temp_period = "monthly", type = "tmean", mon = c(1, 2, 3, 4), years = yr_list[i])
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
    
    
    #write_csv(lfp, "C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/npn_lfp_250709.csv")
    #lfp <- read_csv("C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/npn_lfp_250702.csv")



### get temperature for indiv_leafout for finding outliers ########################################################
    #starting with the raw data and extracting relevant bits
    indiv_leafout_peak <- npn_direct %>% filter(phenophase_id == 467) %>% 
      filter(phenophase_status == 1) %>% 
      group_by(genus, species, individual_id, year_obs) %>% 
      slice(which.min(abs(intensity_p - 0.5))) %>% #select the observations that are closest to mid-leaf expansion
      filter(intensity_p > 0.24 & intensity_p < 0.76)
    
    indiv_leafout_summary <- npn_direct %>% filter(phenophase_id == 467) %>% 
      filter(phenophase_status == 1) %>% 
      filter(intensity_p > 0.24 & intensity_p < 0.76) %>% 
      group_by(genus, species, individual_id, year_obs) %>% 
      summarize(leaf_mean = mean(day_of_year),
                leaf_early = min(day_of_year),
                leaf_late = max(day_of_year),
                leaf_nobs = n()) %>% 
      mutate(leafout_duration = leaf_late - leaf_early) %>% 
      ungroup()
    
    #filter out bad observations and prepare data to extract
    indiv_leafout <- left_join(indiv_leafout_peak, indiv_leafout_summary) %>% 
      filter(n_obs_per_person > 10) %>%  #remove observers that took fewer than 10 observations
      filter(leafout_duration < 15) %>% #remove observations that had leaf out duration of > 2 weeks 
      filter(leaf_late < 173) %>%  #observations where leaf out occurred after the solstice
       mutate(year_obs = year(observation_date),
             too_recent = case_when(year_obs == 2025 & leaf_mean > 170 ~ "too recent",
                                    TRUE ~ "not too recent")) %>% 
      filter(too_recent == "not too recent") %>% #removing observations from this month (NPN has more recent data than PRISM)
      dplyr::select(genus, species, 
                    individual_id, n_obs_per_person, longitude, latitude, year_obs, 
                    leaf_early, leaf_mean, leaf_late, leaf_nobs, leafout_duration, intensity_p #gdd, daylength
                    ) %>% 
      ungroup()
    
    ## extract monthly spring temperature in each year, one year at a time
    t_data_all_years_list <- list()
    t_data_all_years_months <- indiv_leafout
    yr_list <- 2012:2025
    for(j in 1:4){ #start month loop
        for(i in 1:length(yr_list)){ #start year loop
          #create a raster of mean spring temperature for a particular year
          tmean_rast_d <- prism_archive_subset(temp_period = "monthly", type = "tmean", mon = j, years = yr_list[i])
          tmean_rast2_d <- pd_stack(tmean_rast_d)
          tmean_focal_year <- raster::calc(tmean_rast2_d, mean) #raster::plot(r_mean)
          
          #convert the points to extract to sf
          l_yr <- indiv_leafout %>% filter(year_obs == yr_list[i])
          l_yr_sf <- l_yr %>% st_as_sf(coords = c( "longitude", "latitude"), crs = 4326) 
          
          #extract temperature for all points in that year
          t_data_focal_year <- mutate(l_yr, "{paste0('t_month_', j)}" := unlist( #dynamically renaming variable
                                      raster::extract(x = tmean_focal_year, y = l_yr_sf)))
          
          #combine results into a list
          t_data_all_years_list[[i]] <- t_data_focal_year
        }
    
      #merge all the years together
      t_data_all_years_month <- bind_rows(t_data_all_years_list)
      
      #join the different months
      t_data_all_years_months <- left_join(t_data_all_years_months, t_data_all_years_month)
    } #end month loop
    
    
    #join back into the parent dataset and create an average spring temperature variables
    indiv_leafout_t <- left_join(indiv_leafout, t_data_all_years_months) %>% 
    rowwise() %>%  mutate(t_month_1_4 = mean(c(t_month_1, t_month_2, t_month_3, t_month_4))) %>% ungroup() #create all spring temp var
    
    #write_csv(indiv_leafout_t, "C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/indiv_leafout_t_250709.csv")
    #indiv_leafout_t <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/indiv_leafout_t_250709.csv")
       

### get temperature for indiv_flow for finding outliers ########################################################
    #starting with the raw data and extracting relevant bits
      indiv_flow_peak <- npn_direct %>% filter(phenophase_id == 501) %>% 
        filter(phenophase_status == 1) %>% 
        group_by(genus, species, individual_id, year_obs) %>% 
        slice(which.min(abs(intensity_p - 1))) %>% #select the observations that are closest to mid-leaf expansion
        filter(intensity_p > 0.24 ) #"what percentage of all fresh flowers on the plant are open?"
      
      indiv_flow_summary <- npn_direct %>% filter(phenophase_id == 501) %>% 
        filter(phenophase_status == 1) %>% 
        filter(intensity_p > 0.24 ) %>% 
        group_by(genus, species, individual_id, year_obs) %>% 
        summarize(flow_mean = mean(day_of_year),
                  flow_early = min(day_of_year),
                  flow_late = max(day_of_year),
                  flow_nobs = n()) %>% 
        mutate(flow_duration = flow_late - flow_early) %>% 
        ungroup()
    
    #filter out bad observations and prepare data to extract
      indiv_flow <- left_join(indiv_flow_peak, indiv_flow_summary) %>% 
        filter(n_obs_per_person > 10) %>% 
        filter(flow_duration < 15) %>% 
        filter(flow_late < 173) %>%  #solstice
        mutate(year_obs = year(observation_date),
               too_recent = case_when(year_obs == 2025 & flow_mean > 170 ~ "too recent",
                                      TRUE ~ "not too recent")) %>% 
        filter(too_recent == "not too recent") %>% #removing observations from this month
        dplyr::select(genus, species, 
                      individual_id, n_obs_per_person, longitude, latitude, year_obs, 
                      flow_early, flow_mean, flow_late, flow_nobs, flow_duration, intensity_p,
                      gdd, daylength) %>% 
        ungroup()
      
    ## extract monthly spring temperature in each year, one year at a time
    t_data_all_years_list <- list()
    t_data_all_years_months <- indiv_flow
    yr_list <- 2012:2025
    for(j in 1:4){ #start month loop
            for(i in 1:length(yr_list)){ #start year loop
              #create a raster of mean spring temperature for a particular year
              tmean_rast_d <- prism_archive_subset(temp_period = "monthly", type = "tmean", mon = j, years = yr_list[i])
              tmean_rast2_d <- pd_stack(tmean_rast_d)
              tmean_focal_year <- raster::calc(tmean_rast2_d, mean) #raster::plot(r_mean)
              
              #convert the points to extract to sf
              l_yr <- indiv_flow %>% filter(year_obs == yr_list[i])
              l_yr_sf <- l_yr %>% st_as_sf(coords = c( "longitude", "latitude"), crs = 4326) 
              
              #extract temperature for all points in that year
              t_data_focal_year <- mutate(l_yr, "{paste0('t_month_', j)}" := unlist( #dynamically renaming variable
                raster::extract(x = tmean_focal_year, y = l_yr_sf)))
              
              #combine results into a list
              t_data_all_years_list[[i]] <- t_data_focal_year
            } #end year loop
        
        #merge all the years together
        t_data_all_years_month <- bind_rows(t_data_all_years_list)
        
        #join the different months
        t_data_all_years_months <- left_join(t_data_all_years_months, t_data_all_years_month)
    } #end month loop
    
    
    #join back into the parent dataset
     indiv_flow_t <- left_join(indiv_flow, t_data_all_years_months) %>% 
       rowwise() %>%  mutate(t_month_1_4 = mean(c(t_month_1, t_month_2, t_month_3, t_month_4))) %>% ungroup() #create all spring temp var
     
     #write_csv(indiv_flow_t, "C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/indiv_flow_t_250709.csv")
     #indiv_flow_t <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/indiv_flow_t_250709.csv")
 
     
### analyze the difference between leaf out and flowering for each species and create tables and figures in a loop ################
  #load in required datasets for this portion of the script
     lfp <- read_csv("C:/Users/dsk273/Box/NYC projects/NYC flowering and leaf phenology spring 25/npn_lfp_250709.csv")
     indiv_flow_t <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/indiv_flow_t_250709.csv")
     indiv_leafout_t <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/indiv_leafout_t_250709.csv")
     
     
    ## which species to include from npn
    npn_species_to_analyze <- lfp %>% group_by(genus, species) %>% summarize(n = n()) %>% arrange(-n) %>% 
    filter(!(genus == "Quercus" & species == "douglasii"), #not an east coast species
           !(genus == "Populus" & species == "fremontii"),  #not an east coast species
           !(genus == "Celtis" & species == "occidentalis"),  #too few observations
           !(genus == "Fraxinus" & species == "pennsylvanica"),  #only 6 observations with both leaf and flowering
           !(genus == "Acer" & species == "spicatum"),  #only 4 observations with both leaf and flowering
           !(genus == "Corylus" & species == "americana"),  #only 5 observations with both leaf and flowering
           !(genus == "Gleditsia" & species == "triacanthos"),  #only 10 observations with both leaf and flowering
           !(genus == "Quercus" & species == "stellata"),  #too few observations
           !(genus == "Salix" & species == "bebbiana") #too few observations
            ) #%>% print(n = 40) #sum(npn_species_to_analyze$n)
    
  # #loading in the nyc data to get median leaf out date per species for nyc in 2024 for table 1
  #   #this file was generated in the script "NYC_tree_flow.R"
  #   nyc_sos_summary <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/tree_pheno_sp_summary_sos50_2024.csv")


## start species loop
    #creating some empty lists to paste into
    table_lf_all_list <- list()
    table_si_all_list <- list()
    
for(focal_sp_i in 28:28){ # print(npn_species_to_analyze, n = 40)

  #inputs for loop 
    focal_genus <- npn_species_to_analyze$genus[focal_sp_i] #focal_genus <-"Quercus"
    focal_species <- npn_species_to_analyze$species[focal_sp_i] #focal_species <-"rubra"
    print(paste(focal_genus, focal_species, focal_sp_i))


  ### find leaf outliers and extract their weights
    ds_leaf <- indiv_leafout_t %>% # version without Detroit data 
    #ds_leaf <- bind_rows(indiv_leafout_t, d_leafmax) %>%   #add Detroit data from DK indiv_leafout_t %>% 
    filter(genus == focal_genus, species == focal_species) %>% filter(!is.na(t_month_1))  #ggplot(ds, aes(x = gdd, y = leaf_mean)) + geom_point() + theme_bw()
   
  #fit a robust regression using the MM estimator
    #lmrob.control(setting="KS2014") 
    mmfit_leaf <- lmrob(leaf_mean ~  t_month_1_4 + # + t_month_2 + t_month_3 + t_month_4 + 
                          latitude , 
                   data = ds_leaf, method = "MM") # #summary(mmfit_leaf)
      # ggplot(ds_leaf, aes(x = latitude, y = t_month_2)) + geom_point() + theme_bw()
      # ds_leaf %>% dplyr::select(t_month_1, t_month_2, t_month_3, t_month_4, latitude) %>% cor()
      # summary(mmfit_leaf)
      # hist(mmfit_leaf$rweights)
    
        # #if the model doesn't converge, use average spring temperature instead
        # focal_fit_spring_temp_leaf <- lmrob(leaf_mean ~ t_month_1_4 + latitude ,
        #                                     data = ds_leaf, method = "MM") #+ t_month_2 #summary(mmfit_leaf)
        # model_temper_param_leaf <- focal_fit_spring_temp_leaf$converged
        # if(model_temper_param_leaf == FALSE){mmfit_leaf <- focal_fit_spring_temp_leaf}
  
  #extract the weights and residuals from the lmrob and join them back to the leaf data
    ds2_leaf <- ds_leaf %>% mutate(weights_leaf = mmfit_leaf$rweights, 
                         residuals_leaf = mmfit_leaf$residuals) %>% 
      dplyr::select(genus, species, individual_id, year_obs, weights_leaf, residuals_leaf,
             t_month_1, t_month_2, t_month_3, t_month_4) %>% 
      left_join(lfp, .)
  
  ### find flow outliers and extract their weights
   ds_flow <- indiv_flow_t %>% #bind_rows(indiv_flow_t, d_flowmax) %>% #add Detroit data from DK
              filter(genus == focal_genus, species == focal_species) %>% filter(!is.na(t_month_1)) #%>% filter(gdd > 0) #ggplot(ds, aes(x = gdd, y = leaf_mean)) + geom_point() + theme_bw()

  #fit a robust regression using the MM estimator
    mmfit_flow <- lmrob(flow_mean ~  t_month_1_4 + #t_month_2 + t_month_3 + t_month_4 + 
                          latitude ,
                      data = ds_flow, method = "MM") #+ t_month_2 #summary(mmfit_leaf)

      # #if the model doesn't converge, use average spring temperature instead
      # focal_fit_spring_temp_flow <- lmrob(flow_mean ~ t_month_1_4 + latitude ,
      #                                     data = ds_flow, method = "MM") #+ t_month_2 #summary(mmfit_leaf)
      # model_temper_param_flow <- focal_fit_spring_temp_flow$converged
      # if(model_temper_param_flow == FALSE){mmfit_flow <- focal_fit_spring_temp_flow}
    
    
  #passing off the weights and residuals from the lmrob
    ds2_flow <- ds_flow %>% mutate(weights_flow = mmfit_flow$rweights,
                                   residuals_flow = mmfit_flow$residuals) %>%
      dplyr::select(genus, species, individual_id, year_obs, weights_flow, residuals_flow) %>%
      left_join(lfp, .)
  
  ## join flowers and leaf
    ds3 <- left_join(ds2_leaf, ds2_flow) %>% 
      filter(genus == focal_genus, species == focal_species) %>% 
      filter(longitude > -90 & longitude < -60) %>% 
      filter(n_obs_per_person > 10) %>% 
      mutate(dif_leaf_flow = flow_mean - leaf_mean ) %>% 
      mutate(temp_outlier = case_when(weights_leaf < 0.8 ~ "leaf outlier",
                                    weights_flow < 0.8 ~ "flower outlier",
                                    .default = "not an outlier"))
  
  ## add Detroit data
    d_leafmax_join <- dplyr::select(d_leafmax, -intensity, -date_mean)  
    d_flowmax_join <- dplyr::select(d_flowmax, -intensity, -date_mean, -Species)              
    d_leafflowmax_join <- left_join(d_leafmax_join, d_flowmax_join) %>% 
      filter(genus == focal_genus, species == focal_species) %>% 
      filter(!is.na(flow_mean)) %>% 
      mutate(temp_outlier = "direct observations",
             weights_leaf = 1,
             weights_flow = 1,
             dif_leaf_flow = flow_mean - leaf_mean )
    
   # ds3 <- bind_rows(ds3, d_leafflowmax_join)
      
   
    
  #calculate additional variables
  ds4 <- ds3  %>% filter(!is.na(t_month_1)) %>% 
          rowwise() %>%  mutate(t_month_1_4 = mean(c(t_month_1, t_month_2, t_month_3, t_month_4))) %>%  #create all spring temp var
          mutate(weights_lf = min(weights_leaf, weights_flow)) %>% ungroup() 
    
  ## visualize outliers of flowers or leaves vs temperatures and save figure
  flow_leaf_outliers_plot <- ds4 %>% 
    ggplot(aes(x = leaf_mean, y = flow_mean, color = weights_lf)) + geom_point(alpha = 0.8, size = 3) + 
    ggtitle(paste(focal_genus, focal_species)) +
    geom_abline(slope = 1, intercept = 0, lty = 2) + ggthemes::theme_few() + 
    scale_color_viridis_c(name = "weights", direction = -1)+ 
    xlab("peak leaf expansion (day of year)") + ylab("peak flowering (day of year)") +
    geom_point(data = d_leafflowmax_join, aes(x = leaf_mean, y = flow_mean ), color = "darkred", shape = 3)

  flow_leaf_outliers_plot_fig_title <- paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/NPN_flower_leaves/", "fig_flow_leaf_outlier_",
                                              focal_genus, "_", focal_species, ".jpg")
  ggsave(flow_leaf_outliers_plot_fig_title, flow_leaf_outliers_plot)
  
  # #CDF visualization
  # ds3 %>%  # filter(weights_leaf > .8) %>%   filter(weights_flow > .8) %>% 
  #   ggplot(aes(x = dif_leaf_flow, color = temp_outlier))+   stat_ecdf(geom = "step") + theme_bw() #geom_histogram() 
  #  ggplot(ds3, aes(x = t_month_1_4, y = dif_leaf_flow, size = weights_leaf)) + geom_point() + theme_bw()
  
  
  
  focal_fit <- lmrob(dif_leaf_flow ~ leaf_mean +
            t_month_1_4 + # + t_month_2 + t_month_3 + t_month_4 +  
            elevation_in_meters + latitude 
            , data = ds4, weights = weights_lf)
  #summary(focal_fit)
  
      # #if the model doesn't converge, use the average spring temperature instead
      # focal_fit_spring_temp <- lmrob(dif_leaf_flow ~ leaf_mean 
      #                    + t_month_1_4 +  elevation_in_meters + latitude 
      #                    , data = ds4, weights = weights_leaf)
      # 
      # model_temper_param <- focal_fit$converged
      # if(model_temper_param == FALSE){focal_fit <- focal_fit_spring_temp}
  
    ## save the species model for use in the NYC_tree_flow.R script
      write_rds(focal_fit, paste0("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/NPN_flower_leaves/", "NPN_model_",
                                   focal_genus, "_", focal_species, ".rds"))
  
  # quantile(focal_fit$residuals, probs = c(.10, .25, .75, .90))
  # quantile(focal_fit$fitted.values, probs = c(.10, .25, .75, .90))
  
  flow_dif_preds <- predict(object = focal_fit, newdata = ds4, interval = "prediction", level = 0.95,
                               se.fit = TRUE)
    
  ds5 <- ds4 %>% mutate(
    pred_dif_leaf_flow = flow_dif_preds$fit[,1],
    pred_dif_leaf_flow_lwr = flow_dif_preds$fit[,2],
    pred_dif_leaf_flow_upr = flow_dif_preds$fit[,3])
  
  #predicted vs observed
  ggplot(ds5, aes(x =  pred_dif_leaf_flow, y = dif_leaf_flow )) + geom_point() + geom_abline(slope = 1, intercept = 0) + theme_bw() 

  #flowering time vs leaf out time for the full dataset with prediction intervals
  ggplot(ds5, aes(x =  leaf_mean, y = flow_mean )) + geom_point() + geom_abline(slope = 1, intercept = 0) + theme_bw() +
    geom_pointrange(aes(x = leaf_mean, 
                        y = leaf_mean + pred_dif_leaf_flow, 
                        ymin = leaf_mean + pred_dif_leaf_flow_lwr, 
                        ymax = leaf_mean + pred_dif_leaf_flow_upr), col = "red", alpha = 0.3)
  
  

### creating a prediction for an NYC tree in 2024 --------------------------------------------------------------------------
  ## extract temperature data for that tree
  #monthly t_mean raster stack
  tmean_rast_d <- pd_stack(prism_archive_subset(temp_period = "monthly", type = "tmean", mon = c(1,2,3,4), years = 2024))
  
  #coords for nyc to extract from raster
  nyc_coords <- data.frame(longitude = -73.97, latitude = 40.78) %>% 
    st_as_sf(coords = c( "longitude", "latitude"), crs = 4326) 
  
  #extract temperature for nyc for that year in those months
  t_mean_nyc_months <- as.vector(raster::extract(x = tmean_rast_d, y = nyc_coords))
  
  #get leaf_mean median for nyc in 2024 for focal species
  leaf_median_nyc2024 <- as.numeric(nyc_sos_summary %>% filter(species == paste(focal_genus, focal_species)) %>% dplyr::select(med_sos))
    
  #fields required for prediction
  nyc_example_tree_pred_data <- data.frame(leaf_mean = leaf_median_nyc2024, #temporary placeholder
                                      t_month_1 = t_mean_nyc_months[1],
                                      t_month_2 = t_mean_nyc_months[2],
                                      t_month_3 = t_mean_nyc_months[3],
                                      t_month_4 = t_mean_nyc_months[4],
                                      elevation_in_meters = 30, 
                                      latitude = 40.78) %>% 
    rowwise() %>%  mutate(t_month_1_4 = mean(c(t_month_1, t_month_2, t_month_3, t_month_4))) %>% ungroup() #create all spring temp var
  
  #create prediction
  flow_dif_preds_nyc <- predict.lm(object = focal_fit, newdata = nyc_example_tree_pred_data, interval = "prediction", level = 0.5,
                                   se.fit = TRUE)
  
###  table X: empirical differences between flowering times and leaf out times in NPN and predicted in NYC in 2024 ------------
    #note that this is embedded within the species loop
    table_lf_focal <- data.frame(genus = focal_genus, species = focal_species, nobs = nrow(ds5), 
                         lf_dif_global_emp_mean = round(mean(ds5$dif_leaf_flow), 3),
                         lf_dif_gloab_emp_p25 = round(quantile(ds5$dif_leaf_flow, 0.25), 3),
                         lf_dif_gloab_emp_p75 = round(quantile(ds5$dif_leaf_flow, 0.75), 3),
                         lf_dif_r2 = round(summary(focal_fit)$r.squared, 2),
                         lf_dif_pred_nyc_2024_mean = round(flow_dif_preds_nyc$fit[1], 2),
                         lf_dif_pred_nyc_2024_p25 = round(flow_dif_preds_nyc$fit[,2], 2), 
                         lf_dif_pred_nyc_2024_p75 = round(flow_dif_preds_nyc$fit[,3], 2))
  table_lf_all_list[[focal_sp_i]] <- table_lf_focal

### SI table X: model fit details without outliers -----------------------------------------------------------------------------
  #note that this is embedded within the species loop
  table_si_start_focal <- data.frame(
    genus = focal_genus, 
    species = focal_species, 
    nobs = nrow(ds5),
    outliers_leaf_n = length(ds3$temp_outlier[ds3$temp_outlier == "leaf outlier"])
    #outliers_flow_n = length(ds3$temp_outlier[ds3$temp_outlier == "flower outlier"])
    )
  
  tidy_focal_fit <- broom::tidy(focal_fit) %>% 
    mutate(est_se = paste0(round(estimate, 3), " (", round(std.error, 3), ")")) %>% 
    dplyr::select(term, est_se) %>% 
    mutate(term = case_when(term == "(Intercept)" ~ "intercept", 
                            .default = term)) %>% 
    pivot_wider(., names_from = term, values_from = est_se)
  
  table_si_focal <- bind_cols(table_si_start_focal, tidy_focal_fit)
  
  table_si_all_list[[focal_sp_i]] <- table_si_focal
  
} ### end species loop


### combine individual rows from each species into a single dataframe
  table_si_all <- bind_rows(table_si_all_list)
  write_csv(table_si_all, "C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/NPN_flower_leaves/table_si_all_250709_weights_without_D.csv")
 # table_si_all <- read_csv("C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/NPN_flower_leaves/table_si_all_250709_weights_without_D.csv")
  
  table_lf_all <- bind_rows(table_lf_all_list) %>% tibble::remove_rownames()
  write_csv(table_lf_all, "C:/Users/dsk273/Box/Katz lab/NYC/tree_pheno/NPN_flower_leaves/table_lf_all_250709_weights_without_D.csv")
  
  

sum(table_si_all$nobs)













  
  
  
  
  
  

#### ad hoc data visualization ########################################################
# test <- lf %>% mutate(lf_dif = leaf_mean - flow_mean) %>% 
#   arrange(lf_dif) %>% 
#   filter(longitude > -90 & longitude < -60) %>% 
#   filter(tmean_spring > 3 & tmean_spring < 9) 

# test2 <- filter(npn_direct, individual_id == 132692) %>% 
#   filter(phenophase_id == 467)

lfp %>% 
  #filter(latitude > 40.7 - 3 & latitude < 40.7 + 3) %>% 
  #filter(longitude > -74 - 3 & longitude < -74 + 3) %>% 
  #filter(individual_id == 4354) %>% 
  filter(longitude > -90 & longitude < -60) %>% 
  filter(n_obs_per_person > 10) %>% 
  #filter(genus == "Quercus", species == "rubra") %>% 
  ggplot(aes(x = leaf_mean, y = flow_mean, color = species)) + geom_point(alpha = 0.2) +
  facet_wrap(~genus) + geom_abline(slope = 1, intercept = 0, lty = 2) + ggthemes::theme_few() +
   geom_smooth(method = lm, se = TRUE)

lfp_sub <- lfp %>% 
  filter(longitude > -90 & longitude < -60) %>% 
  filter(n_obs_per_person > 10) %>% 
  filter(genus == "Quercus", species == "rubra") 
fit <- lm(flow_mean ~ leaf_mean + tmean_spring  , data = lfp_sub)
summary(fit)

lfp %>% 
 # filter(latitude > 40.7 - 3 & latitude < 40.7 + 3) %>% 
 # filter(longitude > -74 - 3 & longitude < -74 + 3) %>% 
  filter(longitude > -90 & longitude < -60) %>% 
  ggplot(aes(x = flow_early - leaf_early, fill = species)) + geom_histogram() +  facet_wrap(~genus, scales = "free") + theme_bw() + geom_vline(xintercept = 0, col = "red", lty = 2)


lfp %>% 
#  filter(latitude > 40.7 - 3 & latitude < 40.7 + 3) %>% 
 # filter(longitude > -74 - 3 & longitude < -74 + 3) %>% 
  filter(longitude > -90) %>% 
  ggplot(aes(x = longitude, y = latitude, col = flow_early - leaf_early)) + geom_point(size = 3) + theme_bw() + facet_wrap(~genus) +
  scale_color_viridis_c()

indiv_leafout %>% 
  filter(longitude > -90) %>% 
  ggplot(aes(x = tmean_spring , y = leaf_early, color = species)) + geom_point() + theme_bw() + facet_wrap(~genus) +
  #scale_color_viridis_d() + 
  #geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE)
  geom_smooth(method = "lm") 


lfp %>% 
  filter(longitude > -90) %>% 
  ggplot(aes(x = tmean_spring , y = flow_early )) + geom_point() + theme_bw() + facet_wrap(~genus, scales = "free") +
  scale_color_viridis_c() + geom_smooth(method = "lm")

lfp %>% 
  filter(longitude > -90) %>% 
  filter(genus == "Quercus") %>% 
  filter(tmean_spring > 3 & tmean_spring < 9) %>% 
  ggplot(aes(x = leaf_early , y = flow_early, col = tmean_spring )) + geom_point() + theme_bw() + facet_wrap(~species, scales = "free") +
  scale_color_viridis_c() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, lty= 2)

lfp %>% 
  filter(longitude > -90) %>% 
  filter(genus == "Acer") %>% 
   ggplot(aes(x = leaf_early, y = flow_early )) + geom_point() + theme_bw() + facet_wrap(~species) +
  scale_color_viridis_c() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, lty= 2)




lfp %>% 
  filter(longitude > -90) %>% 
  ggplot(aes(x = longitude, y = latitude, col = tmean_spring)) + geom_point(size = 2, alpha = 0.5) + theme_bw() + 
  scale_color_viridis_c()


lfp %>% 
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
  
fitted_models <- lfp %>% 
  filter(n_obs_per_person > 100) %>% 
  mutate(gen_sp = paste(genus, species)) %>% 
  nest(data = -gen_sp) %>% 
  mutate(mod1 = map(data, ~lm(flow_mean ~ leaf_mean, data = .)), tidied = map(mod1, tidy)) %>% unnest(tidied)

glance <- fitted_models %>% mutate(tidy = map(mod1, broom::tidy),
                                glance = map(mod1, broom::glance),
                                augment = map(mod1, broom::augment),
                                rsq = glance %>% map_dbl('r.squared'),
                                mse = glance %>% map_dbl('mse'),
                                nobs = glance %>% map_dbl('nobs'),
                                slope = tidy %>% map_dbl(function(x) x$estimate[2]))
glance %>% 
  filter(gen_sp %in% c("Quercus rubra", "Acer platanoides", "Populus deltoides"))

glance$glance

lfp %>% 
  filter(genus == "Acer" & species == "rubrum") %>% 
  filter(longitude > -90) %>% 
  filter(greenup_0 != -9999) %>% 
  ggplot(aes(x = greenup_0, y = flow_mean)) + geom_point()

fitted_models
fitted_models$model

fitted_models %>% tidy(model)
rowwise(fitted_models) %>% tidy(model)








### old stuff I'm not ready to delete yet ##################################################################

# ### data visualization for finding outliers of leaf out in the broader leaf out dataset ###########################
# # assessing which leaf out observations are unusual for one species
# ds <- indiv_leafout_t %>% filter(genus == "Quercus", species == "rubra") #%>% filter(gdd > 0)
# #ggplot(ds, aes(x = gdd, y = leaf_mean)) + geom_point() + theme_bw()
# 
# #fit a robust regression using the MM estimator
# mmfit <- lmrob(leaf_mean ~  poly(t_month_1, 2) + poly(t_month_2, 2) + poly(t_month_3, 2) + poly(t_month_4, 2) + latitude , 
#                data = ds, method = "MM") #+ t_month_2
# summary(mmfit)
# 
# #some additional visualization
# # mmfit_df <- data.frame(leaf_mean = mmfit$model[,1], t_mean_3 = mmfit$model[,4], # t_mean_2 = mmfit$model[,3],
# #                        resids = mmfit$residuals, weights = mmfit$rweights)
# # ggplot(mmfit_df, aes(x = resids, y = weights, color = weights)) + geom_point() + theme_bw()
# # ggplot(mmfit_df, aes(x = t_mean_3, y = leaf_mean, color = weights)) + geom_point() + theme_bw()
# # hist(mmfit$residuals)
# 
# #passing off the weights and residuals from the lmrob
# ds2 <- ds %>% mutate(weights = mmfit$rweights, 
#               residuals = mmfit$residuals) %>% 
#   select(genus, species, individual_id, year_obs, weights, residuals)
# 
# ds3 <- left_join(lfp, ds2)
# 
# ds3 %>% 
#   filter(longitude > -90 & longitude < -60) %>% 
#   filter(n_obs_per_person > 10) %>% 
#   filter(genus == "Quercus", species == "rubra") %>% 
#   filter(weights > .75) %>% 
#   ggplot(aes(x = leaf_mean, y = flow_mean, color = weights)) + geom_point() +
#   facet_wrap(~genus) + geom_abline(slope = 1, intercept = 0, lty = 2) + ggthemes::theme_few() + scale_color_viridis_c()+
#   geom_smooth(method = lm, se = TRUE)
# 
# ds3 %>% 
#   filter(longitude > -90 & longitude < -60) %>% 
#   filter(n_obs_per_person > 10) %>% 
#   filter(genus == "Quercus", species == "rubra") %>% 
#   mutate(dif_leaf_flow = leaf_mean - flow_mean) %>% 
#   filter(weights > .8) %>% 
#   ggplot(aes(x = dif_leaf_flow))+ #geom_histogram() +
#   stat_ecdf(geom = "step")
# 
# ds4 <- ds3 %>% 
#   filter(longitude > -90 & longitude < -60) %>% 
#   filter(n_obs_per_person > 10) %>% 
#   filter(genus == "Quercus", species == "rubra") %>%  filter(weights > .8) 
#   
# fit <- lm(flow_mean ~ leaf_mean , data = ds4)
# summary(fit)
#   
# quantile(fit$residuals, probs = c(.10, .25, .75, .90))
# 
# 
# 
# 
# 
# 
# ### data visualization for finding outliers of flow in the broader flow dataset ###########################
# # assessing which flow observations are unusual for one species
# ds <- indiv_flow_t %>% filter(genus == "Quercus", species == "rubra") #%>% filter(gdd > 0)
# #ggplot(ds, aes(x = t_month_3, y = flow_mean)) + geom_point() + theme_bw()
# 
# #fit a robust regression using the MM estimator
# mmfit_flow <- lmrob(flow_mean ~  poly(t_month_1, 2) + poly(t_month_2, 2) + poly(t_month_3, 2) + poly(t_month_4, 2) + latitude , 
#                data = ds, method = "MM") #+ t_month_2
# summary(mmfit_flow)
# 
# #passing off the weights and residuals from the lmrob
# ds2 <- ds %>% mutate(weights_flow = mmfit_flow$rweights, 
#                      residuals_flow = mmfit_flow$residuals) %>% 
#   select(genus, species, individual_id, year_obs, weights_flow, residuals_flow)
# 
# ds3 <- left_join(lfp, ds2)
# 
# ds3 %>% 
#   filter(longitude > -90 & longitude < -60) %>% 
#   filter(n_obs_per_person > 10) %>% 
#   filter(genus == "Quercus", species == "rubra") %>% 
#   #filter(weights > .75) %>% 
#   ggplot(aes(x = leaf_mean, y = flow_mean, color = weights_flow)) + geom_point() +
#   facet_wrap(~genus) + geom_abline(slope = 1, intercept = 0, lty = 2) + ggthemes::theme_few() + scale_color_viridis_c()+
#   geom_smooth(method = lm, se = TRUE)
# 
# ds3 %>% 
#   filter(longitude > -90 & longitude < -60) %>% 
#   filter(n_obs_per_person > 10) %>% 
#   filter(genus == "Quercus", species == "rubra") %>% 
#   mutate(dif_leaf_flow = leaf_mean - flow_mean) %>% 
#   #filter(weights_flow > .8) %>% 
#   ggplot(aes(x = dif_leaf_flow))+ #geom_histogram() +
#   stat_ecdf(geom = "step")
# 
# ds4 <- ds3 %>% 
#   filter(longitude > -90 & longitude < -60) %>% 
#   filter(n_obs_per_person > 10) %>% 
#   filter(genus == "Quercus", species == "rubra") %>%  filter(weights_flow > .8) 
# 
# fit <- lm(flow_mean ~ leaf_mean , data = ds4)
# summary(fit)
# 
# quantile(fit$residuals, probs = c(.10, .25, .75, .90))
