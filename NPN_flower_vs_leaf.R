#This script is for 


### select a few top anemophilous taxa from NPN #####################################
#anemophilous angiosperms
acer_species_list <- c(777,1843,59,778,1,2,1591,60,779,780,3,781,61,1199)
alnus_species_list <- c(62,63,319)
betula_species_list <- c(97, 1439, 98, 430, 1850, 1339, 1851, 99, 1805)
fraxinus_species_list <- c(74,872,873,75,1350)
populus_species_list <- c(1361,320,976,977,1188,27,1481)
quercus_species_list <- c(705,100,1365,757,1870,987,1690,1484,988,316,297,1485,1190,765,1486,
                          301,704,101,1691,1212,989,1366,102,1756,1213,1755,1487,1159,305)
ulmus_species_list <- c(1192,1048,1049,1215,1216)

#herbaceous angiosperms
ambrosia_species_list <- c(145,788,146)

#pollen cones
pinus_species_list <- c(1629,1686,965,762,50,295,220,1480,219,51,966,967,52,25,968,1687,53,54)
cupressaceae_species_list <- c(43,1743,1354,289,902,291,290,44, #junipers
                               1723,831,1356,296,1040) #other

list_all_focal_taxa <- c(acer_species_list, alnus_species_list, betula_species_list, fraxinus_species_list, 
                         populus_species_list, quercus_species_list, ulmus_species_list,
                         ambrosia_species_list, pinus_species_list, cupressaceae_species_list)

###download and process data #####################################
npn_direct <- npn_download_status_data(
  request_source = 'Daniel Katz, UT Austin and/or Theresa Crimmins',
  species_ids = list_all_focal_taxa,
  years = c(as.character(2008:2021)), #years to include
  phenophase_ids = c(501, 502,495, 503), #angiosperms: 501 == "Open flowers", 502 == "Pollen release (flowers)" #conifers: 495 ==  503 ==
  additional_fields = c("Observed_Status_Conflict_Flag", "partner_group")
)
#names(npn_direct)

npn_direct <- subset(npn_direct, observed_status_conflict_flag == "-9999")

npn_direct <- filter(npn_direct, !(partner_group %in% c("CSU Chico NSCI 102", "SCMF Naturalists", 
                                                        "Sycamore Canyon Wilderness Park - Riverside",
                                                        "Marist College", "Sam Hughes Neighborhood Association",
                                                        "UNCO BIO 111", "Maricopa Cooperative Extension",
                                                        "Pima County Extension", "Lasell College",
                                                        "UofL campus", "Ursinus College", "U of A Campus Arboretum",
                                                        "RMC Campus Phenology", "SUNY Geneseo", "AZ Project WET")))

npn_flow <- filter(npn_direct, phenophase_id == 501 | phenophase_id == 495)

# looks like flowering intensity value was only entered very rarely
npn_active_flow <- npn_flow %>%
  filter(phenophase_status != -1) %>% #removing observations where the observer was unsure whether the phenophase was occurring
  mutate(flow_prop = case_when(
    phenophase_status == 0  ~ 0,
    intensity_value == "Less than 5%" ~ 0.025,
    intensity_value ==  "5-24%"~ (0.05+0.24)/2,
    intensity_value == "25-49%" ~ (0.25+0.49)/2,
    intensity_value == "50-74%" ~ (0.5+0.74)/2,
    intensity_value == "75-94%" ~ (0.75+0.94)/2,
    intensity_value == "95% or more" ~ 0.97,
    phenophase_status == 1  ~ 0.5)) #assuming that when intensity value isn't given, a tree is halfway through flowering
hist(npn_active_flow$flow_prop)

#how many observations are for pollen release?                    
filter(npn_direct, phenophase_id == 502 | phenophase_id == 503) %>% group_by(phenophase_status) %>% summarize(n_obs = n())
