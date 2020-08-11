##########################################
####  Run additional models for maxent
##########################################
#### | Project name: Bioclim envelopes
#### | Creator: Mirza Cengic
#### | Contact: mirzaceng@gmail.com
##########################################

# Script setup ------------------------------------------------------------

# Load packages
pacman::p_load(Rahat, biomod3, raster, tidyverse, hypervolume, sf, sp, fasterize, janitor, glue, reshape2)

# Source additional scripts with functions
source("hypervolume_overlap_percent.R")
source("niche_overlap.R")
source("modeling_functions.R")

# Load data ---------------------------------------------------------------

# Load data with 300 species to model (100 species of fish, amphibians, and mammals). 
processing_species_files <- "Projects/Bioclim_envelopes/Data/Species_point_data" %>%
  milkunize2("archive") %>%
  list.files(full.names = TRUE)

# Load global realms raster. This will be used to limit presence/absence selection (see methods).
# Spatial resolution is 5 arc-minutes.
biomes_raster <- "Biomes_5m.tif" %>% 
  milkunize2() %>% 
  raster()
# Load global basins raster. This is used to delinate species data subsetting regions
basins_5m <- raster("basins_5min.tif") %>% 
  setMinMax()

# Load raster predictors at 5 arc degree spatial resolution -- set with 4 layers
bioclim_stack <- "Data_RAW/Current_climate/5m" %>%
  milkunize2("archive") %>%
  list.files(pattern = "tif", full.names = TRUE) %>%
  stack()

# Change layer names
names(bioclim_stack) <- names(bioclim_stack) %>% 
  str_replace("bio10_", "bio_")

# Load normalized raster predictors at 5 arc degree spatial resolution
# To be used for niche overlap stuff
bioclim_stack_norm <- "Data_RAW/Current_climate/5m/Normalized" %>%
  milkunize2("archive") %>%
  list.files(pattern = "tif", full.names = TRUE) %>%
  stack()

# Run model ####

i <- as.numeric(commandArgs(trailingOnly = TRUE))

species_list <- "Projects/Bioclim_envelopes/Data/All_species_list.csv" %>% 
  milkunize2("archive") %>% 
  read_csv()

input_species_range <- st_read(processing_species_files[i])


# Get species name
species_name <- as.character(unique(input_species_range$binomial))
# Get taxonomic group
taxonomic_group <- species_list %>% 
  filter(species == species_name) %>% 
  pull(group)

cat(glue("Creating data for {species_name}. {length(processing_species_files) - i} more species to go."), "\n")


####alms_speci#### Get predictor set
clim_var_sets <- c("2var", "4var", "10var", "allvar")

bioclim_2var <- paste0("bio_", str_pad(c(1, 12), width = 2, pad = "0", side = "left"))
bioclim_4var <- paste0("bio_", str_pad(c(1, 4, 12, 15), width = 2, pad = "0", side = "left"))
bioclim_10var <- paste0("bio_", str_pad(c(2:4, 8:9, 13:15, 18:19), width = 2, pad = "0", side = "left"))
bioclim_allvar <- paste0("bio_", str_pad(1:19, width = 2, pad = "0", side = "left"))

# Set folder for species

species_data_path <- glue("Projects/Bioclim_envelopes/Output_2019/{taxonomic_group}/{species_name}") %>% 
  milkunize2("archive") 

dir.create(species_data_path, recursive = TRUE)
dir.create(glue("{species_data_path}/Outputs"), recursive = TRUE)


# Count the number of files in the folder
files_number <- species_data_path %>% 
  list.files(recursive = TRUE, pattern = "assessment*.*.csv$") %>% 
  length()

# Main loop ---------------------------------------------------------------
for (mod_eval in c("ex", "cv"))
{
  for (pseudo_method in c("PA1", "PA2", "PA3"))
  {
    for (predictor_set in c("2var", "4var", "10var", "allvar"))
    {
      
      cat(glue("{predictor_set}, {mod_eval}, {pseudo_method}"), "\n")
      
      scenario_name <- glue("{predictor_set}_{pseudo_method}_{mod_eval}")
      
      model_assessment_name <- glue("{species_data_path}/Outputs/{species_name}_{scenario_name}_assessment_maxent.csv")
      
      if (file.exists(model_assessment_name))
      {
        next()
      }
      
      setwd(species_data_path)
      
      fitted_model <- run_model(data_sf = input_species_range, 
                                evaluation_mode = mod_eval,
                                variable_set = predictor_set,
                                pseudoabs_method = pseudo_method,
                                climate_stack = bioclim_stack,
                                which_model = "maxent_bckg")
      
      model_evaluation <- evaluate_model(model = fitted_model,
                                         evaluation_mode = mod_eval,
                                         variable_set = predictor_set,
                                         pseudoabs_method = pseudo_method) %>% 
        mutate(
          Algorithm = "MAXENT.Background"
        )
      
      
      write_csv(model_evaluation, model_assessment_name)
      
    }
  }
}