##########################################
####  Run additional species-specific models
##########################################
#### | Project name: Bioclimatic envelopes
#### | Creator: Mirza Cengic
#### | Contact: mirzaceng@gmail.com
##########################################

# Script setup ------------------------------------------------------------
# Load packages
pacman::p_load(Rahat,
               biomod3,
               raster, tidyverse, usdm,
               # hypervolume,
               sf, sp, fasterize, janitor, glue, reshape2)

# Source additional scripts with functions
source("hypervolume_overlap_percent.R")
source("niche_overlap.R")
source("modeling_functions.R")

# Load data ---------------------------------------------------------------

# Load csv that contains list of 300 species to model (100 species of fish, amphibians, and mammals). 
# Species list is taken from manuscript Supplementary Info
processing_species_files <- "Projects/Bioclim_envelopes/Data/Species_point_data_all" %>%
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

# Fix layer names
names(bioclim_stack) <- names(bioclim_stack) %>%  
  str_replace("bio10_", "bio_")

# Do stuff ----------------------------------------------------------------

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


# Niche overlap -----------------------------------------------------------


#### Get predictor set
# clim_var_sets <- c("2var", "4var", "10var", "allvar")
# Run 5th predictor set instead of 4 regular ones
clim_var_sets <- "vifvar"


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



#####################################################################
############ Develop here VIF subsetting per species ###############

# cat(glue("Calculating overlap for {species_name}, {variable_set}iables, for {evaluation_mode}"), "\n")

species_name <- unique(input_species_range$binomial)

#### Create dataset for model fitting
# Subset presences and absences, sample x number of absences and merge back together
# Subset data for training
species_presences_train <- input_species_range %>% 
  filter(PA == 1)

species_absences_train <- input_species_range %>% 
  filter(PA == 0)


basefolder_vif <- "Projects/Bioclimatic_envelopes/Output/VIF_vars" %>% 
  milkunize2()


for (pseudoabs_method in c("PA1", "PA2", "PA3"))
{
  print(str_glue("Running {species_name} for {pseudoabs_method}"))
  
  if (pseudoabs_method == "PA1")
  {
    pa_number <- 1000
  }
  if (pseudoabs_method == "PA2")
  {
    pa_number <- 10000
  }
  if (pseudoabs_method == "PA3")
  {
    pa_number <- nrow(species_presences_train)
  }
  

  
  #outname_vif3 <- str_glue("{basefolder_vif}/VIF_th3_{pseudoabs_method}_{species_name}.csv")
  outname_vif10 <- str_glue("{basefolder_vif}/VIF_th10_{pseudoabs_method}_{species_name}.csv")
  
  if (!file.exists(outname_vif10))
  {
    
    set.seed(666)
    species_absences_sampled <- sample_n(species_absences_train, pa_number, replace = FALSE)
    # species_absences_sampled <- sample_n(species_absences_train, pa_number, replace = FALSE)
    
    
    species_data_merged <- rbind(species_presences_train, 
                                 species_absences_sampled)
    
    
    
    species_data_bioclim <- raster::extract(bioclim_stack, species_data_merged, sp = TRUE)
    
    # bioclim_stack_subset <- subset(bioclim_stack, get(layers_to_subset))
    
    
    my_predictors <- species_data_bioclim %>% 
      st_as_sf() %>% 
      st_set_geometry(NULL) %>% 
      dplyr::select(-area_id, -n, -Set, -PA, -binomial)
    
    
    #set.seed(666)
    #my_vifs_3 <- usdm::vifstep(my_predictors, th = 3)
    set.seed(666)
    my_vifs_10 <- usdm::vifstep(my_predictors, th = 10)
    
    
    
    #my_df_vif3 <- my_vifs_3@results
    
    #my_df_vif3 <- my_df_vif3 %>% 
    #  mutate(
    #    binomial = species_name,
    #    group = taxonomic_group,
    #    vars_num = nrow(my_df_vif3),
    #    VIF_thresh = 3, 
    #    pa_meth = pseudoabs_method
    #  )
    
    my_df_vif10 <- my_vifs_10@results
    
    my_df_vif10 <- my_df_vif10 %>% 
      mutate(
        binomial = species_name,
        group = taxonomic_group,
        vars_num = nrow(my_df_vif10),
        VIF_thresh = 10, 
        pa_meth = pseudoabs_method
      )
    
    
    
    write_csv(my_df_vif10, outname_vif10)
    #write_csv(my_df_vif3, outname_vif3)
  } else {
    my_df_vif10 <- read_csv(outname_vif10)
    #my_df_vif3 <- read_csv(outname_vif3)
  }
  
  
}

#####


bioclim_vifvar <- unique(as.character(my_df_vif10$Variables))


# Set folder for species

species_data_path <- glue("Projects/Bioclim_envelopes/Output_2019/{taxonomic_group}/{species_name}_VIF10") %>% 
  milkunize2("archive") 

dir.create(species_data_path, recursive = TRUE)
dir.create(glue("{species_data_path}/Outputs"), recursive = TRUE)


# Count the number of files in the folder
 files_number <- species_data_path %>% 
  list.files(recursive = TRUE, pattern = "assessment*.*.csv$") %>% 
  length()


if (files_number == 48)
{
  print("All done")
} else {
  # Niche overlap calc
  # Main loop ---------------------------------------------------------------
  predictor_set = "vifvarth10"
  
  # for (predictor_set in c("2var", "4var"))
  
    for (mod_eval in c("ex", "cv"))
    {
      # for (pseudo_method in "PA1")
      for (pseudo_method in c("PA1", "PA2", "PA3"))
      {
        # for (predictor_set in c("allvar"))
        cat(glue("{predictor_set}, {mod_eval}, {pseudo_method}"), "\n")
        
        scenario_name <- glue("{predictor_set}_{pseudo_method}_{mod_eval}")
        
        model_assessment_name <- glue("{species_data_path}/Outputs/{species_name}_{scenario_name}_assessment.csv")
        model_assessment_ensemble_name <- glue("{species_data_path}/Outputs/{species_name}_{scenario_name}_assessment_ensemble.csv")
        
        if (file.exists(model_assessment_ensemble_name))
        {
          cat("File exists", "\n")
          next()
        }
        
        setwd(species_data_path)
        
        fitted_model <- run_model(data_sf = input_species_range, 
                                  evaluation_mode = mod_eval,
                                  variable_set = bioclim_vifvar,
                                  mod_name = str_glue("{species_name}_VIF10"),
                                  pseudoabs_method = pseudo_method,
                                  climate_stack = bioclim_stack)
        
        ####
        model_evaluation <- evaluate_model(model = fitted_model,
                                           evaluation_mode = mod_eval,
                                           variable_set = predictor_set,
                                           pseudoabs_method = pseudo_method)
        
        
        write_csv(model_evaluation, model_assessment_name)
        
        # Ensemble model
        ensemble_model_out <- BIOMOD_EnsembleModeling(modeling.output = fitted_model,
                                                      chosen.models = "all",
                                                      em.by = "PA_dataset+repet",
                                                      eval.metric = c("CSI", "TSS"), #metric used to scale the ensamble
                                                      eval.metric.quality.threshold = c(-1, -1),
                                                      models.eval.meth = c("CSI", "TSS", "ROC"),
                                                      prob.mean = TRUE,
                                                      prob.mean.weight = TRUE)
        
        ####
        ensemble_model_evaluation <- evaluate_ensemble(model = fitted_model, 
                                                       ensemble_model = ensemble_model_out, 
                                                       evaluation_mode = mod_eval,
                                                       variable_set = predictor_set,
                                                       pseudoabs_method = pseudo_method)
        
        
        write_csv(ensemble_model_evaluation, model_assessment_ensemble_name)  
        
      
    }
  }  
}

