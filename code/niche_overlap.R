##########################################
####  Niche overlap function
##########################################
#### | Project name: Bioclimatic envelopes
#### | Creator: Mirza Cengic
#### | Contact: mirzaceng@gmail.com
##########################################

#### Define function that calculates hypervolume overlap
#' Niche overlap
#'
#' @param species_data 
#' @param var_set 
#' @param mode 
#'
#' @return
#' @export
#'
#' @examples
#' 
niche_overlap <- function(species_data, var_set = "4var", mode)
{
  stopifnot(var_set %in% c("2var", "4var", "10var", "allvar"), mode %in% c("cv", "ex"), !missing(species_data))
  
  bioclim_2var <- paste0("bio_", str_pad(c(1, 12), width = 2, pad = "0", side = "left"))
  bioclim_4var <- paste0("bio_", str_pad(c(1, 4, 12, 15), width = 2, pad = "0", side = "left"))
  bioclim_10var <- paste0("bio_", str_pad(c(2:4, 8:9, 13:15, 18:19), width = 2, pad = "0", side = "left"))
  bioclim_allvar <- paste0("bio_", str_pad(1:19, width = 2, pad = "0", side = "left"))
  
  cat(glue("Calculating overlap percentage on {str_replace(var_set, 'var', '')} axes for {unique(species_data$binomial)} - {mode}."), "\n")
  
  my_cols <- paste0("bioclim_", var_set)
  
  if (mode == "cv")
  {
    
    species_train_data <- species_data %>% 
      filter(Set == "train")
    
    set.seed(666)
    my_sample <- sample(1:nrow(species_train_data), round(0.8 * nrow(species_train_data)), replace = FALSE)
    
    train_data_cv <- species_train_data[my_sample, ]
    test_data_cv  <- species_train_data[-my_sample, ]
    #
    hypervolume_train_presences <- train_data_cv %>% 
      filter(PA == 1) %>% 
      select(get(my_cols)) %>% 
      hypervolume_box(name = paste0("train", var_set), samples.per.point = 10)
    
    hypervolume_test_presences <- test_data_cv %>% 
      filter(PA == 1) %>% 
      select(get(my_cols)) %>% 
      hypervolume_box(name = paste0("test", var_set), samples.per.point = 10)
      
  }
  if (mode == "ex")
  {
  
    # Presences
    hypervolume_train_presences <- species_data %>% 
      filter(Set == "train") %>% 
      filter(PA == 1) %>% 
      select(get(my_cols)) %>% 
      hypervolume_box(name = paste0("train", var_set), samples.per.point = 10)
    
    hypervolume_test_presences <- species_data %>% 
      filter(Set == "test") %>% 
      filter(PA == 1) %>% 
      select(get(my_cols)) %>% 
      hypervolume_box(name = paste0("test", var_set), samples.per.point = 10)
    
  }
  
 # Calculate values
  hv_comparison_presences <- hypervolume_set(hypervolume_train_presences, 
                                             hypervolume_test_presences, 
                                             check.memory = FALSE) %>% 
    hypervolume_overlap_statistics() %>% 
    melt() %>% 
    rownames_to_column() %>% 
    transmute(
      index = rowname,
      value,
      PA_set = "presences"
    )
  return(hv_comparison_presences)
  
}