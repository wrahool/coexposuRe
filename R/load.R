sample_atleast_once <- function(x, n){
  
  # Only consider unique items
  if(length(unique(x)) > n){
    stop("Not enough unique items in input to give at least one of each")
  }
  
  # Get values for vector - force each item in at least once
  # then randomly select values to get the remaining.
  vals <- c(unique(x),
            sample(unique(x), n - length(unique(x)), 
                   replace = TRUE))
  
  # Now shuffle them
  sample(vals)
}

# function to calculate the NMI for community structure c
get_NMI <- function(c, outlet_types, v) {
  
  # a warning message here such as
  # Warning message:
  # In if (!is.na(c)) { :
  #    the condition has length > 1 and only the first element will be used
  # is perfectly fine, by design.
  # if there is an error in community detection, the value defaults to NA
  # which this if condition catches.
  # if the community detection works, then it returns a list of things, none of
  # which are NA.
  
  if(!is.na(c)) {
    confusion_tbl <- outlet_types %>%
      merge(tibble(
        outlet_name = c$names,
        pred_type = c$membership
      ))
    
    NMI_score <- aricode::NMI(confusion_tbl$outlet_type,
                              confusion_tbl$pred_type, variant = v)
    
    if(is.nan(NMI_score)) # happens with variant = min, and variant = sqrt when algorithm only produces 1 community
      NMI_score <- 0
  } else
    NMI_score <- NA
  
  return(NMI_score)
  
}

# function to calculate the AMI for community structure c
get_AMI <- function(c, outlet_types) {
  
  # a warning message here such as
  # Warning message:
  # In if (!is.na(c)) { :
  #    the condition has length > 1 and only the first element will be used
  # is perfectly fine, by design.
  # if there is an error in community detection, the value defaults to NA
  # which this if condition catches.
  # if the community detection works, then it returns a list of things, none of
  # which are NA.
  
  if(!is.na(c)) {
    confusion_tbl <- outlet_types %>%
      merge(tibble(
        outlet_name = c$names,
        pred_type = c$membership
      ))
    
    AMI_score <- aricode::AMI(confusion_tbl$outlet_type,
                              confusion_tbl$pred_type)
  } else
    AMI_score <- NA
  
  return(AMI_score)
  
}

# function to calculate the scaling factor
get_scalingfactor <- function(c, outlet_types) {
  
  # scaled NMI using a scaling factor to penalize algorithms that produce a large number of communities
  # a warning message here such as
  # Warning message:
  # In if (!is.na(c)) { :
  #    the condition has length > 1 and only the first element will be used
  # is perfectly fine, by design.
  # if there is an error in community detection, the value defaults to NA
  # which this if condition catches.
  # if the community detection works, then it returns a list of things, none of
  # which are NA.
  
  
  if(!is.na(c)) {
    confusion_tbl <- outlet_types %>%
      merge(tibble(
        outlet_name = c$names,
        pred_type = c$membership
      ))
    
    NMI_score <- aricode::NMI(confusion_tbl$outlet_type,
                              confusion_tbl$pred_type)
    
    R <- confusion_tbl %>% pull(outlet_type) %>% unique() %>% length()
    S <- confusion_tbl %>% pull(pred_type) %>% unique() %>% length()
    
    scaling_factor <- exp(-(abs(R-S))/R)
    # SNMI_score <- scaling_factor * NMI_score
    
  } else
    scaling_factor <- NA
  
  return(scaling_factor)
  
}

# function to calculate the scaled NMI for community structure c
get_SNMI <- function(c, outlet_types) {
  
  # scaled NMI using a scaling factor to penalize algorithms that produce a large number of communities
  # a warning message here such as
  # Warning message:
  # In if (!is.na(c)) { :
  #    the condition has length > 1 and only the first element will be used
  # is perfectly fine, by design.
  # if there is an error in community detection, the value defaults to NA
  # which this if condition catches.
  # if the community detection works, then it returns a list of things, none of
  # which are NA.
  
  
  if(!is.na(c)) {
    confusion_tbl <- outlet_types %>%
      merge(tibble(
        outlet_name = c$names,
        pred_type = c$membership
      ))
    
    NMI_score <- aricode::NMI(confusion_tbl$outlet_type,
                              confusion_tbl$pred_type)
    
    R <- confusion_tbl %>% pull(outlet_type) %>% unique() %>% length()
    S <- confusion_tbl %>% pull(pred_type) %>% unique() %>% length()
    
    scaling_factor <- exp(-(abs(R-S))/R)
    SNMI_score <- scaling_factor * NMI_score
    
  } else
    SNMI_score <- NA
  
  return(SNMI_score)
  
}