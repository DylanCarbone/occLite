# Helper function that makes sure occupancy values are greater than 0 and less than 1
bound_zero_one = function(num){

  num = ifelse(num >= 1, 1-1e-6, ifelse(num <= 0, 1e-6, num))

  return(num)
}

# Helper function that makes sure occupancy values are greater than 0 and less than 1
bound_for_logit = function(num){

  num = ifelse(num == 1, 1-1e-6, ifelse(num == 0, 1e-6, num))
  
  if(any(num > 1 | num < 0, na.rm = TRUE)){
    stop("Values are greater than 1 or less than 0")
  }

  return(num)
}