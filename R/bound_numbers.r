# Helper function that makes sure occupancy values are not bounded between 0 and 1 for later logit transformation
bound_numbers = function(num){

  num = ifelse(num == 1, 1-1e-6, ifelse(num == 0, 1e-6, num))
  
  if(any(num > 1 | num < 0)){
    stop("Values are greater than 1 or less than 0")
  }

  return(num)
}