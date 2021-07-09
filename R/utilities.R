
filter_lowly_abundant = function(.data, .level, .value, .estimate, theshold){
  
  .value = enquo(.value)
  .estimate= enquo(.estimate)
  
  .data %>% filter(level == .level) %>% filter(`median_proportion`>theshold | !is.na(Estimate))
  
}

prop_passes_threshold = function(.data, .level, .value, .estimate, theshold){
  
  .value = enquo(.value)
  .estimate= enquo(.estimate)
  
  .data %>% filter_lowly_abundant(.level, !!.value, !!.estimate, theshold) %>% nrow() > 0
  
}


# Greater than
gt = function(a, b){	a > b }

# Smaller than
st = function(a, b){	a < b }

# Negation
not = function(is){	!is }

# Raise to the power
pow = function(a,b){	a^b }
