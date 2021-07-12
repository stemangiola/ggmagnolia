
filter_lowly_abundant = function(.data, .level, y, value, theshold){
  
  y = enquo(y)
  value= enquo(value)
  
  .data %>% filter(level == .level) %>% filter(!!y>theshold | !is.na(!!value))
  
}

prop_passes_threshold = function(.data, .level, y, value, theshold){
  
  y = enquo(y)
  value= enquo(value)
  
  .data %>% filter_lowly_abundant(.level, !!y, !!value, theshold) %>% nrow() > 0
  
}


# Greater than
gt = function(a, b){	a > b }

# Smaller than
st = function(a, b){	a < b }

# Negation
not = function(is){	!is }

# Raise to the power
pow = function(a,b){	a^b }
