calculate_x_for_polar = function(.data){
  # Create annotation
  internal_branch_length = 40
  external_branch_length = 10
  
  
  # Integrate data
  tree_df_source = 
    ARMET::tree %>%
    data.tree::ToDataFrameTree("name", "isLeaf", "level", "leafCount", pruneFun = function(x)	x$level <= 5) %>%
    as_tibble() %>%
    rename(`Cell type category` = name) %>%
    mutate(level = level -1)
  
  # Calculate x
  map_df(
    0:4,
    ~ tree_df_source %>%
      filter(level == .x | (level < .x & isLeaf)) %>%
      mutate(leafCount_norm = leafCount/sum(leafCount)) %>%
      mutate(leafCount_norm_cum = cumsum(leafCount_norm)) %>%
      mutate(length_error_bar = leafCount_norm - 0.005) %>%
      mutate(x = leafCount_norm_cum - 0.5 * leafCount_norm) 	
  ) %>%
    
    # Attach data
    left_join(.data, by = c("Cell type category", "level")) %>%
    
    # process
    mutate(Estimate = ifelse(significant, fold_change, NA)) %>%
    mutate(branch_length = ifelse(isLeaf, 0.1, 2)) %>%
    
    # Correct branch length
    mutate(branch_length = ifelse(!isLeaf, internal_branch_length,	external_branch_length) ) %>%
    mutate(branch_length =  ifelse(	isLeaf, branch_length + ((max(level) - level) * internal_branch_length),	branch_length	)) 
  
}

magnolia_layer = function( gg, xx, .level, prop_filter, size_geom_text, linetype = "solid"){
   
  gg + 
  geom_bar(
    data = xx %>% filter(level == .level),
    aes(width = leafCount_norm, y = `median_proportion` ), # / leafCount_norm),
    color = "grey20", stat = "identity"
  )  +
    {
      # Make plotting robust if no level 3 cell types were detected
      xx %>% 
        when(
        prop_passes_threshold(., .level, `median_proportion`, Estimate, prop_filter) ~ 
          geom_errorbar(
            data = filter_lowly_abundant(., .level, `median_proportion`, Estimate, prop_filter),
            aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[.level:(.level+1)])), # / leafCount_norm),
            color = "grey20",  stat = "identity",linetype=linetype
          ) ,
        
        # If there are NOT cell types
        ~ geom_errorbar(ymin=0, ymax=0)
      )
    } +
    {
      # Make plotting robust if no level 3 cell types were detected
      xx %>% 
        when(
          prop_passes_threshold(., .level, `median_proportion`, Estimate, prop_filter) ~ 
        
        geom_text(
          data = filter_lowly_abundant(., .level, `median_proportion`, Estimate, prop_filter) ,
          aes(label=`formatted`, y = mean(soc[.level:(.level+1)]) , angle= angle) ,size =size_geom_text
        ),
        
        # If there are NOT cell types
       ~ geom_errorbar(ymin=0, ymax=0)
      )
    }
  
  
}
