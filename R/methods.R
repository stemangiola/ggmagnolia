plot_polar_old = function(	.data,
                       size_geom_text = 3.5,
                       my_breaks=c(0, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.7, 1),
                       prop_filter = 0.005,
                       barwidth = 0.5,
                       barheight = 2,
                       legend_justification = 0.67,
                       fill_direction = 1){
  
  xx  = 
    .data %>%
    select(-one_of("draws", "rng", "draws_cens")) %>%
    calculate_x_for_polar %>%
    
    # Add angle text
    mutate(angle = x * -360) %>%
    mutate(angle = ifelse(angle < - 240 | angle > -120, angle, angle - 180)) %>%
    mutate(angle = ifelse(angle > - 360, angle, angle + 360)) %>%
    #mutate(median_proportion = ifelse(level<max(level), median_proportion + (0.008*level), median_proportion  )) %>%
    filter(level > 0) %>%
    
    # Calculate median proportions
    filter(.value_alpha1 %>% is.na %>% `!`) %>%
    mutate(summarised_proportion =
             proportions %>% map(~ .x %>% summarise(median_proportion = .value %>% median))
    ) %>%
    unnest(summarised_proportion) %>%
    #mutate(median_proportion = median_proportion / max(median_proportion)) %>%
    
    distinct() %>%
    left_join(
      ct_names_polar, by=c("Cell type category" = "name")
    )
  
  my_rescale = function(x) { as.character( format(round( x * xx %>% pull(median_proportion) %>% max, 2), nsmall = 2)) }
  cusotm_root_trans = function() scales::trans_new("cusotm_root",function(x) x^(1/4), function(x) x^(4))
  DTC_scale =  xx %>% pull(Estimate) %>% abs() %>% { switch( (!all(is.na(.))) + 1, 1, .) } %>% max(na.rm = T)
  y_text_out = 8.5
  y_text_in = 1.4641
  
  # Size outer circles
  soc = (c(0, 0.2, 0.4, 0.6, 0.8)*0.5 + 1)^4  
  
  xx %>%
    {
      # Case if none is significant
      switch(
        (!length(na.omit(xx$Estimate))>0) + 1,
        
        # If there are significant
        ggplot(data=(.), aes(x = x,fill = Estimate,size = 1/sqrt(level))),
        
        # If there are not
        ggplot(data=(.),aes(x = x,fill = "Non significant",size = 1/sqrt(level)))
      )
    }	+
    annotate("rect", xmin=0, xmax=1, ymin=soc[1], ymax= soc[2], fill="grey95") +
    annotate("rect", xmin=0, xmax=1, ymin=soc[2], ymax=soc[3], fill="grey90") +
    annotate("rect", xmin=0, xmax=1, ymin=soc[3], ymax=soc[4], fill="grey87") +
    annotate("rect", xmin=0, xmax=1, ymin=soc[4], ymax=soc[5], fill="grey83") +
    
    # Lv 1
    geom_bar(
      data = xx %>% filter(level == 1),
      aes(width = leafCount_norm, y = `median_proportion`), # / leafCount_norm),
      color = "grey20",  stat = "identity"
    ) +
    geom_errorbar(
      data = xx %>% filter(level == 1), # %>% mutate(x = ifelse(`Cell type category`=="immune_cell", 0.5, x)),
      aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[1:2])), # / leafCount_norm),
      color = "grey20",  stat = "identity"
    ) +
    geom_text(
      data =
        xx %>%
        filter(level == 1) %>%
        # mutate(x = ifelse(`Cell type category`=="immune_cell", 0.5, x)) %>%
        # mutate(angle = ifelse(`Cell type category`=="immune_cell", -0, angle)) %>%
        mutate(`formatted` = sub("\\s+$", "", `formatted`)),
      aes(label=`formatted`, y = mean(soc[1:2]), angle= angle  ) ,size =size_geom_text ) +
    #scale_x_continuous(labels = xx %>% filter(level == 1) %>% pull(`formatted`), breaks = xx %>% filter(level == 1) %>% pull(leafCount_norm_cum) - 0.5 * xx %>% filter(level == 1) %>% pull(leafCount_norm)) +
    
    # Lv 2
    geom_bar(
      data = xx %>% filter(level == 2),
      aes(width = leafCount_norm, y = `median_proportion` ), # / leafCount_norm),
      color = "grey20",  stat = "identity"
    ) +
    geom_errorbar(
      data = xx %>% filter(level == 2),
      aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[2:3])), # / leafCount_norm),
      color = "grey20",  stat = "identity",linetype="dotted"
    ) +
    geom_text(
      data =
        xx %>%
        filter(level == 2) %>%
        mutate(`formatted` = sub("\\s+$", "", `formatted`)),
      aes(label=`formatted`, y = mean(soc[2:3]), angle= angle) ,size =size_geom_text) +
    #scale_x_continuous(labels = xx %>% filter(level == 2) %>% pull(`formatted`), breaks = xx %>% filter(level == 2) %>% pull(leafCount_norm_cum) - 0.5 * xx %>% filter(level == 2) %>% pull(leafCount_norm)) +
    
    # Lv 3
    geom_bar(
      data = xx %>% filter(level == 3),
      aes(width = leafCount_norm, y = `median_proportion` ), # / leafCount_norm),
      color = "grey20", stat = "identity"
    )  +
    {
      # Make plotting robust if no level 3 cell types were detected
      switch(
        (!  xx %>% filter(level == 3) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
        
        # If there are cell types
        geom_errorbar(
          data = xx %>% filter(level == 3) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)),
          aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[3:4])), # / leafCount_norm),
          color = "grey20",  stat = "identity",linetype="dotted"
        ) ,
        
        # If there are NOT cell types
        geom_errorbar(ymin=0, ymax=0)
      )
    } +
    {
      # Make plotting robust if no level 3 cell types were detected
      switch(
        (!  xx %>% filter(level == 3) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
        
        # If there are cell types
        
        geom_text(
          data =
            xx %>%
            filter(level == 3)  %>%
            filter(
              `median_proportion`>prop_filter |
                !is.na(Estimate)
            ) %>%
            mutate(`formatted` = sub("\\s+$", "", `formatted`)),
          aes(label=`formatted`, y = mean(soc[3:4]) , angle= angle) ,size =size_geom_text
        ),
        
        # If there are NOT cell types
        geom_errorbar(ymin=0, ymax=0)
      )
    } +
    
    
    # Lv 4
    geom_bar(
      data = xx %>% filter(level == 4),
      aes(width = leafCount_norm, y = `median_proportion` ), # / leafCount_norm),
      color = "grey20", stat = "identity"
    )  +
    {
      # Make plotting robust if no level 4 cell types were detected
      switch(
        (!  xx %>% filter(level == 4) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
        
        # If there are cell types
        geom_errorbar(
          data = xx %>% filter(level == 4) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)),
          aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[4:5])), # / leafCount_norm),
          color = "grey20",  stat = "identity",linetype="dotted"
        ) ,
        
        # If there are NOT cell types
        geom_errorbar(ymin=0, ymax=0)
      )
    } +
    {
      # Make plotting robust if no level 4 cell types were detected
      switch(
        (!  xx %>% filter(level == 4) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
        
        # If there are cell types
        
        geom_text(
          data =
            xx %>%
            filter(level == 4)  %>%
            filter(
              `median_proportion`>prop_filter |
                !is.na(Estimate)
            ) %>%
            mutate(`formatted` = sub("\\s+$", "", `formatted`)),
          aes(label=`formatted`, y = mean(soc[4:5]) , angle= angle) ,size =size_geom_text
        ),
        
        # If there are NOT cell types
        geom_errorbar(ymin=0, ymax=0)
      )
    } +
    
    {
      # Case if none is significant
      switch(
        (!length(na.omit(xx$Estimate))>0) + 1,
        
        # If there are significant
        scale_fill_distiller(
          palette = "Spectral",
          na.value = 'white',
          direction = fill_direction, name = "Trend",
          limits=c(-DTC_scale, DTC_scale)
        ) ,
        
        # If there are not
        scale_fill_manual(values= c("Non significant" = "white" ), guide=FALSE)
      )
      
    } +
    {
      # Case if none is significant
      switch(
        (!length(na.omit(xx$Estimate))>0) + 1,
        
        # If there are significant
        guides(
          fill = guide_colorbar(
            label.position = "left",
            title.position = "left",
            title.hjust = 0.5,
            override.aes=list(fill=NA),
            ticks.colour = "black",
            barwidth = barwidth,
            barheight = barheight
          )
        ) ,
        
        # If there are not
        scale_fill_manual(values= c("Non significant" = "white" ), guide=FALSE)
      )
      
    } +
    scale_y_continuous( breaks=my_breaks, trans = cusotm_root_trans(), labels = my_rescale ) +
    scale_size(range = c(0.3, 0.8), guide=FALSE) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x.top = element_blank(),
      axis.title.x=element_blank(),
      axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), hjust = 0.65, vjust = 1),
      panel.border = element_blank(),
      axis.line.y = element_line(),
      panel.grid  = element_blank(),
      legend.position=c(0,0.05),
      legend.justification=c(legend_justification, 0),
      legend.title=element_text(angle = 90),
      legend.background = element_rect(colour = "transparent", fill = alpha("red", 0))
    )	+ ylab("Cell type proportion") +
    coord_polar(theta = "x")
}


#' plot_polar
#' 
#' @description Create polar plot of results
#' @rdname plot_polar
#'
#'
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @import data.tree
#' @importFrom ape as.phylo
#' @importFrom scales trans_new
#' @importFrom purrr map_df
#'
#'
#' @param .data ARMET-tc object
#' @param size_geom_text A double
#' @param my_breaks An integer
#' @param prop_filter A double
#' @param barwidth A double
#' @param barheight A double
#' @param legend_justification A double
#' @param fill_direction An integer
#'
#' @return a ggplot
#'
#' @export
plot_polar = function(	.data, 
                       size_geom_text = 3.5,
                       my_breaks=c(0, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.7, 1),
                       prop_filter = 0.005,
                       barwidth = 0.5,
                       barheight = 2,
                       legend_justification = 0.67,
                       fill_direction = 1){
  
  xx  =    
    .data %>%
    select(-one_of("draws", "draws_cens")) %>%
    calculate_x_for_polar %>%
    
    # Add angle text
    mutate(angle = x * -360) %>%
    mutate(angle = ifelse(angle < - 240 | angle > -120, angle, angle - 180)) %>%
    mutate(angle = ifelse(angle > - 360, angle, angle + 360)) %>%
    #mutate(median_proportion = ifelse(level<max(level), median_proportion + (0.008*level), median_proportion  )) %>%
    filter(level > 0) %>%
    
    # Calculate median proportions
    filter(.value_alpha1 %>% is.na %>% `!`) %>%
    mutate(summarised_proportion =
             proportions %>% map(~ .x %>% summarise(median_proportion = .value %>% median))
    ) %>%
    unnest(summarised_proportion) %>%
    #mutate(median_proportion = median_proportion / max(median_proportion)) %>%
    
    distinct() %>%
    left_join(
      ct_names_polar,
      by=c("Cell type category" = "name")
    )
  
  my_rescale = function(x) { as.character( format(round( x * xx %>% pull(median_proportion) %>% max, 2), nsmall = 2)) }
  cusotm_root_trans = function() scales::trans_new("cusotm_root",function(x) x^(1/4), function(x) x^(4))
  DTC_scale =  xx %>% pull(Estimate) %>% abs() %>% { switch( (!all(is.na(.))) + 1, 1, .) } %>% max(na.rm = T)
  y_text_out = 8.5
  y_text_in = 1.4641
  
  # Size outer circles
  soc = (c(0, 0.2, 0.4, 0.6, 0.8)*0.5 + 1)^4  
  
  (
    xx %>%
    {
      # Case if none is significant
      switch(
        (!length(na.omit(xx$Estimate))>0) + 1,
        
        # If there are significant
        ggplot(data=(.), aes(x = x,fill = Estimate,size = 1/sqrt(level))),
        
        # If there are not
        ggplot(data=(.),aes(x = x,fill = "Non significant",size = 1/sqrt(level)))
      )
    }	+
    annotate("rect", xmin=0, xmax=1, ymin=soc[1], ymax= soc[2], fill="grey95") +
    annotate("rect", xmin=0, xmax=1, ymin=soc[2], ymax=soc[3], fill="grey90") +
    annotate("rect", xmin=0, xmax=1, ymin=soc[3], ymax=soc[4], fill="grey87") +
    annotate("rect", xmin=0, xmax=1, ymin=soc[4], ymax=soc[5], fill="grey83")
  ) %>% 
    
    # Lv 1
    magnolia_layer(xx, 1, prop_filter, size_geom_text) %>% 
    
    # Lv 2
    magnolia_layer(xx, 2, prop_filter, size_geom_text)%>%  
    
    # Lv 3
    magnolia_layer(xx, 3, prop_filter, size_geom_text, linetype="dotted") %>% 
    
    magnolia_layer(xx, 4, prop_filter, size_geom_text, linetype="dotted") +
    
    # geom_bar(
    #   data = xx %>% filter(level == 3),
    #   aes(width = leafCount_norm, y = `median_proportion` ), # / leafCount_norm),
    #   color = "grey20", stat = "identity"
    # )  +
    # {
    #   # Make plotting robust if no level 3 cell types were detected
    #   switch(
    #     (!  xx %>% filter(level == 3) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
    #     
    #     # If there are cell types
    #     geom_errorbar(
    #       data = xx %>% filter(level == 3) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)),
    #       aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[3:4])), # / leafCount_norm),
    #       color = "grey20",  stat = "identity",linetype="dotted"
    #     ) ,
    #     
    #     # If there are NOT cell types
    #     geom_errorbar(ymin=0, ymax=0)
    #   )
    # } +
    # {
    #   # Make plotting robust if no level 3 cell types were detected
    #   switch(
    #     (!  xx %>% filter(level == 3) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
    #     
    #     # If there are cell types
    #     
    #     geom_text(
    #       data =
    #         xx %>%
    #         filter(level == 3)  %>%
    #         filter(
    #           `median_proportion`>prop_filter |
    #             !is.na(Estimate)
    #         ) %>%
    #         mutate(`formatted` = sub("\\s+$", "", `formatted`)),
    #       aes(label=`formatted`, y = mean(soc[3:4]) , angle= angle) ,size =size_geom_text
    #     ),
    #     
    #     # If there are NOT cell types
    #     geom_errorbar(ymin=0, ymax=0)
    #   )
    # } +
    # 
    # 
    # # Lv 4
    # geom_bar(
    #   data = xx %>% filter(level == 4),
    #   aes(width = leafCount_norm, y = `median_proportion` ), # / leafCount_norm),
    #   color = "grey20", stat = "identity"
    # )  +
    # {
    #   # Make plotting robust if no level 4 cell types were detected
    #   switch(
    #     (!  xx %>% filter(level == 4) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
    #     
    #     # If there are cell types
    #     geom_errorbar(
    #       data = xx %>% filter(level == 4) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)),
    #       aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[4:5])), # / leafCount_norm),
    #       color = "grey20",  stat = "identity",linetype="dotted"
    #     ) ,
    #     
    #     # If there are NOT cell types
    #     geom_errorbar(ymin=0, ymax=0)
    #   )
    # } +
    # {
    #   # Make plotting robust if no level 4 cell types were detected
    #   switch(
    #     (!  xx %>% filter(level == 4) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
    #     
    #     # If there are cell types
    #     
    #     geom_text(
    #       data =
    #         xx %>%
    #         filter(level == 4)  %>%
    #         filter(
    #           `median_proportion`>prop_filter |
    #             !is.na(Estimate)
    #         ) %>%
    #         mutate(`formatted` = sub("\\s+$", "", `formatted`)),
    #       aes(label=`formatted`, y = mean(soc[4:5]) , angle= angle) ,size =size_geom_text
    #     ),
    #     
    #     # If there are NOT cell types
    #     geom_errorbar(ymin=0, ymax=0)
    #   )
    # } +
    
    {
      # Case if none is significant
      switch(
        (!length(na.omit(xx$Estimate))>0) + 1,
        
        # If there are significant
        scale_fill_distiller(
          palette = "Spectral",
          na.value = 'white',
          direction = fill_direction, name = "Trend",
          limits=c(-DTC_scale, DTC_scale)
        ) ,
        
        # If there are not
        scale_fill_manual(values= c("Non significant" = "white" ), guide=FALSE)
      )
      
    } +
    {
      # Case if none is significant
      switch(
        (!length(na.omit(xx$Estimate))>0) + 1,
        
        # If there are significant
        guides(
          fill = guide_colorbar(
            label.position = "left",
            title.position = "left",
            title.hjust = 0.5,
            override.aes=list(fill=NA),
            ticks.colour = "black",
            barwidth = barwidth,
            barheight = barheight
          )
        ) ,
        
        # If there are not
        scale_fill_manual(values= c("Non significant" = "white" ), guide=FALSE)
      )
      
    } +
    scale_y_continuous( breaks=my_breaks, trans = cusotm_root_trans(), labels = my_rescale ) +
    scale_size(range = c(0.3, 0.8), guide=FALSE) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x.top = element_blank(),
      axis.title.x=element_blank(),
      axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), hjust = 0.65, vjust = 1),
      panel.border = element_blank(),
      axis.line.y = element_line(),
      panel.grid  = element_blank(),
      legend.position=c(0,0.05),
      legend.justification=c(legend_justification, 0),
      legend.title=element_text(angle = 90),
      legend.background = element_rect(colour = "transparent", fill = alpha("red", 0))
    )	+ ylab("Cell type proportion") +
    coord_polar(theta = "x")
}

