


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
#' @importFrom purrr when
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
    calculate_x_for_polar %>%
    distinct() %>% 
    # Add angle text
    mutate(angle = x * -360) %>%
    mutate(angle = ifelse(angle < - 240 | angle > -120, angle, angle - 180)) %>%
    mutate(angle = ifelse(angle > - 360, angle, angle + 360)) %>%
    filter(level > 0) 
  
  my_rescale = function(x) { as.character( format(round( x * xx %>% pull(y) %>% max, 2), nsmall = 2)) }
  cusotm_root_trans = function() scales::trans_new("cusotm_root",function(x) x^(1/4), function(x) x^(4))
  DTC_scale =  xx %>% pull(value) %>% abs() %>% { switch( (!all(is.na(.))) + 1, 1, .) } %>% max(na.rm = T)
  y_text_out = 8.5
  y_text_in = 1.4641
  
  # magnolia_input = xx %>% replace_na(list(ancestor = "tissue")) %>% select(child = `Cell type category`, parent = ancestor, do_color = significant, y = median_proportion, value=Estimate)  %>%
  # left_join(
  #   ct_names_polar, by=c(child = "name")
  # ) %>% select(-child) %>% rename(child = formatted) %>% select(child, parent, everything()) %>% left_join(
  #   ct_names_polar, by=c(parent = "name")
  # ) %>% select(-parent) %>% rename(parent = formatted) %>% select(child, parent, everything()) %>% replace_na(list(parent = "tissue")) 
  # Browse[2]> save(magnolia_input, file="data/magnolia_input.rda", compress = "xz")
  # FromDataFrameNetwork(xx %>% select(ancestor, `child`) %>% replace_na(list(ancestor = "tissue")))
  
  
  # Size outer circles
  soc = (c(0, 0.2, 0.4, 0.6, 0.8)*0.5 + 1)^4  
  
  (
    xx %>%
    {
      # Case if none is significant
      switch(
        (!length(na.omit(xx$value))>0) + 1,
        
        # If there are significant
        ggplot(data=(.), aes(x = x,fill = value,size = 1/sqrt(level))),
        
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
    magnolia_layer(xx, 1, prop_filter, size_geom_text, soc = soc) %>% 
    
    # Lv 2
    magnolia_layer(xx, 2, prop_filter, size_geom_text, soc = soc)%>%  
    
    # Lv 3
    magnolia_layer(xx, 3, prop_filter, size_geom_text, linetype="dotted", soc = soc) %>% 
    
    # Lv 4
    magnolia_layer(xx, 4, prop_filter, size_geom_text, linetype="dotted", soc = soc) +
    

    
    {
      # Case if none is significant
      switch(
        (!length(na.omit(xx$value))>0) + 1,
        
        # If there are significant
        scale_fill_distiller(
          palette = "Spectral",
          na.value = 'white',
          direction = fill_direction, 
          name = "Trend",
          limits=c(-DTC_scale, DTC_scale)
        ) ,
        
        # If there are not
        scale_fill_manual(values= c("Non significant" = "white" ), guide=FALSE)
      )
      
    } +
    {
      # Case if none is significant
      switch(
        (!length(na.omit(xx$value))>0) + 1,
        
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

