
bargraph_jitter <- function(data, stats = NULL, x, y, group = NULL,
                                  ylab = "y",  errorbar = "SEM", dotsize = 1,
                            dot_fill = "white", fill_colors = NULL,
                                  title = element_blank(), subtitle = element_blank()){
  
  
  # enquote variables
  x = enquo(x)
  y = enquo(y)
  group = enquo(group)
  
  if(rlang::quo_is_null(group)){
    
    stop("group variable is missing!")
    
  }
  
  plot_out <- data %>%
    ggplot(aes(x = !!x, y = !!y , group = !!group)) +
    # Plot bars
    geom_bar(stat ="summary",
             fun = "mean",
             aes(fill = !!group)) +
    
    {if(errorbar == "SD"){ 
      # Plot errorbars as Mean+- sd
      stat_summary(fun.data="mean_sdl",           # compute mean+-sd
                   fun.args = list(mult = 1),     # set multiplication of sd as 1
                   geom="errorbar", color="black", width=0.2)
      
    } else if(errorbar == "SEM"){
      # Plot errorbars as Mean+- sem
      stat_summary(fun.data="mean_se",           # compute mean+-sem
                   geom="errorbar", color="black", width=0.2)
    }} +
    
    # Add jitter
    geom_jitter(width = 0.1,
                shape = 21,
                colour = "black",
                fill = dot_fill,
                size = dotsize) +
    
    # Add labels
    labs(x = element_blank(), y = ylab, title = title,
         subtitle = subtitle) +
    
    # Adjust theme
    theme(
      plot.title = element_text(hjust = 0.5),     # Center title
      plot.subtitle = element_text(hjust = 0.5)   # Center subtitle
    ) +
    
    # Change color of bars
    {if(!is.null(fill_colors)){
      scale_fill_manual(values = fill_colors)
    }} +
    
    # Adjust axes
    #scale_x_discrete(expand = c(0, 0.5)) +
    scale_y_continuous(expand = c(0,0))
  
  # Add p-values if there are any significant differences
  
  if( !is.null(stats)){
    if (sum(stats[,"p.adj.signif"] != "ns") > 1){ 
      plot_out <- plot_out + stat_pvalue_manual(
        stats,  label = "p.adj.signif", tip.length = 0,
        hide.ns = TRUE)}
    
    # Adjust the y scale so the points or signif. values are not cut
    ymax <- max(ggplot_build(
      plot_out)$layout$panel_params[[1]]$y.range) # Extract the maximum y-value from the plot
    #ymaxAdj <- ceiling((ymax * 2))/2   # Round to the nearset 0.5 upwards
    ymaxAdj <- ceiling(ymax)   # Round upwards
    
    plot_out <- plot_out + scale_y_continuous(expand = c(0,0), limits = c(0, ymaxAdj))
  }
    plot_out + coord_cartesian(clip = "off")
  
}
