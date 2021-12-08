get_ChIP_enrichment = function(tib, data_col='median') {
  chip_data = tib %>% dplyr::filter(ChIP==1) %>% .[,data_col]
  chip_data = chip_data[[1]][[1]]
  input_data = tib %>% dplyr::filter(ChIP==0) %>% .[,data_col]
  input_data = input_data[[1]][[1]]
  
  chip_enrich = chip_data/input_data
  return(chip_enrich)
}


convert_to_robustz = function(input_vec){
  finite <- is.finite(input_vec)
  this_median <- median(input_vec[finite])
  mad <- median(abs(input_vec[finite] - this_median))
  output <- (input_vec - this_median)/(mad*1.4826)
  output
}

getPercent = function(x) {
  x/sum(x)*100
}

getCPM = function(x) {
  x/sum(x)*1e6
}

getLog2CPM = function(x) {
  log2((x+1e-5)/sum(x)*1e6)
}

applyGetCPM = function(mat) {
  apply(mat, MARGIN=2, FUN=getCPM)
}

applyGetLog2CPM = function(mat) {
  apply(mat, MARGIN=2, FUN=getLog2CPM)
}

#plotLocus----
plotLocus = function(CDSdf, ChIPdf, plotStart, plotEnd, 
                     lineAlpha=1, lineSize=1.25, ylabel="log2 (fold ChIP enrichment)",
                     ylims="detect", useLog=TRUE, plotFeatures=TRUE) {
  
  plotCDSs = CDSdf %>% filter(end > plotStart, start < plotEnd)
  if (useLog) {
    plotChIP = ChIPdf %>% filter(position<plotEnd, position>plotStart)
  } else {
    plotChIP = ChIPdf %>% filter(position<plotEnd, position>plotStart) %>% mutate(log2Enrich=2**log2Enrich)
  }
  
  
  geneTop = min(plotChIP$log2Enrich) - 0.5
  geneBottom = geneTop-0.25
  geneText = geneTop + 0.25
  
  plot = ggplot()
  
  if (plotFeatures) {
    plot = plot + 
      geom_rect(data=plotCDSs, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneBottom, ymax=geneTop, fill=strand),
                color="black") +
                scale_fill_manual(values=strand_colors) +
      geom_text(data=plotCDSs, aes(x=midpoint/1e6, y=geneText, label=paste0("italic(",gene,")")), parse=T)
  }
  plot = plot + 
    geom_line(data=plotChIP, aes(x=position/1e6, y=log2Enrich, color=Description), size=lineSize, alpha=lineAlpha) +
    theme_classic() +
    theme(text = element_text(size=15),
          axis.line = element_line(size=1.5),
          axis.text = element_text(size=12, color="black"),
          axis.ticks = element_line(color="black")) +
    labs(y=ylabel, x="Genome position (Mb)") +
    coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6)) +
  
  if (!(ylims == "detect")) {
    plot = plot + coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6), ylim=ylims)
  }
  
  if (plotPurBox) {
    #print(plotBox$start)
    if(length(plotBox$start)==0) {
      plot=plot
    }
    if(!(length(plotBox$start)==0)) {
      plot = plot +
        geom_rect(data=plotBox, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneTop, ymax=Inf), alpha=0.4)
    }
  }
  
  return(plot)
}

#plotEnrich----
plotEnrich = function(CDSdf, ChIPdf, plotStart, plotEnd, 
                     lineAlpha=1, lineSize=1.25, ylabel="log2 (fold ChIP enrichment)",
                     ylims="detect", useLog=TRUE, plotFeatures=TRUE, plotPurBox=FALSE, plotError=FALSE) {
  
  plotCDSs = CDSdf %>% filter(end > plotStart, start < plotEnd)
  if (useLog) {
    plotChIP = ChIPdf %>% filter(position<plotEnd, position>plotStart)
  } else {
    plotChIP = ChIPdf %>% filter(position<plotEnd, position>plotStart) %>% mutate(mean_log2_enrich=2**mean_log2_enrich)
  }
  
  #print(plotChIP %>% .$mean_log2_enrich)
  min_chip = min(plotChIP %>% .$mean_log2_enrich)
  max_chip = max(plotChIP %>% .$mean_log2_enrich)
  chip_range = max_chip - min_chip
  #geneTop = min(plotChIP$mean_log2_enrich) - 0.5
  geneTop = min_chip - 0.05*chip_range
  geneBottom = geneTop - 0.075*chip_range
  geneText = geneBottom - 0.05 * chip_range
  plotBox = purbox_df %>% dplyr::filter(end > plotStart, start < plotEnd)
  
  plot = ggplot()
  
  if (plotFeatures) {
    plot = plot + 
      geom_rect(data=plotCDSs, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneBottom, ymax=geneTop, fill=strand),
                color="black", show.legend=FALSE) +
                scale_fill_manual(values=strand_colors) +
      geom_text(data=plotCDSs, aes(x=midpoint/1e6, y=geneText, label=paste0("italic(",gene,")")), check_overlap=T, parse=T)
  }
  if (plotError) {
    plot = plot +
      new_scale_fill() +
      geom_ribbon(data=plotChIP, aes(x=position/1e6, ymin=mean_log2_enrich-sd_log2_enrich, ymax=mean_log2_enrich+sd_log2_enrich, fill=group), alpha=0.4) +
      scale_fill_manual(values=sample_colors)
  }
  plot = plot + 
    geom_line(data=plotChIP, aes(x=position/1e6, y=mean_log2_enrich, color=group), size=lineSize, alpha=lineAlpha) +
    scale_color_manual(values=sample_colors) +
    theme_classic() +
    theme(text = element_text(size=15),
          axis.line = element_line(size=1.5),
          axis.text = element_text(size=12, color="black"),
          axis.ticks = element_line(color="black")) +
    labs(y=ylabel, x="Genome position (Mb)") +
    coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    #aes(group=rev(group)) +
    
    if (!(ylims == "detect")) {
      plot = plot + coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6), ylim=ylims)
    }
  
  if (plotPurBox) {
    #print(plotBox$start)
    if(length(plotBox$start)==0) {
      plot=plot
    }
    if(!(length(plotBox$start)==0)) {
      plot = plot +
        geom_rect(data=plotBox, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneTop, ymax=Inf), alpha=0.4)
    }
  }
  return(plot)
}

#plotEnrichFigure: plot linear enrichment for figures----
plotEnrichFigure = function(CDSdf, ChIPdf, plotStart, plotEnd, 
                      lineAlpha=1, lineSize=0.3, ylabel="ChIP enrichment",
                      ylims="detect", useLog=TRUE, plotFeatures=TRUE, plotPurBox=FALSE, plotError=FALSE) {
  
  plotCDSs = CDSdf %>% filter(end > plotStart, start < plotEnd)
  if (useLog) {
    plotChIP = ChIPdf %>% filter(position<plotEnd, position>plotStart)
  } else {
    plotChIP = ChIPdf %>% filter(position<plotEnd, position>plotStart) %>% mutate(mean_enrich=2**mean_enrich)
  }
  
  #print(plotChIP %>% .$mean_log2_enrich)
  min_chip = min(plotChIP %>% .$mean_enrich)
  max_chip = max(plotChIP %>% .$mean_enrich)
  chip_range = max_chip - min_chip
  #print(chip_range)
  #geneTop = min(plotChIP$mean_log2_enrich) - 0.5
  geneTop = min_chip - 0.05*chip_range
  geneBottom = geneTop - 0.075*chip_range
  geneText = geneBottom - 0.05 * chip_range
  plotBox = purbox_df %>% dplyr::filter(end > plotStart, start < plotEnd)
  
  plot = ggplot()
  
  if (plotFeatures) {
    plot = plot + 
      geom_rect(data=plotCDSs, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneBottom, ymax=geneTop, fill=strand),
                color="black", show.legend=FALSE, size=0.15) +
      scale_fill_manual(values=strand_colors) +
      geom_text(data=plotCDSs, aes(x=midpoint/1e6, y=geneText, label=paste0("italic(",gene,")")), 
                check_overlap=T, parse=T, size=1.65)
  }
  if (plotError) {
    plot = plot +
      new_scale_fill() +
      geom_ribbon(data=plotChIP, aes(x=position/1e6, ymin=mean_enrich-sd_enrich, ymax=mean_enrich+sd_enrich, fill=group), alpha=0.4) +
      scale_fill_manual(values=sample_colors)
  }
  plot = plot + 
    geom_line(data=plotChIP, aes(x=position/1e6, y=mean_enrich, color=group), size=lineSize, alpha=lineAlpha) +
    scale_color_manual(values=sample_colors) +
    theme_classic() +
    theme(text = element_text(size=8),
          axis.line = element_line(size=0.2),
          axis.text = element_text(size=6, color="black"),
          axis.ticks = element_line(color="black", size=0.2),
          legend.position='none') +
    labs(y=ylabel, x="Genome position (Mb)") +
    coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
    #aes(group=rev(group)) +
    
    if (!(ylims == "detect")) {
      plot = plot + coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6), ylim=ylims)
    }
  
  if (plotPurBox) {
    #print(plotBox$start)
    if(length(plotBox$start)==0) {
      plot=plot
    }
    if(!(length(plotBox$start)==0)) {
      plot = plot +
        geom_rect(data=plotBox, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneTop, ymax=Inf), alpha=0.4)
    }
  }
  return(plot)
}

#plotLocusCPM----
plotLocusCPM = function(CDSdf, ChIPdf, plotStart, plotEnd, 
                     lineAlpha=1, lineSize=1.25, ylabel="Reads per millions mapped",
                     ylims="detect", plotFeatures=TRUE, arrowLength=1.5, plotPurBox=FALSE) {
  
  plotCDSs = CDSdf %>% dplyr::filter(end > plotStart, start < plotEnd)
  plotChIP = ChIPdf %>% dplyr::filter(position<plotEnd, position>plotStart)
  plotBox = purbox_df %>% dplyr::filter(end > plotStart, start < plotEnd)
  
  min_chip = min(plotChIP$cpm)
  max_chip = max(plotChIP$cpm)
  chip_range = max_chip - min_chip
  geneTop = min_chip - 0.05*chip_range
  geneBottom = geneTop - 0.075*chip_range
  geneText = geneBottom - 0.05 * chip_range
  
  plot = ggplot()
  
  if (plotFeatures) {
    # print(plotCDSs %>% head)
    plot = plot + 
      geom_rect(data=plotCDSs, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneBottom, ymax=geneTop, fill=strand), 
                color="black") +
                scale_fill_manual(values=strand_colors) +
      geom_text(data=plotCDSs, aes(x=midpoint/1e6, y=geneText, label=paste0("italic(",gene,")")), parse=T)
  }
  # print(plotChIP %>% head)
  plot = plot + 
    geom_line(data=plotChIP, aes(x=position/1e6, y=cpm, color=Description), size=lineSize, alpha=lineAlpha) +
    scale_color_manual(values=sample_colors) +
    theme_classic() +
    theme(text = element_text(size=18),
          axis.line = element_line(size=2),
          axis.text = element_text(size=14, color="black"),
          axis.ticks = element_line(color="black")) +
    labs(y=ylabel, x="Genome position (Mb)") +
    coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    #aes(group=rev(Description)) +
  
  if (!(ylims == "detect")) {
    plot = plot + coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6), ylim=ylims)
  }
  
  if (plotPurBox) {
    #print(plotBox$start)
    if(length(plotBox$start)==0) {
      plot=plot
    }
    if(!(length(plotBox$start)==0)) {
      plot = plot +
        geom_rect(data=plotBox, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneTop, ymax=Inf), alpha=0.4)
    }
  }
  
  return(plot)
}

#plotEnrichLinear: plot linear enrichment----
plotEnrichLinear = function(CDSdf, ChIPdf, plotStart, plotEnd, 
                            lineAlpha=1, lineSize=1, ylabel="ChIP enrichment",
                            ylims='detect', useLog=TRUE, plotFeatures=TRUE, plotPurBox=FALSE, plotError=FALSE) {
  
  plotCDSs = CDSdf %>% filter(end > plotStart, start < plotEnd)
  if (useLog) {
    plotChIP = ChIPdf %>% filter(position<plotEnd, position>plotStart)
  } else {
    plotChIP = ChIPdf %>% filter(position<plotEnd, position>plotStart) %>% mutate(mean_enrich=2**mean_log2_enrich)
  }
  
  #print(plotStart)
  #print(plotChIP %>% .$mean_log2_enrich)
  min_chip = min(plotChIP %>% .$mean_enrich)
  max_chip = max(plotChIP %>% .$mean_enrich)
  chip_range = max_chip - min_chip
  #geneTop = min(plotChIP$mean_log2_enrich) - 0.5
  geneTop = min_chip - 0.05*chip_range
  geneBottom = geneTop - 0.075*chip_range
  geneText = geneBottom - 0.05 * chip_range
  plotBox = purbox_df %>% dplyr::filter(end > plotStart, start < plotEnd)
  
  plot = ggplot()
  
  if (plotFeatures) {
    plot = plot + 
      geom_rect(data=plotCDSs, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneBottom, ymax=geneTop, fill=strand),
                color="black", show.legend=FALSE, size=0.25) +
      scale_fill_manual(values=strand_colors) +
      geom_text(data=plotCDSs, aes(x=midpoint/1e6, y=geneText, label=paste0("italic(",gene,")")), 
                check_overlap=T, parse=T, size=4.5)
  }
  if (plotError) {
    plot = plot +
      new_scale_fill() +
      geom_ribbon(data=plotChIP, aes(x=position/1e6, ymin=mean_enrich-sd_enrich, ymax=mean_enrich+sd_enrich, fill=group), alpha=0.4) +
      scale_fill_manual(values=sample_colors)
  }
  plot = plot + 
    geom_line(data=plotChIP, aes(x=position/1e6, y=mean_enrich, color=group), size=lineSize, alpha=lineAlpha) +
    scale_color_manual(values=sample_colors) +
    theme_classic() +
    theme(text = element_text(size=14),
          axis.line = element_line(size=1),
          axis.text = element_text(size=12, color="black"),
          axis.ticks = element_line(color="black", size=1),
          legend.position='none') +
    labs(y=ylabel, x="Genome position (Mb)") +
    coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
    #aes(group=rev(group)) +
    
    if (!(ylims == "detect")) {
      plot = plot + coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6), ylim=ylims)
    }
  
  if (plotPurBox) {
    #print(plotBox$start)
    if(length(plotBox$start)==0) {
      plot=plot
    }
    if(!(length(plotBox$start)==0)) {
      plot = plot +
        geom_rect(data=plotBox, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneTop, ymax=Inf), alpha=0.4)
    }
  }
  return(plot)
}

#plotLocusCPM----
plotLocusCPM = function(CDSdf, ChIPdf, plotStart, plotEnd, 
                        lineAlpha=1, lineSize=1.25, ylabel="Reads per millions mapped",
                        ylims="detect", plotFeatures=TRUE, arrowLength=1.5, plotPurBox=FALSE) {
  
  plotCDSs = CDSdf %>% dplyr::filter(end > plotStart, start < plotEnd)
  plotChIP = ChIPdf %>% dplyr::filter(position<plotEnd, position>plotStart)
  plotBox = purbox_df %>% dplyr::filter(end > plotStart, start < plotEnd)
  
  min_chip = min(plotChIP$cpm)
  max_chip = max(plotChIP$cpm)
  chip_range = max_chip - min_chip
  geneTop = min_chip - 0.05*chip_range
  geneBottom = geneTop - 0.075*chip_range
  geneText = geneBottom - 0.05 * chip_range
  
  plot = ggplot()
  
  if (plotFeatures) {
    # print(plotCDSs %>% head)
    plot = plot + 
      geom_rect(data=plotCDSs, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneBottom, ymax=geneTop, fill=strand), 
                color="black") +
      scale_fill_manual(values=strand_colors) +
      geom_text(data=plotCDSs, aes(x=midpoint/1e6, y=geneText, label=paste0("italic(",gene,")")), parse=T)
  }
  # print(plotChIP %>% head)
  plot = plot + 
    geom_line(data=plotChIP, aes(x=position/1e6, y=cpm, color=Description), size=lineSize, alpha=lineAlpha) +
    scale_color_manual(values=sample_colors) +
    theme_classic() +
    theme(text = element_text(size=18),
          axis.line = element_line(size=2),
          axis.text = element_text(size=14, color="black"),
          axis.ticks = element_line(color="black")) +
    labs(y=ylabel, x="Genome position (Mb)") +
    coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    #aes(group=rev(Description)) +
    
    if (!(ylims == "detect")) {
      plot = plot + coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6), ylim=ylims)
    }
  
  if (plotPurBox) {
    #print(plotBox$start)
    if(length(plotBox$start)==0) {
      plot=plot
    }
    if(!(length(plotBox$start)==0)) {
      plot = plot +
        geom_rect(data=plotBox, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneTop, ymax=Inf), alpha=0.4)
    }
  }
  
  return(plot)
}

#plotLocusCPM_avg----
plotLocusCPM_avg = function(CDSdf, ChIPdf, plotStart, plotEnd, 
                        lineAlpha=1, lineSize=1.25, ylabel=expression(atop("Reads per", paste("millions mapped"))),
                        ylims='detect', plotFeatures=TRUE, arrowLength=1.5, plotPurBox=FALSE, plotError=FALSE) {
  
  plotCDSs = CDSdf %>% dplyr::filter(end > plotStart, start < plotEnd)
  plotChIP = ChIPdf %>% dplyr::filter(position<plotEnd, position>plotStart)
  plotBox = purbox_df %>% dplyr::filter(end > plotStart, start < plotEnd)
  
  #print(plotChIP %>% .$mean_cpm)
  min_chip = min(plotChIP %>% .$mean_cpm)
  max_chip = max(plotChIP %>% .$mean_cpm)
  chip_range = max_chip - min_chip
  geneTop = min_chip - 0.05*chip_range
  geneBottom = geneTop - 0.075*chip_range
  geneText = geneBottom - 0.05 * chip_range
  
  plot = ggplot()
  
  if (plotFeatures) {
    # print(plotCDSs %>% head)
    plot = plot + 
      geom_rect(data=plotCDSs, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneBottom, ymax=geneTop, fill=strand), 
                color="black") +
      scale_fill_manual(values=strand_colors) +
      geom_text(data=plotCDSs, aes(x=midpoint/1e6, y=geneText, label=paste0("italic(",gene,")")), parse=T)
  }
  
  if (plotError) {
    plot = plot +
      new_scale_fill() +
      geom_ribbon(data=plotChIP, aes(x=position/1e6, ymin=mean_cpm-sd_cpm, ymax=mean_cpm+sd_cpm, fill=Description), alpha=0.4) +
      scale_fill_manual(values=sample_colors)
  }
  
  #print(plotChIP %>% head)
  plot = plot + 
    geom_line(data=plotChIP, aes(x=position/1e6, y=mean_cpm, color=Description), size=lineSize, alpha=lineAlpha) +
    scale_color_manual(values=sample_colors) +
    theme_classic() +
    theme(text = element_text(size=8),
          axis.line = element_line(size=0.3),
          axis.text = element_text(size=5.5, color="black"),
          axis.ticks = element_line(color="black", size=0.3),
          legend.position = 'none') +
    labs(y=ylabel, x="Genome position (Mb)") +
    coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
    #aes(group=rev(Description)) +
    #ylim(0,180)
    
    if (!(ylims == "detect")) {
      plot = plot + coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6), ylim=ylims)
    }
  
  if (plotPurBox) {
    #print(plotBox$start)
    if(length(plotBox$start)==0) {
      plot=plot
    }
    if(!(length(plotBox$start)==0)) {
      plot = plot +
        geom_rect(data=plotBox, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneTop, ymax=Inf), alpha=0.4)
    }
  }
  
  return(plot)
}


#plotEffect----
plotEffect = function(CDSdf, ChIPdf, plotStart, plotEnd, 
                      lineAlpha=1, lineSize=1.25, ylabel="(p)ppGpp effect (AU)",
                      ylims="detect", useLog=TRUE, plotFeatures=TRUE, plotPurBox=FALSE, plotError=FALSE) {
  
  plotCDSs = CDSdf %>% filter(end > plotStart, start < plotEnd)
  if (useLog) {
    plotChIP = ChIPdf %>% filter(position<plotEnd, position>plotStart)
  } else {
    plotChIP = ChIPdf %>% filter(position<plotEnd, position>plotStart) %>% mutate(mean=2**mean, sd=2**sd)
  }
  
  #print(plotChIP %>% .$mean_log2_enrich)
  min_chip = min(plotChIP %>% .$mean)
  max_chip = max(plotChIP %>% .$mean)
  chip_range = max_chip - min_chip
  #geneTop = min(plotChIP$mean_log2_enrich) - 0.5
  geneTop = min_chip - 0.2*chip_range
  geneBottom = geneTop - 0.225*chip_range
  geneText = geneBottom - 0.1 * chip_range
  plotBox = purbox_df %>% dplyr::filter(end > plotStart, start < plotEnd)
  
  plot = ggplot()
  
  if (plotFeatures) {
    plot = plot + 
      geom_rect(data=plotCDSs, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneBottom, ymax=geneTop, fill=strand),
                color="black", show.legend=FALSE) +
      scale_fill_manual(values=strand_colors) +
      geom_text(data=plotCDSs, aes(x=midpoint/1e6, y=geneText, label=paste0("italic(",gene,")")), check_overlap=T, parse=T)
  }
  if (plotError) {
    plot = plot +
      new_scale_fill() +
      geom_ribbon(data=plotChIP, aes(x=position/1e6, ymin=mean-sd, ymax=mean+sd), alpha=0.4,
                  fill=(values=line_color)) 
      #scale_fill_manual(values=line_color)
  }
  plot = plot + 
    geom_line(data=plotChIP, aes(x=position/1e6, y=mean), size=lineSize, alpha=lineAlpha,
              color=(values=line_color)) +
    #scale_color_manual(values=line_color) +
    theme_classic() +
    theme(text = element_text(size=15),
          axis.line = element_line(size=1.5),
          axis.text = element_text(size=12, color="black"),
          axis.ticks = element_line(color="black")) +
    labs(y=ylabel, x="Genome position (Mb)") +
    coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    #aes(group=rev(group))
    
    if (!(ylims == "detect")) {
      plot = plot + coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6), ylim=ylims)
    }
  
  if (plotPurBox) {
    #print(plotBox$start)
    if(length(plotBox$start)==0) {
      plot=plot
    }
    if(!(length(plotBox$start)==0)) {
    plot = plot +
      geom_rect(data=plotBox, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneTop, ymax=Inf), alpha=0.4)
    }
  }
  return(plot)
}


#plotEffectFigure----
plotEffectFigure = function(CDSdf, ChIPdf, plotStart, plotEnd, 
                            lineAlpha=1, lineSize=0.3, ylabel="(p)ppGpp effect (AU)",
                            ylims="detect", useLog=TRUE, plotFeatures=TRUE, plotPurBox=FALSE, plotError=FALSE) {
  
  plotCDSs = CDSdf %>% filter(end > plotStart, start < plotEnd)
  if (useLog) {
    plotChIP = ChIPdf %>% filter(position<plotEnd, position>plotStart)
  } else {
    plotChIP = ChIPdf %>% filter(position<plotEnd, position>plotStart) %>% mutate(mean=2**mean) %>% mutate(sd=2**sd)
  }
  
  #print(plotChIP %>% .$mean_log2_enrich)
  min_chip = min(plotChIP %>% .$mean)
  max_chip = max(plotChIP %>% .$mean)
  chip_range = max_chip - min_chip
  #print(chip_range)
  #geneTop = min(plotChIP$mean_log2_enrich) - 0.5
  geneTop = min_chip - 0.2*chip_range
  geneBottom = geneTop - 0.225*chip_range
  geneText = geneBottom - 0.1 * chip_range
  plotBox = purbox_df %>% dplyr::filter(end > plotStart, start < plotEnd)
  
  plot = ggplot()
  
  if (plotFeatures) {
    plot = plot + 
      geom_rect(data=plotCDSs, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneBottom, ymax=geneTop, fill=strand),
                color="black", show.legend=FALSE, size=0.15) +
      scale_fill_manual(values=strand_colors) +
      geom_text(data=plotCDSs, aes(x=midpoint/1e6, y=geneText, label=paste0("italic(",gene,")")), 
                check_overlap=T, parse=T, size=1.65)
  }
  if (plotError) {
    plot = plot +
      new_scale_fill() +
      geom_ribbon(data=plotChIP, aes(x=position/1e6, ymin=mean-sd, ymax=mean+sd), alpha=0.4,
                  fill=(values=line_color))
      #scale_fill_manual(values=sample_colors)
  }
  plot = plot + 
    geom_line(data=plotChIP, aes(x=position/1e6, y=mean), size=lineSize, alpha=lineAlpha,
              color=(values=line_color)) +
    #scale_color_manual(values=sample_colors) +
    theme_classic() +
    theme(text = element_text(size=8),
          axis.line = element_line(size=0.2),
          axis.text = element_text(size=6, color="black"),
          axis.ticks = element_line(color="black", size=0.2),
          legend.position='none') +
    labs(y=ylabel, x="Genome position (Mb)") +
    coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  #aes(group=rev(group)) +
  
  if (!(ylims == "detect")) {
    plot = plot + coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6), ylim=ylims)
  }
  
  if (plotPurBox) {
    #print(plotBox$start)
    if(length(plotBox$start)==0) {
      plot=plot
    }
    if(!(length(plotBox$start)==0)) {
      plot = plot +
        geom_rect(data=plotBox, aes(xmin=start/1e6, xmax=end/1e6, ymin=geneTop, ymax=Inf), alpha=0.4)
    }
  }
  return(plot)
}


#Other stuff ----

getRPPNumerator = function(ChIPpos, ChIPsig) {
  products = ChIPpos * ChIPsig
  rppNumerator = sum(products)
  return(rppNumerator)  
}

getOperonRPP = function(ChIPdf, operonsDF, operonID) {
  operon = operonsDF %>% filter(OperonID == operonID)
  opStart = operon$start[1]
  opEnd = operon$end[1]
  strand = operon$strand[1]
  ChIPData = ChIPdf %>% filter(position>=opStart, position<=opEnd)
  
  ChIPsignal = ChIPData$log2Enrich
  if (strand == "-") {
    ChIPsignal = rev(ChIPsignal)
  }
  
  ChIPposition = ChIPData$position
  ChIPposition = ChIPposition - min(ChIPposition) + 1
  
  ChIPSum = sum(ChIPsignal)
  denom = ChIPSum*max(ChIPposition)
  
  numer = getRPPNumerator(ChIPposition, ChIPsignal)
  
  RPP = numer/denom
  return(RPP)
}

getAllRPPs = function(df) {
  RPPvec = NA
  operonIDvec = NA
  for (i in 1:nrow(operonsDF)) {
    operonIDvec[i] = operonsDF$OperonID[i]
    RPPvec[i] = getOperonRPP(df, operonsDF, operonIDvec[i])
  }
  return(tibble(OperonID=operonIDvec, RPP=RPPvec))
}

getLociRPPs = function(locus_tags, enrichmentDF, preInversionFeatures) {
  rpps = NA
  for (i in 1:length(locus_tags)) {
    rpps[i] = getGeneRPP(.locus_tag=locus_tags[i], ChIPdf=enrichmentDF, featureDF=preInversionFeatures)
  }
  return(rpps)
}

getGeneRPP = function(.locus_tag, ChIPdf, featureDF) {
  feature = featureDF %>% filter(locus_tag == .locus_tag)
  featStart = feature$start[1]
  featEnd = feature$end[1]
  strand = feature$strand[1]
  ChIPData = ChIPdf %>% filter(position>=featStart, position<=featEnd)
  
  ChIPsignal = ChIPData$log2Enrich
  if (strand == "-") {
    ChIPsignal = rev(ChIPsignal)
  }
  
  ChIPposition = ChIPData$position
  ChIPposition = ChIPposition - min(ChIPposition) + 1
  
  ChIPSum = sum(ChIPsignal)
  denom = ChIPSum*max(ChIPposition)
  
  numer = getRPPNumerator(ChIPposition, ChIPsignal)
  
  RPP = numer/denom
  return(RPP)
}

getLocusCV = function(ChIPdf, startPosition, endPosition) {
  signal = ChIPdf$log2Enrich[ChIPdf$position>=startPosition & ChIPdf$position<=endPosition]
  return(sd(signal)/mean(signal))
}

getAllCVs = function(df, genesDF) {
  CVvec = NA
  locusIDvec = NA
  for (i in 1:nrow(genesDF)) {
    locTag = genesDF$locus_tag[i]
    locusIDvec[i] = locTag
    startPos = genesDF$start[genesDF$locus_tag==locTag]
    endPos = genesDF$end[genesDF$locus_tag==locTag]
    CVvec[i] = getLocusCV(df, startPos, endPos)
  }
  return(tibble(locus_tag=locusIDvec, CV=CVvec))
}

# dyGraphGeneAnnotation <- function(dygraph, midpoint, text) {
#   dygraph = dygraph %>%
#     dyAnnotation(midpoint, text, attachAtBottom = TRUE, width = 60)
#   return(dygraph)
# }

