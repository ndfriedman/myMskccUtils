#the first R script Im making for work actually... after 7 months i caved
library(ggplot2)
library(grid)
require(cowplot)
library(egg)

#utils Copy

#given a df generates a ggplotbarplot as specified here
generate_ggplot_bar <- function(df, title, hideLegend = TRUE) {
  
  barColorPalette = c("#D3D3D3", #other signatures light gray 
                      "#FF0000", #apobec is red
                      "#FFA500", #smoking is orange
                      "#FFF600", #UV is
                      "#ADFF2F", #Pole is
                      "#267574", #mmr is blue-green
                      "#2A52BE", #aging is blue 
                      "#FF1493" #BRCA sig deep pink
  )
  
  yLabel <- 'Signature fractions'
  plt <- ggplot(df, aes(x = reorder(dmp_sample, -brca_signature), #order x axis by the brca signature
                        y = mean, #the dataframe has a different mean for each signature we give a shit about
                        fill=factor(signatureName, levels=c("Other", "signature_APOBEC", "Smoking", 'signature_UV', 'signature_POLE', "signature_MMR/MSI", "Age",  "Brca")))) + #order bar charts how I want 
    geom_bar(stat = "identity")+
    scale_fill_manual(values=barColorPalette) +
    theme( #rotate axes, change size
      axis.title.x=element_blank(), #make there be no x title
      axis.title.y=element_text(size=4),
      legend.title=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
      axis.text.x=element_blank(), axis.text.y=element_blank()
    ) + 
    ggtitle(title)+ 
    guides(shape = guide_legend(override.aes = list(size = .1))) +
    ylab("Signature Means") +
    theme(plot.title = element_text(size = 5))
  
  if(hideLegend == TRUE){plt <- plt + theme(legend.position="none")} #kill the legend conditionally based on function params
  else{
    plt <- plt+ theme(legend.text=element_text(size=4),
                      legend.title = element_text(size=5, ))+
      guides(color = guide_legend(override.aes = list(size=.1)))+
      guides(fill=guide_legend(title="Signature"))
  }
  return(plt)
}

get_empty_theme <- function(){
  return(
    theme( #rotate axes, change size
      axis.title.x=element_blank(), #make there be no x title
      axis.title.y=element_text(size=4),
      legend.title=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
      axis.text.x=element_blank(), axis.text.y=element_blank()
    )
  )
}

make_bar_comparison <- function(df, title,  hideLegend = TRUE, textSize=1) {
  
  #ordering
  #c('brca', 'Age', 'signature_MMR/MSI', 'signature_UV', 'signature_POLE', 'Smoking', 'signature_5', 'signature_8', 'other', 'signature_APOBEC')
  barColorPalette = c(
    "#FF1493", #brca
    "#00FFFF", #age  
    "#267574", #mmr
    "#FFF600", #Uv
    "#ADFF2F", #POLE
    "#FFA500", #smoking 
    "#2A52BE", #sig5
    "#551A8B", #sig8
    "#D3D3D3", #OTHER
    "#FF0000" #APOBEC
  )
  
  plt <- ggplot(df)+
    geom_bar(aes(x = reorder(dmp_sample, -brca_signature), y=-otherSigMagnitude, fill=factor(otherPredominantSignature,
                                                                                             levels=c('brca', 'Age', 'signature_MMR/MSI', 'signature_UV', 'signature_POLE', 'Smoking', 'signature_5', 'signature_8', 'other', 'signature_APOBEC'))),
             stat = "identity")+
    geom_bar(aes(x = reorder(dmp_sample, -brca_signature), y=brca_signature), stat="identity",fill="#FF1493")+
    geom_hline(yintercept=0, color = "black", size=.25)+
    
    scale_fill_manual(values=barColorPalette, drop = FALSE)+ 
    ggtitle(title)+
    theme(plot.title = element_text(size = 5)) +
    get_empty_theme()
  #theme(axis.text.x = element_text(angle = 60, hjust = 1, size=textSize))
  if(hideLegend == TRUE){
    plt <- plt + theme(legend.position="none")} #kill the legend conditionally based on function params
  else{
    plt <- plt + guides(fill=guide_legend(title="Signature Type"))
  }
  plt <- plt + labs(y = "Signature Magnitude")+ 
    theme(axis.title.y = element_text(angle = 0))
  return(plt)
}

make_purity_bar <- function(df, hideLegend = TRUE){
  plt <- ggplot(df, aes(x = reorder(dmp_sample, -brca_signature), y=0))+
    geom_tile(aes(fill=purity,
                  width=0.7, height=0.7), size=5)+
    scale_fill_gradient(low = "white", high = "black")+
    get_empty_theme()
  if(hideLegend == TRUE){
    plt <- plt + theme(legend.position="none")}
  else{
    plt <- plt + guides(fill=guide_legend(title="Purity"))
  }
  plt <- plt + labs(y = "Purity")+ 
    theme(axis.title.y = element_text(angle = 0))
  return(plt)
}

#noob R way to do this
generate_ordered_color_mapping <- function(df){
  colorMapping <- vector(mode="list", length=10)
  names(colorMapping) <- c("Bladder Cancer", "Breast Cancer", "Colorectal Cancer", 'Endometrial Cancer', 'Melanoma', "Non-Small Cell Lung Cancer", "Other",  "Ovarian Cancer", "Pancreatic Cancer", "Prostate Cancer")
  colors = c("#FAE5D3", "#3498DB", "#E8DAEF", "#D0ECE7", "#FCDC3B", "#D98880", "#D5DBDB", "#1ABC9C", "#2ECC71", "#7D3C98")
  colorMapping[[1]] <- "#FAE5D3" #bladder is light red 
  colorMapping[[2]] <- "#3498DB" #breast is dark purple
  colorMapping[[3]] <- "#E8DAEF" #colorectal is raw siena
  colorMapping[[4]] <- "#D0ECE7" #endometrial is mars orange
  colorMapping[[5]] <- "#FCDC3B" #melanoma is pineapple
  colorMapping[[6]] <- "#D98880" #lung is light red
  colorMapping[[7]] <- "#D5DBDB" #other is bright gray
  colorMapping[[8]] <- "#1ABC9C" #ovarian is dark green
  colorMapping[[9]] <- "#2ECC71" #pancreas is dark blue
  colorMapping[[10]] <- "#7D3C98" #prostate is dark red)
  order<- unique(c(as.character(df$cancerTypeAdjusted)))
  order<-sort(order)
  x <- c()
  for (i in order) {
    x <- append(x,colorMapping[[i]])
  }
  return(x)
}

#TODO fix weirdness for legends in the new way I am doing things
generate_ggplot_tiles <- function(df, fillArg, hideLegend = TRUE, mode = "cancerTypeTiles", textSizeParam=1) {
  tileColorPalette <- NULL 
  if(mode == "cancerTypeTiles"){
    tileColorPalette = generate_ordered_color_mapping(df)  
  }
  if(mode == "allelicStatus"){
    tileColorPalette = c("#0000FF", "#CCCCFF", "#FF9933", "#FFCCCC", "#E6E6E6")
  }
  if(mode == "brca"){
    tileColorPalette = c("#D3D3D3", "#98FB98", "#008000")
  }
  textSize <- textSizeParam
  plt <- ggplot(df, aes(x = reorder(dmp_sample, -brca_signature), y=0)) + 
    geom_tile(aes_string(fill=fillArg,
                         width=0.7, height=0.7), size=5)+
    scale_fill_manual(values=tileColorPalette) +
    theme(
      axis.text.y=element_blank(), #try to remove everything to make it a pure geometric shape
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.line = element_blank(),
      panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
      axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
  theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    theme_void()
  
  #the mode affects our text:
  if(mode == "brca"){plt <- plt + theme(axis.text.x=element_blank())}
  else{plt <- plt + theme(axis.text.x = element_text(angle = 60, hjust = 1, size=1))}
  if(hideLegend == TRUE){plt <- plt + theme(legend.position="none")} #kill the legend conditionally based on function params
  else{
    plt <- plt+ theme(legend.text=element_text(size=4),
                      legend.title = element_text(size=5))+
      guides(color = guide_legend(override.aes = list(size=.1)))+
      guides(fill=guide_legend(title=fillArg))
  }
  plt <- plt + labs(y = fillArg)+ 
    theme(axis.title.y = element_text(angle = 0, size=5))
  return(plt)
}

generate_mut_burden_bar <- function(df){
  plt <- ggplot(df, aes(x = reorder(dmp_sample, -brca_signature),
                        y=nMutationsAdj, fill=nMutClipped)
  )+
    geom_bar(stat = "identity")+
    scale_fill_manual(values=c("#2F4F4F","#000000"))+
    theme(
      axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
      axis.text.x=element_blank(), axis.text.y=element_blank(),
      axis.title.x=element_blank(),axis.title.y=element_blank(),
      axis.line = element_blank(),
      panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank())
  theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    theme_void()
  plt <- plt + theme(legend.position="none")
  plt <- plt + labs(y = "mutation burden")+ 
    theme(axis.title.y = element_text(angle = 0, size=5))
  return(plt)
}

#omnibus plotting function
make_full_plot <- function(df1,df2,df3,df4,df5, #working with the 5 dataframe version
                           df1Title, df2Title, df3Title, df4Title, df5Title,
                           tilePanelMode, saveFilename
){
  
  fillArg <- ""
  if(tilePanelMode == "cancerTypeTiles"){
    fillArg <- "cancerTypeAdjusted"
  }
  if(tilePanelMode == "allelicStatus"){
    fillArg <- "biallelic_class"
  }
  
  topMargin <- unit(c(0,-.5,-.5,-.5), "lines")
  myStandardMargin <- unit(c(.3,-.5,-.5,-.5), "lines")
  
  tilePanelBar1 <- generate_ggplot_tiles(df1, fillArg, TRUE, tilePanelMode)
  tilePanelBar2 <- generate_ggplot_tiles(df2, fillArg, TRUE, tilePanelMode)
  tilePanelBar3 <- generate_ggplot_tiles(df3, fillArg, TRUE, tilePanelMode)
  tilePanelBar4 <- generate_ggplot_tiles(df4, fillArg, TRUE, tilePanelMode)
  tilePanelBar5 <- generate_ggplot_tiles(df5, fillArg, TRUE, tilePanelMode)
  
  finalPlotSansLegend <- ggarrange(
    
    make_bar_comparison(df1, df1Title) + theme(plot.margin = topMargin),
    generate_mut_burden_bar(df1) + theme(plot.margin = myStandardMargin),
    make_purity_bar(df1) + theme(plot.margin = myStandardMargin),
    tilePanelBar1,
    
    make_bar_comparison(df2, df2Title) + theme(plot.margin =topMargin),
    generate_mut_burden_bar(df2) + theme(plot.margin = myStandardMargin),
    make_purity_bar(df2) + theme(plot.margin = myStandardMargin),
    tilePanelBar2,
    
    make_bar_comparison(df3, df3Title) + theme(plot.margin =topMargin),
    generate_mut_burden_bar(df3) + theme(plot.margin = myStandardMargin),
    make_purity_bar(df3) + theme(plot.margin = myStandardMargin),
    tilePanelBar3,
    
    make_bar_comparison(df4, df4Title) + theme(plot.margin =topMargin),
    generate_mut_burden_bar(df4) + theme(plot.margin = myStandardMargin),
    make_purity_bar(df4) + theme(plot.margin = myStandardMargin),
    tilePanelBar4,
    
    make_bar_comparison(df5, df5Title) + theme(plot.margin =topMargin),
    generate_mut_burden_bar(df5) + theme(plot.margin = myStandardMargin),
    make_purity_bar(df5) + theme(plot.margin = myStandardMargin),
    tilePanelBar5,
    
    heights = c(5,.5,.5,.5,
                5,.5,.5,.5,
                5,.5,.5,.5,
                5,.5,.5,.5,
                5,.5,.5,.5)
  )
  
  legendSigs <- get_legend(make_bar_comparison(df1, df1Title, FALSE))
  legendPurity <- get_legend(make_purity_bar(df1, FALSE))
  legendClassOrCancerType <- get_legend(generate_ggplot_tiles(df3, fillArg, FALSE, tilePanelMode)) 
  
  finalPlot <- plot_grid(finalPlotSansLegend, legendSigs, legendPurity, legendClassOrCancerType,
                         align='hv', ncol=4,
                         rel_widths = c(1,.2, .2, .2), scale = 0.9)
  ggsave(saveFilename, plot=finalPlot)
}


