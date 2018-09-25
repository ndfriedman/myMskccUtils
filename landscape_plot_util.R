#the first R script Im making for work actually... after 7 months i caved
library(ggplot2)
library(grid)
require(cowplot)
library(egg)
require(plyr)

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

#makes a plot comparing the signature of interest of the predominant signature
make_bar_comparison <- function(df, 
                                xAxisValParam, yAxisValParam1, yAxisValParam2, yAxisFillParam, orderingValParam,
                                mainColor="#FF0000", coloringSpecCols, barColorPalette, #used for coloring the bar chart
                                title,  hideLegend = TRUE, textSize=1) {
  #first make the ordering column before any of the rest of this bullshit
  df$orderingVal <-df[[orderingValParam]]
  #CHANGE THE DF Columns
  df <- my.rename(df, xAxisValParam, "xAxisVal")
  df <- my.rename(df, yAxisValParam1, "yAxisVal1")
  df <- my.rename(df, yAxisValParam2, "yAxisVal2")
  df <- my.rename(df, yAxisFillParam, "yAxisFill")
  
  plt <- ggplot(df)+
    
    #The first bar of the predominant signature in the positive direction
    geom_bar(aes(x = reorder(xAxisVal, -orderingVal), y=yAxisVal1), stat="identity",fill=mainColor)+
    
    #The second bar of the second predominant signature in the negative direction
    geom_bar(aes(x = reorder(xAxisVal, -orderingVal), y=-yAxisVal2, 
    fill = factor(yAxisFill, levels=coloringSpecCols) #color the other signature column by which signature it is                                                                                     #levels=c('brca', 'Age', 'signature_MMR/MSI', 'signature_UV', 'signature_POLE', 'Smoking', 'signature_5', 'signature_8', 'other', 'signature_APOBEC')
    ),
    stat = "identity")+
    
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

make_percentage_bar <- function(#makes a bar chart for a 0-100% percentage (ie purity)
                                df, 
                                xAxisValParam, xAxisOrderingParam, fillValParam,
                                titleVal, colorLow="white", colorHigh="#8B0000",
                                hideLegend = TRUE){ #makes a bar chart for a 0-100% percentage (ie purity)
  
  df <- my.rename(df, xAxisValParam, "xAxisVal")
  df <- my.rename(df, xAxisOrderingParam, "xAxisOrdering")
  df <- my.rename(df, fillValParam, "fillVal")
  
  plt <- ggplot(df, aes(x = reorder(xAxisVal, -xAxisOrdering), y=0))+
    geom_tile(aes(fill=fillVal,
                  width=0.7, height=0.7), size=5)+
    scale_fill_gradient(low = colorLow, high = colorHigh)+
    get_empty_theme()
  if(hideLegend == TRUE){
    plt <- plt + theme(legend.position="none")}
  else{
    plt <- plt + guides(fill=guide_legend(title=titleVal))
  }
  plt <- plt + labs(y = titleVal)+ 
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

my.rename <- function(df, old.name, new.name){ #R SUCKS!!!!!! heres my renaming function cause god forbid there would be an easy or intiuitive way to do this with R
  names(df)[names(df) == old.name] <- new.name 
  return(df)
}

generate_mut_burden_bar <- function(df, xAxisValParam, xAxisOrderingParam, yAxisValParam, yAxisFillParam){
  
  #Rename columns as needed
  df <- my.rename(df, xAxisValParam, "xAxisVal")
  df <- my.rename(df, xAxisOrderingParam, "xAxisOrdering")
  df <- my.rename(df, yAxisValParam, "yAxisVal")
  df <- my.rename(df, yAxisFillParam, "yAxisFill")
  
  plt <- ggplot(df, aes(x = reorder(xAxisVal, -xAxisOrdering),
                        y=yAxisVal, fill=yAxisFill)
  )+
    geom_bar(stat = "identity")+
    #scale_y_log10()+
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

plot_proportion_bar <- function(df, xAxisValParam, xAxisOrderingParam, yAxisValParam, label){
  #Rename columns as needed
  df <- my.rename(df, xAxisValParam, "xAxisVal")
  df <- my.rename(df, xAxisOrderingParam, "xAxisOrdering")
  df <- my.rename(df, yAxisValParam, "yAxisVal")
  
  plt <- ggplot(df, aes(x = reorder(xAxisVal, -xAxisOrdering),
                        y=yAxisVal)
  )+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=c("#000000"))+
  geom_hline(aes(yintercept=1), size=0.1)+
  geom_hline(aes(yintercept=.5), size=0.1, colour="#D3D3D3")+
  get_empty_theme()
  plt <- plt + labs(y = label)+ 
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


