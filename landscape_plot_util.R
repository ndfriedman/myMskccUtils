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
    get_adjusted_theme()+
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

get_adjusted_theme <- function(){
  return(
    theme( #rotate axes, change size
      #panel.grid.major.x = element_line(colour = "grey"),
      #panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      #axis.line = element_line(colour = 'black', size = 1),
      axis.title.x=element_blank(), #make there be no x title
      axis.title.y=element_text(size=2),
      legend.title=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y = element_text(size=2,angle=45, hjust=1))
    )
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
                                xAxisValParam,
                                yAxisValParam1, yAxisValParam2, #the parameters that define the upper and lower bar respectively
                                yAxisFillParam, #the parameter we use to choose the color for the lower bar
                                orderingValParam, #the parameter for ordering the values on the X axis
                                mainColor="#FF0000", #the color for the bar of the main signature that goes upwards
                                coloringSpecCols, barColorPalette, #used for coloring the bar chart
                                title,  hideLegend = TRUE, textSize=1, 
                                pointPlotValParam = NA #parameter for what value to use for ggpoints if we want them
                                ) {
  
  pointPlotSize = 20/nrow(df) #parameter for the size of the points if we include them
  #first make the ordering column before any of the rest of this bullshit
  #CHANGE THE DF Columns
  df <- my.rename(df, xAxisValParam, "xAxisVal")
  df <- my.rename(df, yAxisValParam1, "yAxisVal1")
  df <- my.rename(df, yAxisValParam2, "yAxisVal2")
  df <- my.rename(df, yAxisFillParam, "yAxisFill")
  if(!is.na(pointPlotValParam)){ #rename the parameter for the point plot in the graph
    df <- my.rename(df, pointPlotValParam, "pointPlotVal")
  }
  plt <- ggplot(df)+
    
    #The first bar of the predominant signature in the positive direction
    geom_bar(aes(x = reorder(xAxisVal, -orderingVal), y=yAxisVal1), stat="identity",fill=mainColor)+
    #The second bar of the second predominant signature in the negative direction
    geom_bar(aes(x = reorder(xAxisVal, -orderingVal), y=-yAxisVal2, 
    fill = factor(yAxisFill, levels=coloringSpecCols) #color the other signature column by which signature it is                                                                                     #levels=c('brca', 'Age', 'signature_MMR/MSI', 'signature_UV', 'signature_POLE', 'Smoking', 'signature_5', 'signature_8', 'other', 'signature_APOBEC')
    ),
    stat = "identity")+
    ylim(-1,1)+  
    
    geom_hline(yintercept=0, color = "black", size=.25)+
    scale_fill_manual(values=barColorPalette, drop = FALSE)+ 
    ggtitle(paste(title, "// n cases: ", nrow(df)))+
    get_adjusted_theme()+
    theme(plot.title = element_text(size = 5)) +
    theme(axis.line.x = element_line(color="white", size = 2))
  if(hideLegend == TRUE){
    plt <- plt + theme(legend.position="none")} #kill the legend conditionally based on function params
  else{
    plt <- plt + guides(fill=guide_legend(title="Signature Type"))
  }
  plt <- plt + labs(y = "Signature Magnitude")+ 
    theme(axis.title.y = element_text(angle = 0))
  if(!is.na(pointPlotValParam)){
    plt <- plt + geom_point(shape=4, size=pointPlotSize, colour = "#967BB6", aes(x = reorder(xAxisVal, -orderingVal), y=pointPlotVal))
  }
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
#TODO fix the fill arg not found error
generate_ggplot_tiles <- function(df, xAxisValParam, xAxisOrderingParam, fillArgParam, hideLegend = TRUE,
                                  mode = "cancerTypeTiles", textSizeParam=1,tileColorPalette=NA, orderSpecCols=NA) {
  
  df <- my.rename(df, xAxisValParam, "xAxisVal")
  df <- my.rename(df, xAxisOrderingParam, "xAxisOrdering")
  df <- my.rename(df, fillArgParam, "fillArg")
  
  textSize <- textSizeParam
  
  plt <- ggplot(df, aes(x = reorder(xAxisVal, -xAxisOrdering), y=0)) 
  if(!is.na(orderSpecCols)){
    plt <- plt + geom_tile(aes(
                         fill= factor(fillArg, levels=orderSpecCols),
                         width=0.7, height=0.7), size=5)
  }
  else{
    plt <- plt + geom_tile(aes(
      fill= fillArg,
      width=0.7, height=0.7), size=5)
  }
  plt <- plt + theme(
      axis.text.y=element_blank(), #try to remove everything to make it a pure geometric shape
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.line = element_blank(),
      panel.background=element_blank(),panel.border=element_blank(),
      #panel.grid.major=element_blank(),
      #panel.grid.major.x = element_line(colour = "grey"),
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
      guides(fill=guide_legend(title=fillArgParam))
  }
  plt <- plt + labs(y = fillArgParam)+ 
    theme(axis.title.y = element_text(angle = 0, size=5))
  if(!is.na(tileColorPalette)){
    plt <- plt + scale_fill_brewer(palette=tileColorPalette)
  }
  return(plt)
}

generate_gradient_tiles <- function(df, xAxisValParam, xAxisOrderingParam, fillArgParam, hideLegend = TRUE) {
    #a function designed to genrate a row of red to blue tiles
    df <- my.rename(df, xAxisValParam, "xAxisVal")
    df <- my.rename(df, xAxisOrderingParam, "xAxisOrdering")
    df <- my.rename(df, fillArgParam, "fillArg")
    
    plt <- ggplot(df, aes(x = reorder(xAxisVal, -xAxisOrdering), y=0)) + 
      geom_tile(aes(fill= fillArg,
                    width=0.7, height=0.7), size=5)+
      scale_fill_gradient(low = "blue", high = "red")+
      get_empty_theme()
    plt <- plt + labs(y = "MSI Score")+ 
      theme(axis.title.y = element_text(angle = 0, size=5))
    if(hideLegend == TRUE){plt <- plt + theme(legend.position="none")}
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
  #df <- my.rename(df, yAxisFillParam, "yAxisFill")
  
  plt <- ggplot(df, aes(x = reorder(xAxisVal, -xAxisOrdering),
                        y=yAxisVal, fill="#000000")
  )+
    geom_bar(stat = "identity")+
    scale_y_log10()+
    scale_fill_manual(values=c("#2F4F4F","#000000"))+
    get_adjusted_theme()+
    theme(
      axis.ticks.x=element_blank(), 
      axis.text.x=element_blank(),
      axis.title.x=element_blank(),
      panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank())
  theme_bw() +
    theme(axis.line = element_line(colour = "black"))
  plt <- plt + theme(legend.position="none")
  plt <- plt + labs(y = "log mutation burden")+ 
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


plot_panel <- function(df,
                         #parameters for plotting the signature comparisson bars
                         xAxisValParam_=NA,
                         barParam1_=NA,
                         barParam2_=NA,
                         yAxisFillParam_=NA,
                         orderingValParam_=NA,
                         coloringSpecCols_=NA, barColorPalette_=NA,
                         title_=NA,
                         primaryBarColor_=NA,
                         pointPlotValParam_ = NA,
                         mutBurdenBarParam_=NA,
                         gradientBarParam_=NA, #parameter for gradient bar
                         tileFillCancerParam_=NA,
                         cancerTypeColorPalette_=NA, 
                         tileFillAllelicStatusParam_=NA,
                         allelicTypeColorPalette_=NA,
                         tileMode_=NA,
                         legendMode=FALSE
                         ){
  hideLegendParam=TRUE
  if(legendMode == TRUE){hideLegendParam=FALSE}
  panelReturnList <- c() #These lists are returned by the function to build the panel
  returnListCntr <- 1
  hs <- list()
  if(!is.na(barParam1_)){
    compBar <- make_bar_comparison(
      df, 
      xAxisValParam=xAxisValParam_,
      yAxisValParam1=barParam1_,
      yAxisValParam2=barParam2_,
      yAxisFillParam=yAxisFillParam_,
      orderingValParam=orderingValParam_,
      coloringSpecCols=coloringSpecCols_, barColorPalette=barColorPalette_,
      title=title_,
      mainColor=primaryBarColor_,
      hideLegend = hideLegendParam,
      pointPlotValParam = pointPlotValParam_) #todo make the point size be a function of df size
    if(legendMode == TRUE){panelReturnList[[returnListCntr]] <- get_legend(compBar)}
    else{panelReturnList[[returnListCntr]] <- compBar}
    returnListCntr <- returnListCntr + 1
  }
  if(!is.na(mutBurdenBarParam_)){
    mutBurdenBar <- generate_mut_burden_bar(
      df, 
      xAxisValParam=xAxisValParam_, 
      xAxisOrderingParam=orderingValParam_, 
      yAxisValParam=mutBurdenBarParam_,
      yAxisFillParam="#000000") #todo change if theres ever a reason we will be expressing the mut burden bar in different colors
    panelReturnList[[returnListCntr]] <- mutBurdenBar
    returnListCntr <- returnListCntr + 1  
  }
  if(!is.na(gradientBarParam_)){
    gradientBar <- generate_gradient_tiles(df, 
                            xAxisValParam=xAxisValParam_, 
                            xAxisOrderingParam=orderingValParam_, 
                            fillArgParam=gradientBarParam_, 
                            hideLegend = hideLegendParam) 
    
    if(legendMode == TRUE){panelReturnList[[returnListCntr]] <- get_legend(gradientBar)}
    else{panelReturnList[[returnListCntr]] <- gradientBar}
    returnListCntr <- returnListCntr + 1
  }
  if(!is.na(tileFillAllelicStatusParam_)){                    
    tileBarAllelic <- generate_ggplot_tiles(
      df, 
      xAxisValParam=xAxisValParam_, 
      xAxisOrderingParam=orderingValParam_,
      fillArgParam=tileFillAllelicStatusParam_, 
      hideLegend = hideLegendParam,
      mode = "brca", 
      textSizeParam=.5,
      tileColorPalette=allelicTypeColorPalette_)
    if(legendMode == TRUE){panelReturnList[[returnListCntr]] <- get_legend(tileBarAllelic)}
    else{panelReturnList[[returnListCntr]] <- tileBarAllelic}
    returnListCntr <- returnListCntr + 1
  }
  if(!is.na(tileFillCancerParam_)){                    
    tileBar <- generate_ggplot_tiles(
                          df, 
                          xAxisValParam=xAxisValParam_, 
                          xAxisOrderingParam=orderingValParam_,
                          fillArgParam=tileFillCancerParam_, 
                          hideLegend = FALSE,
                          mode = "cancerTypeTiles", 
                          textSizeParam=.5,
                          tileColorPalette=cancerTypeColorPalette_)
    if(legendMode == TRUE){panelReturnList[[returnListCntr]] <- get_legend(tileBar)}
    else{panelReturnList[[returnListCntr]] <- tileBar}
    returnListCntr <- returnListCntr + 1
  }
  return(panelReturnList)  
}


