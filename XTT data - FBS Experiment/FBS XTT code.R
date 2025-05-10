# Setup ####

library(tidyverse)
library(patchwork)

overwrite <- T
plotWindow <- F


# Data Wrangling ####

absorb <- read_csv('all.csv') |> 
  select(!BLK) |> 
  pivot_longer(!Hours,names_to = 'Serum condition',values_to = 'Abs') |> 
  mutate('Serum condition'=paste0(sprintf("%.1f", round(as.numeric(`Serum condition`)*100,1)),'%'))|>
  mutate('Serum condition'=str_replace_all(`Serum condition`,'^0.0%','0.0% + HSA')) |> 
  mutate('Serum condition'=factor(`Serum condition`,levels=c('0.0% + HSA', '0.1%', '0.5%', '1.0%', '2.0%', '10.0%'))) |> 
  filter(Abs!=min(Abs)) # manually removing 1 outlier

conds <-unique(absorb$`Serum condition`)
hrs <- absorb |> 
  filter(Hours!=2) |> 
  select(Hours) |> 
  pull() |> 
  unique()


# Absorbance plots ####

abs2 <- list()
for(c in 1:length(conds)){
  abs2[[c]] <- absorb |> 
    filter(Hours==2) |> 
    mutate('Serum condition'=conds[c])
}
abs2 <- bind_rows(abs2)

absorbPlotting <- absorb |> 
  filter(Hours!=2) |> 
  bind_rows(abs2) |> 
  rowwise() |> 
  mutate(plot=any(Hours!=2,`Serum condition`=='10.0%'))

## All absorbances ####

absPlot <- absorbPlotting |> 
  ggplot(aes(x=Hours,y=`Abs`,color=`Serum condition`))+
  # facet_wrap(~`Serum condition`)+
  geom_point(aes(shape=plot))+
  scale_x_continuous(breaks=c(0,hrs),limits=c(0,max(hrs)+2))+
  scale_y_continuous(limits=c(0,max(absorbPlotting$Abs)))+
  geom_smooth(method='lm',formula=y~x,se=F)+
  scale_shape_manual(values=c(NA,16),guide='none')+
  labs(y='Corrected absorbance',title='All corrected absorbances')+
  theme_bw()
if(overwrite) {
  absPlot <- absPlot+scale_color_discrete(guide='none')
}

if(plotWindow) absPlot


## Mean absorbances ####

meanFunk <- function(x,i){mean(x[i])} # for boot
se <- function(data,times=2000){
  b <- boot::boot(data,meanFunk,times)
  return(sd(b$t))
}

seFC <- function(dataA,dataA0,times=100000){
  est <- rep(NA,times)
  for(b in 1:times){
    A <- sample(dataA,replace=T)
    A0 <- sample(dataA0,replace=T)
    est[b] <- mean(A)/mean(A0)
  }
  return(sd(est))
}

means <- absorbPlotting |> 
  group_by(Hours,`Serum condition`) |> 
  summarise('Mean Abs'=mean(Abs),'Mean Abs SE'=se(Abs)) |> 
  rowwise() |> 
  mutate(plot=any(Hours!=2,`Serum condition`=='10.0%'))

absMeansPlot <- means |>
  ggplot(aes(x=Hours,y=`Mean Abs`,color=`Serum condition`))+
  geom_point(size=2)+
  geom_path(linewidth=1.1)+
  geom_errorbar(aes(ymin=`Mean Abs`-`Mean Abs SE`, ymax=`Mean Abs`+`Mean Abs SE`),width=max(hrs)/10)+
  labs(y='Corrected absorbance',title='Mean corrected absorbances')+
  scale_x_continuous(breaks=c(0,hrs))+
  scale_y_continuous(limits=c(0,max(absorbPlotting$Abs)))+
  theme_bw()

if(plotWindow) absMeansPlot

# Fold changes plots ####

init <- absorb |> 
  filter(Hours==2) |> 
  select(Abs) |> 
  pull()

folds <- absorb |> 
  filter(Hours!=2) |> 
  mutate('Fold change'=Abs/mean(init))

folds2 <- list()
for(c in 1:length(conds)){
  folds2[[c]] <- absorb |> 
    filter(Hours==2) |> 
    mutate('Fold change'=Abs/mean(init)) |> 
    mutate('Serum condition'=conds[c])
}
folds2 <- bind_rows(folds2)

foldsPlotting <- folds |> 
  bind_rows(folds2) |> 
  rowwise() |> 
  mutate(plot=any(Hours!=2,`Serum condition`=='10.0%'))


## All log fold changes ####

logFoldsPlotting <- folds |> 
  bind_rows(folds2) |> 
  mutate('Log fold change'=log2(`Fold change`)) |> 
  select(!`Fold change`) |> 
  rowwise() |> 
  mutate(plot=any(Hours!=2,`Serum condition`=='10.0%'))

logFoldsPlot <- logFoldsPlotting |> 
  ggplot(aes(x=Hours,y=`Log fold change`,color=`Serum condition`))+
  # facet_wrap(~`Serum condition`)+
  geom_point(aes(shape=plot))+
  scale_x_continuous(breaks=c(0,hrs),limits=c(-2,max(hrs)+4))+
  scale_y_continuous(limits=c(min(logFoldsPlotting$`Log fold change`),
                              max(logFoldsPlotting$`Log fold change`)+0.2))+
  geom_smooth(method='lm',formula=y~x,se=F)+
  scale_shape_manual(values=c(NA,16),guide='none')+
  labs(y=bquote(log[2]-'fold change'),title=bquote('All'~log[2]-'fold changes'))+
  theme_bw()
if(overwrite){
  logFoldsPlot <- logFoldsPlot+scale_color_discrete(guide='none')
}

if(plotWindow) logFoldsPlot


## Mean log fold changes ####

meanLogFolds <- logFoldsPlotting |> 
  group_by(Hours,`Serum condition`) |> 
  summarise('Mean log fold change'=mean(`Log fold change`),
            'CI min'=`Mean log fold change` - se(`Log fold change`),
            'CI max'=`Mean log fold change` + se(`Log fold change`)) |> 
  rowwise() |> 
  mutate(plot=any(Hours!=2,`Serum condition`=='10.0%'))

logFoldMeansPlot <- meanLogFolds |> 
  ggplot(aes(x=Hours,y=`Mean log fold change`,color=`Serum condition`))+
  geom_point(size=2)+
  geom_path(linewidth=1.1)+
  geom_errorbar(aes(ymin=`CI min`, ymax=`CI max`),width=max(hrs)/10)+
  scale_x_continuous(breaks=c(0,hrs),limits=c(-2,max(hrs)+4))+
  scale_y_continuous(limits=c(min(logFoldsPlotting$`Log fold change`),
                              max(logFoldsPlotting$`Log fold change`)+0.2))+
  labs(y=bquote(log[2]-'fold change'),title=bquote('Mean'~log[2]-'fold changes'))+
  theme_bw()

if(plotWindow) logFoldMeansPlot


## All fold changes ####

foldsPlot <- foldsPlotting |> 
  ggplot(aes(x=Hours,y=`Fold change`,color=`Serum condition`))+
  # facet_wrap(~`Serum condition`)+
  geom_point(aes(shape=plot))+
  scale_x_continuous(breaks=c(0,hrs),limits=c(-2,max(hrs)+4))+
  scale_y_continuous(limits=c(min(foldsPlotting$`Fold change`),
                              max(foldsPlotting$`Fold change`)+0.2))+
  geom_smooth(method='lm',formula=y~x,se=F)+
  scale_shape_manual(values=c(NA,16),guide='none')+
  labs(y='Fold change',title='All fold changes')+
  theme_bw()
if(overwrite){
  foldsPlot <- foldsPlot+scale_color_discrete(guide='none')
}

if(plotWindow) foldsPlot


## Mean fold changes ####

meanFolds <- foldsPlotting |> 
  group_by(Hours,`Serum condition`) |> 
  summarise('Mean fold change'=mean(`Fold change`),
            'CI min'=`Mean fold change` - seFC(`Fold change`,init),
            'CI max'=`Mean fold change` + seFC(`Fold change`,init)) |> 
  rowwise() |> 
  mutate(plot=any(Hours!=2,`Serum condition`=='10.0%'))

foldMeansPlot <- meanFolds |> 
  ggplot(aes(x=Hours,y=`Mean fold change`,color=`Serum condition`))+
  geom_point(size=2)+
  geom_path(linewidth=1.1)+
  geom_errorbar(aes(ymin=`CI min`, ymax=`CI max`),width=max(hrs)/10)+
  scale_x_continuous(breaks=c(0,hrs),limits=c(-2,max(hrs)+4))+
  scale_y_continuous(limits=c(min(foldsPlotting$`Fold change`),
                              max(foldsPlotting$`Fold change`)+0.2))+
  labs(y='Fold change',title='Mean fold changes')+
  theme_bw()

if(plotWindow) foldMeansPlot

foldMeansPlot0.1 <- meanFolds |> 
  filter(`Serum condition`=='0.1%') |> 
  ggplot(aes(x=Hours,y=`Mean fold change`))+
  geom_point(size=2)+
  geom_path(linewidth=1.1)+
  geom_errorbar(aes(ymin=`CI min`, ymax=`CI max`),width=max(hrs)/10)+
  scale_x_continuous(breaks=c(0,hrs),limits=c(-2,max(hrs)+4))+
  scale_y_continuous(limits=c(0,2.2))+
  labs(y='Fold change (relative to hour 2)',title='Mean fold changes in cell populations\nafter start of starvation conditions',subtitle='Prior experiment')+
  theme_bw()

if(plotWindow) foldMeansPlot0.1


# Cumulative population doublings ####

cPDs <- list()
ind <- 1
for(c1 in 1:length(conds)){
  x <- logFoldsPlotting |> 
    filter(Hours==max(hrs),`Serum condition`==conds[c1]) |> 
    ungroup() |> 
    select(`Log fold change`) |> 
    pull()
  xbar <- meanLogFolds |> 
    filter(Hours==max(hrs),`Serum condition`==conds[c1]) |> 
    ungroup() |> 
    select(!c(Hours,`Serum condition`,plot))
  
  for(c2 in 1:length(conds)){
    if(c2==c1) test <- NA
    else{
      y <- logFoldsPlotting |> 
        filter(Hours==max(hrs),`Serum condition`==conds[c2]) |> 
        ungroup() |> 
        select(`Log fold change`) |> 
        pull()
      
      test <- wilcox.test(x,y)$p.value
    }
    
    cPDs[[ind]] <- tibble('Serum condition'=conds[c1],
                           'MLFC'=xbar$`Mean log fold change`,
                           'CI min'=xbar$`CI min`,
                           'CI max'=xbar$`CI max`,
                           'CI alpha'=xbar$`CI alpha`,
                           'Serum condition 2'=conds[c2],
                           'p-value'=test)
    ind <- ind+1
  }
}

cPDs <- bind_rows(cPDs) |> 
  nest('Significance tests'=c(`Serum condition 2`,`p-value`))

## Plotting ####

coords <- list(tibble(x=c(3,5),
                      y=c(2,2)),
               tibble(x=c(1,2),
                      y=c(1.25,1.25)))
labels <- tibble(x=c(4,1.5),
                 y=c(2.1,1.35),
                 lab=c('ns','ns'))

doublingsPlot <- cPDs |> 
  ggplot(aes(x=`Serum condition`,fill=`Serum condition`))+
  geom_col(aes(y=MLFC))+
  geom_errorbar(aes(ymin=`CI min`, ymax=`CI max`),width=0.5)+
  geom_line(data=coords[[1]],aes(x,y),inherit.aes=F)+
  geom_line(data=coords[[2]],aes(x,y),inherit.aes=F)+
  geom_text(data=labels,aes(x,y,label = lab),inherit.aes=F)+
  scale_fill_discrete(guide='none')+
  labs(y='Mean cumulative population doublings',title='Cumulative population doublings after 72 hours',
       subtitle=bquote('All comparisons significant at'~alpha==0.05~'except those marked \"ns\"'))

if(plotWindow) doublingsPlot


# Writing plots to file ####

if(overwrite){
  meanLogFolds |> 
    select(!plot) |> 
    write_csv('logFolds.csv')
  
  absPlots <- absPlot + absMeansPlot
  logFoldsPlots <- logFoldsPlot + logFoldMeansPlot
  foldsPlots <- foldsPlot + foldMeansPlot
  
  figs <- c('absPlots',
            'logFoldsPlots',
            'foldsPlots',
            'foldMeansPlot0.1',
            'doublingsPlot')
  
  fignames <- c('Absorbances',
                'Log2-fold changes',
                'Fold changes',
                'Fold changes 0.1pct',
                'Cumulative population doublings')
  
  folderpath <- paste(getwd(),'/figures/', sep='')
  dir.create(folderpath)
  
  widths <- c(rep(3000,3),900,1400)
  heights <- c(rep(1120,3),1000,1120)
  
  for (i in 1:length(figs)) {
    fig <- figs[i]
    png(filename=paste(folderpath, fignames[i], '.png', sep=''), width=widths[i], height=heights[i], res=250)
    print(get(fig))
    dev.off()
  }
}
