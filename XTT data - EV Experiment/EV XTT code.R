# Setup ####

library(stringr)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(pacman)
library(knitr)

overwrite <- F # Change to TRUE to overwrite plots in /../figures/ with new plots
plotWindow <- F # Change to TRUE to show plots in RStudio window

lkw <- 1.5 # legend key width

# Data Wrangling ####

absorb <- read_csv('all.csv') |> drop_na() |> 
  pivot_longer(!Hours,names_to = 'Treatment',values_to = 'Abs') |> 
  mutate(Dose=Treatment,.before=3) |> 
  mutate(Treatment=str_replace_all(Treatment,"\\s\\d+","")) |> 
  mutate(Dose=str_replace_all(Dose,"^.*?(?=\\d{2,3}$)","")) |> 
  mutate(Dose=str_replace_all(Dose,"PBS","0")) |> 
  mutate(Treatment=str_replace_all(Treatment,"PBS","No EVs")) |> 
  mutate(Treatment=str_replace_all(Treatment,"Ctrl","Control")) |> 
  mutate(Treatment=str_replace_all(Treatment,"Pt3","PKD2 truncating")) |> 
  mutate(Treatment=str_replace_all(Treatment,"Pt9","PKD1 truncating")) |> 
  mutate(Treatment=factor(Treatment,levels=c('No EVs','Control','PKD1 truncating','PKD2 truncating'))) |> 
  mutate(Dose=factor(Dose,levels=c('0','10','50','100')))

absorbOrig <- absorb

absorb0s <- absorb |> 
  filter(Abs<0) |> 
  mutate(Abs=0)

absorb <- absorb |> 
  filter(Abs>=0) |> 
  bind_rows(absorb0s)

conds <-unique(absorb$Treatment)
lineCols <- palette.colors(length(conds),palette='Dark2')
names(lineCols) <- conds

hrs <- absorb |> 
  select(Hours) |> 
  pull() |> 
  unique()

doses <- absorb |> 
  select(Dose) |> 
  pull() |> 
  unique()


# Absorbance plots ####

absorbPlotting <- absorb |> 
  rowwise() |> 
  mutate(plot=any(Hours!=0,Treatment=='No EVs'))


## All absorbances ####

absPlot <- absorbPlotting |> 
  ggplot(aes(x=Hours,y=`Abs`,color=Treatment))+
  geom_point(aes(shape=plot))+
  scale_x_continuous(breaks=hrs,limits=c(0,max(hrs)+2))+
  scale_y_continuous(limits=c(0,max(absorbPlotting$Abs)))+
  geom_smooth(aes(linetype=Dose),method='lm',formula=y~x,se=F)+  
  scale_linetype_manual(values=c('solid','dotted','dotdash','longdash'),name=bquote('Dose'~(mu*g/mL)))+ 
  scale_shape_manual(values=c(NA,16),guide='none')+
  scale_color_discrete(type=lineCols)+
  labs(y='Corrected absorbance',title='All corrected absorbances')+
  theme_bw()+
  theme(legend.key.width=unit(lkw,'cm'))

if(overwrite) {
  absPlot <- absPlot+scale_color_discrete(type=lineCols,guide='none')+scale_linetype_discrete(guide='none')
}

if(plotWindow) absPlot


## Mean absorbances ####

meanFunk <- function(x,i){mean(x[i])} # for boot
se <- function(data,times=2000){
  b <- boot::boot(data,meanFunk,times)
  return(sd(b$t))
}

means <- absorbPlotting |> 
  group_by(Hours,Treatment,Dose) |> 
  summarise('Mean Abs'=mean(Abs),'Mean Abs SE'=se(Abs)) |> 
  rowwise() |> 
  mutate(plot=any(Hours!=0,Treatment=='No EVs'))

absMeansPlot <- means |>
  ggplot(aes(x=Hours,y=`Mean Abs`,color=Treatment))+ 
  geom_point(size=2)+
  geom_path(aes(linetype=Dose),linewidth=1.1)+
  geom_errorbar(aes(ymin=`Mean Abs`-`Mean Abs SE`, ymax=`Mean Abs`+`Mean Abs SE`),width=max(hrs)/10)+
  scale_linetype_manual(values=c('solid','dotted','dotdash','longdash'),name=bquote('Dose'~(mu*g/mL)))+ 
  scale_color_discrete(type=lineCols)+
  labs(y='Corrected absorbance',title='Mean corrected absorbances')+
  scale_x_continuous(breaks=hrs)+
  scale_y_continuous(limits=c(0,max(absorbPlotting$Abs)))+
  theme_bw()+
  theme(legend.key.width=unit(lkw,'cm'))

if(plotWindow) absMeansPlot


doseNamer <- function(variable, value){
  return(
    list('0'=bquote('0'~mu*g/mL),
         '10'=bquote('10'~mu*g/mL),
         '50'=bquote('50'~mu*g/mL),
         '100'=bquote('100'~mu*g/mL))[value]
  )
}

absMeansPlotSepTreat <- means |>
  ggplot(aes(x=Hours,y=`Mean Abs`,color=Treatment))+ 
  facet_grid(cols=vars(Treatment))+
  geom_point(size=2)+
  geom_path(aes(linetype=Dose),linewidth=1.1)+
  geom_errorbar(aes(ymin=`Mean Abs`-`Mean Abs SE`, ymax=`Mean Abs`+`Mean Abs SE`),width=max(hrs)/10)+
  scale_linetype_manual(values=c('solid','dotted','dotdash','longdash'),name=bquote('Dose'~(mu*g/mL)))+ 
  scale_color_discrete(type=lineCols,guide='none')+
  labs(y='Corrected absorbance',title='Mean corrected absorbances',subtitle='Separated by treatment')+
  scale_x_continuous(breaks=hrs)+
  scale_y_continuous(limits=c(0,max(absorbPlotting$Abs)))+
  theme_bw()+
  theme(legend.key.width=unit(lkw,'cm'))

if(plotWindow) absMeansPlotSepTreat

absMeansPlotSepDose <- means |>
  ggplot(aes(x=Hours,y=`Mean Abs`,color=Treatment))+ 
  facet_grid(cols=vars(Dose),labeller='doseNamer')+
  geom_point(size=2)+
  geom_path(aes(linetype=Dose),linewidth=1.1)+
  geom_errorbar(aes(ymin=`Mean Abs`-`Mean Abs SE`, ymax=`Mean Abs`+`Mean Abs SE`),width=max(hrs)/10)+
  scale_linetype_manual(values=c('solid','dotted','dotdash','longdash'),name=bquote('Dose'~(mu*g/mL)),guide='none')+ 
  scale_color_discrete(type=lineCols)+
  labs(y='Corrected absorbance',title='Mean corrected absorbances',subtitle='Separated by dose')+
  scale_x_continuous(breaks=hrs)+
  scale_y_continuous(limits=c(0,max(absorbPlotting$Abs)))+
  theme_bw()+
  theme(legend.key.width=unit(lkw,'cm'))

if(plotWindow) absMeansPlotSepDose


# Fold changes plots ####

init <- absorb |> 
  filter(Hours==0 & Treatment=='No EVs') |> 
  select(Abs) |> 
  pull() |> 
  mean()

folds <- absorb |> 
  mutate('Fold change'=Abs/init) |> 
  select(!Abs)

logFoldsPlotting <- folds |> 
  mutate('Log fold change'=log2(`Fold change`)) |> 
  select(!`Fold change`) |> 
  rowwise() |> 
  mutate(plot=any(Hours!=0,Treatment=='No EVs'))

meanLogFolds <- logFoldsPlotting |> 
  group_by(Hours,Treatment,Dose) |> 
  summarise('Mean log fold change'=mean(`Log fold change`),
            'errBar min'=`Mean log fold change` - se(`Log fold change`),
            'errBar max'=`Mean log fold change` + se(`Log fold change`)) |> 
  rowwise() |> 
  mutate(plot=any(Hours!=0,Treatment=='No EVs'))

## All fold changes ####

foldsPlot <- logFoldsPlotting |> 
  ggplot(aes(x=Hours,y=2^`Log fold change`,color=Treatment))+
  geom_point(aes(shape=plot))+
  scale_x_continuous(breaks=hrs,limits=c(-2,max(hrs)+4))+
  scale_y_continuous(limits=c(min(2^logFoldsPlotting$`Log fold change`),
                              max(2^logFoldsPlotting$`Log fold change`)+0.2))+
  geom_smooth(aes(linetype=Dose),method='lm',formula=y~x,se=F)+ 
  scale_linetype_manual(values=c('solid','dotted','dotdash','longdash'),name=bquote('Dose'~(mu*g/mL)))+ 
  scale_shape_manual(values=c(NA,16),guide='none')+
  scale_color_discrete(type=lineCols)+
  labs(y='Fold change',title='All fold changes')+
  theme_bw()+
  theme(legend.key.width=unit(lkw,'cm'))

if(overwrite){
  foldsPlot <- foldsPlot+scale_color_discrete(type=lineCols,guide='none')+scale_linetype_discrete(guide='none')
}

if(plotWindow) foldsPlot


## Mean fold changes ####

lineCols2 <- lineCols
names(lineCols2)[2] <- 'No mutation'

foldMeansPlot <- meanLogFolds |> 
  mutate(Treatment=str_replace_all(Treatment,'Control','No mutation')) |> 
  mutate(Treatment=factor(Treatment,levels=c('No EVs','No mutation','PKD1 truncating','PKD2 truncating'))) |> 
  ggplot(aes(x=Hours,y=2^`Mean log fold change`,color=Treatment))+ 
  geom_point(size=2)+
  geom_path(aes(linetype=Dose),linewidth=1.1)+
  geom_errorbar(aes(ymin=2^`errBar min`, ymax=2^`errBar max`),width=max(hrs)/10)+
  scale_linetype_manual(values=c('solid','dotted','dotdash','longdash'),name=bquote('Dose'~(mu*g/mL)))+ 
  scale_color_discrete(type=lineCols2)+
  scale_x_continuous(breaks=hrs,limits=c(-2,max(hrs)+4))+
  scale_y_continuous(limits=c(min(2^logFoldsPlotting$`Log fold change`),
                              max(2^logFoldsPlotting$`Log fold change`)+0.2))+
  labs(y='Fold change',title='Mean fold changes',color='EV mutation')+
  theme_bw()+
  theme(legend.key.width=unit(lkw,'cm'))

if(plotWindow) foldMeansPlot


foldMeansPlotSepTreat <- meanLogFolds |> 
  mutate(Treatment=str_replace_all(Treatment,'Control','No mutation')) |> 
  mutate(Treatment=factor(Treatment,levels=c('No EVs','No mutation','PKD1 truncating','PKD2 truncating'))) |> 
  ggplot(aes(x=Hours,y=2^`Mean log fold change`,color=Treatment))+ 
  facet_grid(cols=vars(Treatment))+
  geom_point(size=2)+
  geom_path(aes(linetype=Dose),linewidth=1.1)+
  geom_errorbar(aes(ymin=2^`errBar min`, ymax=2^`errBar max`),width=max(hrs)/10)+
  scale_linetype_manual(values=c('solid','dotted','dotdash','longdash'),name=bquote('Dose'~(mu*g/mL)))+ 
  scale_color_discrete(type=lineCols2,guide='none')+
  scale_x_continuous(breaks=hrs,limits=c(-4,max(hrs)+4))+
  labs(y='Fold change',title='Mean fold changes',subtitle='Separated by treatment',color='EV mutation')+
  theme_bw()+
  theme(legend.key.width=unit(lkw,'cm'))

if(plotWindow) foldMeansPlotSepTreat


foldMeansPlotSepDose <- meanLogFolds |> 
  mutate(Treatment=str_replace_all(Treatment,'Control','No mutation')) |> 
  mutate(Treatment=factor(Treatment,levels=c('No EVs','No mutation','PKD1 truncating','PKD2 truncating'))) |> 
  ggplot(aes(x=Hours,y=2^`Mean log fold change`,color=Treatment))+ 
  facet_grid(cols=vars(Dose),labeller='doseNamer')+
  geom_point(size=2)+
  geom_path(aes(linetype=Dose),linewidth=1.1)+
  geom_errorbar(aes(ymin=2^`errBar min`, ymax=2^`errBar max`),width=max(hrs)/10)+
  scale_linetype_manual(values=c('solid','dotted','dotdash','longdash'),name=bquote('Dose'~(mu*g/mL)),guide='none')+ 
  scale_color_discrete(type=lineCols2)+
  scale_x_continuous(breaks=hrs,limits=c(-4,max(hrs)+4))+
  labs(y='Fold change',title='Mean fold changes',subtitle='Separated by dose',color='EV mutation')+
  theme_bw()

if(plotWindow) foldMeansPlotSepDose

foldMeansPlotSepDoseSingle <- meanLogFolds |> 
  mutate(Treatment=str_replace_all(Treatment,'Control','No mutation')) |> 
  mutate(Treatment=factor(Treatment,levels=c('No EVs','No mutation','PKD1 truncating','PKD2 truncating'))) |> 
  ggplot(aes(x=Hours,y=2^`Mean log fold change`,color=Treatment))+ 
  facet_grid(cols=vars(Dose),labeller='doseNamer')+
  geom_point(size=2)+
  geom_path(linewidth=1.1)+
  geom_errorbar(aes(ymin=2^`errBar min`, ymax=2^`errBar max`),width=max(hrs)/10)+
  scale_x_continuous(breaks=hrs,limits=c(-4,max(hrs)+4))+
  scale_color_discrete(type=lineCols2)+
  labs(y='Fold change (relative to hour 0)',title='Mean fold changes of cell populations after EV treatment',subtitle='Separated by dose',color='EV mutation')+
  theme_bw()

if(plotWindow) foldMeansPlotSepDoseSingle

foldMeansPlotSepDoseSingleExp <- meanLogFolds |> 
  mutate(Treatment=str_replace_all(Treatment,'Control','No mutation')) |> 
  mutate(Treatment=factor(Treatment,levels=c('No EVs','No mutation','PKD1 truncating','PKD2 truncating'))) |> 
  ggplot(aes(x=Hours,y=2^`Mean log fold change`,color=Treatment))+ 
  facet_grid(cols=vars(Dose),labeller='doseNamer')+
  geom_point(size=2)+
  geom_path(linewidth=1.1)+
  geom_errorbar(aes(ymin=2^`errBar min`, ymax=2^`errBar max`),width=max(hrs)/10)+
  scale_x_continuous(breaks=hrs,limits=c(-4,max(hrs)+4))+
  scale_y_continuous(limits=c(0,2.2))+
  scale_color_discrete(type=lineCols2)+
  labs(y='Fold change (relative to hour 0)',title='Mean fold changes of cell populations after EV treatment',subtitle='Separated by dose',color='EV mutation')+
  theme_bw()

if(plotWindow) foldMeansPlotSepDoseSingleExp


# Relative-survival (to No EVs) ####

meanAbsNoEVs <- absorb |>
  filter(Treatment=='No EVs') |> 
  select(!Dose) |> 
  group_by(Hours,Treatment) |> 
  summarise(`Mean Abs`=mean(Abs))

relAbs <- list()
i <- 1
for(h in 1:length(hrs)){
  baseline <- meanAbsNoEVs |> 
    ungroup(Hours) |> 
    filter(Hours==hrs[h]) |> 
    select(`Mean Abs`) |> 
    pull()
  if(baseline==0) next
  for(c in 2:length(conds)){
    for(d in 2:length(doses)){
      relAbs[[i]] <- absorb |> 
        filter(Hours==hrs[h] & Treatment==conds[c] & Dose==doses[d]) |> 
        mutate(RelPercent=100*Abs/baseline) |> 
        select(!Abs)
      i <- i+1
    }
  }
}
relAbs <- bind_rows(relAbs)

meanRelAbs <- relAbs |> 
  group_by(Hours,Treatment,Dose) |> 
  summarise('Mean Abs'=mean(RelPercent),'Mean Abs SE'=se(RelPercent))

## Comparing treatments ####

sigTests <- list()
i <- 1
for(h in 2:3){
  for(d in 2:4){
    ctrl <- relAbs |> 
      filter(Hours==hrs[h] & Dose==doses[d] & Treatment=='Control') |> 
      ungroup() |> 
      select(RelPercent) |> 
      pull()
    
    for(c in 3:4){
      treated <- relAbs |> 
        filter(Hours==hrs[h] & Dose==doses[d] & Treatment==conds[c]) |> 
        ungroup() |> 
        select(RelPercent) |> 
        pull()
      
      test <- wilcox.test(treated,ctrl,alternative='greater')$p.value
      
      sigTests[[i]] <- tibble(Hours=hrs[h],
                              Treatment=conds[c],
                              Dose=doses[d],
                              mean=mean(treated),
                              'p-value'=test)
      i <- i+1
    }
  }
}

sigTests <- bind_rows(sigTests)

sigLabs <- sigTests |> 
  filter(`p-value`<0.05 & `p-value`>=0.005) |> 
  mutate(lab='*')
sigLabs <- sigTests |> 
  filter(`p-value`<0.005) |> 
  mutate(lab='**') |> 
  bind_rows(sigLabs)
sigTests <- sigTests |> 
  filter(`p-value`>0.05) |> 
  mutate(lab='') |> 
  bind_rows(sigLabs)

labels <- sigTests |> 
  mutate(y=(mean-3)/100) |> 
  mutate(x=Hours+7) |> 
  select(!c(Hours,mean,`p-value`))



sigTestsP2 <- list()
i <- 1
for(h in 2:3){
  for(d in 2:4){
    PKD1 <- relAbs |> 
      filter(Hours==hrs[h] & Dose==doses[d] & Treatment=='PKD1 truncating') |> 
      ungroup() |> 
      select(RelPercent) |> 
      pull()
    
    PKD2 <- relAbs |> 
      filter(Hours==hrs[h] & Dose==doses[d] & Treatment=='PKD2 truncating') |> 
      ungroup() |> 
      select(RelPercent) |> 
      pull()
    
    test <- wilcox.test(PKD2,PKD1)$p.value
    
    sigTestsP2[[i]] <- tibble(Hours=hrs[h],
                            Dose=doses[d],
                            mean=mean(PKD2),
                            'p-value'=test)
    i <- i+1
  }
}

sigTestsP2 <- bind_rows(sigTestsP2)

sigLabsP2 <- sigTestsP2 |> 
  filter(`p-value`<0.05 & `p-value`>=0.005) |> 
  mutate(lab='†')
sigLabsP2 <- sigTestsP2 |> 
  filter(`p-value`<0.005) |> 
  mutate(lab='‡') |> 
  bind_rows(sigLabsP2)
sigTestsP2 <- sigTestsP2 |> 
  filter(`p-value`>0.05) |> 
  mutate(lab='') |> 
  bind_rows(sigLabsP2)

labelsP2 <- sigTestsP2 |> 
  mutate(y=(mean+3)/100) |> 
  mutate(x=Hours-7) |> 
  select(!c(Hours,mean,`p-value`))



meanRelAbsPlotSepDose <- meanRelAbs |>
  mutate(Treatment=str_replace_all(Treatment,'Control','No mutation')) |> 
  mutate(Treatment=factor(Treatment,levels=c('No mutation','PKD1 truncating','PKD2 truncating'))) |> 
  ggplot(aes(x=Hours,y=`Mean Abs`/100,color=Treatment))+ 
  facet_grid(cols=vars(Dose),labeller = 'doseNamer')+
  geom_point(size=2)+
  geom_path(aes(linewidth=Dose))+
  geom_errorbar(aes(ymin=(`Mean Abs`-`Mean Abs SE`)/100, ymax=(`Mean Abs`+`Mean Abs SE`)/100),width=max(hrs)/10)+
  scale_linewidth_manual(values=c(0.2,1,2),name=bquote('Dose'~(mu*g/mL)),guide='none')+ 
  labs(y='Surviving cells (relative to untreated)',title='Relative-survival of cells after EV treatment',color='EV mutation',subtitle='Separated by dose')+
  geom_text(data=labels,mapping=aes(x,y,label = lab),inherit.aes=F)+
  geom_text(data=labelsP2,mapping=aes(x,y,label = lab),inherit.aes=F)+
  scale_x_continuous(breaks=hrs,limits=c(0,60))+
  scale_y_continuous(labels=scales::percent)+
  scale_color_manual(values=lineCols2)+
  theme_bw()

if(plotWindow) meanRelAbsPlotSepDose


## Comparing doses ####

sigTestsDose <- list()
i <- 1
for(h in 2:3){
  for(c in 2:4){
    for(ld in 2:3){
      lowDose <- relAbs |> 
        filter(Hours==hrs[h] & Dose==doses[ld] & Treatment==conds[c]) |> 
        ungroup() |> 
        select(RelPercent) |> 
        pull()
      
      for(hd in 3:4){
        if (ld==hd) next
        highDose <- relAbs |> 
          filter(Hours==hrs[h] & Dose==doses[hd] & Treatment==conds[c]) |> 
          ungroup() |> 
          select(RelPercent) |> 
          pull()
        
        test <- wilcox.test(highDose,lowDose,alternative='greater')$p.value
        
        sigTestsDose[[i]] <- tibble(Hours=hrs[h],
                                Treatment=conds[c],
                                'High Dose'=doses[hd],
                                'Low dose'=doses[ld],
                                mean=mean(highDose),
                                'p-value'=test)
        i <- i+1
      }
    }
  }
}

sigTestsDose <- bind_rows(sigTestsDose)

sigLabsDose <- sigTestsDose |> 
  filter(`p-value`<0.05 & `p-value`>=0.005 & `Low dose`=='10') |> 
  mutate(lab='¤')
sigLabsDose <- sigTestsDose |> 
  filter(`p-value`<0.005 & `Low dose`=='10') |> 
  mutate(lab='¤¤') |> 
  bind_rows(sigLabsDose)
sigLabsDose <- sigTestsDose |> 
  filter(`p-value`<0.05 & `p-value`>=0.005 & `Low dose`=='50') |> 
  mutate(lab='#') |> 
  bind_rows(sigLabsDose)
sigLabsDose <- sigTestsDose |> 
  filter(`p-value`<0.005 & `Low dose`=='50') |> 
  mutate(lab='##') |> 
  bind_rows(sigLabsDose)
sigTestsDose <- sigTestsDose |> 
  filter(`p-value`>0.05) |> 
  mutate(lab='') |> 
  bind_rows(sigLabsDose)

sigTestsDose <- sigTestsDose |> 
  mutate(Treatment=str_replace_all(Treatment,'Control','No mutation')) |> 
  mutate(Treatment=factor(Treatment,levels=c('No mutation','PKD1 truncating','PKD2 truncating')))

labelsDose1 <- sigTestsDose |> 
  filter(`Low dose`=='10') |> 
  mutate(y=(mean)/100) |> 
  mutate(x=Hours+8)

labelsDose2 <- sigTestsDose |> 
  filter(`Low dose`=='50') |> 
  mutate(y=(mean)/100) |> 
  mutate(x=Hours-8)

meanRelAbsPlotSepTreat <- meanRelAbs |>
  mutate(Treatment=str_replace_all(Treatment,'Control','No mutation')) |> 
  mutate(Treatment=factor(Treatment,levels=c('No mutation','PKD1 truncating','PKD2 truncating'))) |> 
  ggplot(aes(x=Hours,y=`Mean Abs`/100,color=Treatment))+ 
  facet_grid(cols=vars(Treatment))+
  geom_point(size=2)+
  geom_path(aes(linewidth=Dose))+
  geom_errorbar(aes(ymin=(`Mean Abs`-`Mean Abs SE`)/100, ymax=(`Mean Abs`+`Mean Abs SE`)/100),width=max(hrs)/10)+
  scale_linewidth_manual(values=c(0.2,1,2),name=bquote('Dose'~(mu*g/mL)))+ 
  labs(y='Surviving cells (relative to untreated)',title='Relative-survival of cells after EV treatment',subtitle='Separated by EV mutation')+
  geom_text(data=labelsDose1,mapping=aes(x,y,label = lab),inherit.aes=F)+
  geom_text(data=labelsDose2,mapping=aes(x,y,label = lab),size=2.7,inherit.aes=F)+
  scale_x_continuous(breaks=hrs,limits=c(0,60))+
  scale_y_continuous(labels=scales::percent)+
  scale_color_manual(values=lineCols2,guide='none')+
  theme_bw()

if(plotWindow) meanRelAbsPlotSepTreat

# Writing plots to file ####

if(overwrite){
  absPlots <- absPlot + absMeansPlot
  foldsPlots <- foldsPlot + foldMeansPlot
  absPlotsSep <- absMeansPlotSepTreat / absMeansPlotSepDose
  foldsPlotsSep <- foldMeansPlotSepTreat / foldMeansPlotSepDose
  meanRelAbsPlots <- meanRelAbsPlotSepDose / meanRelAbsPlotSepTreat
  
  figs <- c('absPlots',
            'foldsPlots',
            'absPlotsSep',
            'foldsPlotsSep',
            'foldMeansPlotSepDoseSingle',
            'foldMeansPlotSepDoseSingleExp',
            'meanRelAbsPlotSepDose',
            'meanRelAbsPlotSepTreat',
            'meanRelAbsPlots')
  
  fignames <- c('Absorbances',
                'Fold changes',
                'Absorbances separated',
                'Fold changes separated',
                'Fold changes by dose',
                'Fold changes by dose expanded',
                'Population relative to untreated by dose',
                'Population relative to untreated by treatment',
                'Population relative to untreated')
  
  folderpath <- paste(getwd(),'/figures/', sep='')
  dir.create(folderpath)
  
  widths <- c(rep(3000,2),rep(3000,4),rep(1700,3))
  heights <- c(rep(1120,2),rep(2240,2),rep(1000,4),2000)
  
  for (i in 1:length(figs)) {
    fig <- figs[i]
    png(filename=paste(folderpath, fignames[i], '.png', sep=''), width=widths[i], height=heights[i], res=250)
    print(get(fig))
    dev.off()
  }
}

# Package Citations ####

if(overwrite){
  write_bib(file="Bibliography of packages.bib")
  
  packages <- tibble('Package name' = character(),
                                  Version = character(),
                                  Maintainer = character())
  
  for (pkg in p_loaded()){
    packages <- packages %>%
      add_row(
      'Package name' = pkg,
      Version = as.character(packageVersion(pkg)),
      Maintainer = maintainer(pkg)
    )
  }
  
  write_file(kable(packages,'latex'),'packages LaTeX table.txt')
}
