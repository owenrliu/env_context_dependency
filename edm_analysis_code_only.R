## ----global_options, include=F,message=F---------------------------------
library(knitr)
knitr::opts_chunk$set(echo = FALSE,warning=FALSE,message = FALSE,results = 'asis')

## ----setup, include=FALSE,message=F--------------------------------------
# Load required packages (if these are not installed on your system, use `install.packages()`)
library(tidyverse)
library(here)
library(rEDM)
library(igraph)
library(quantreg)
library(knitr)
library(kableExtra)
library(gridExtra)
library(RANN)
library(plot3D)
library(ggsci)

# plot theme
plot_theme <-   theme_minimal()+
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=14),
        axis.title=element_text(family="sans",size=14),
        axis.text=element_text(family="sans",size=8),
        strip.background = element_rect(colour="black"),
        panel.border = element_rect(color="black",fill=NA))

theme_set(plot_theme)

## ----load data, message=FALSE,echo=F-------------------------------------
source("data/data preparation scripts/produce_final_data.R")

## ----list of variables,echo=F,fig.cap="List of variables"----------------
tibble(Variable=colnames(westend),Description=c("Monitoring Transect","Monitoring Period","Laminaria Normalized Density",
                                            "Macrocystis pyrifera >1m Normalized Density","Pterygophora Normalized Density",
                                            "Strongylocentrotus purpuratus Normalized Density","Mesocentrotus franciscanus Normalized Density",
                                            "Macrocystis pyrifera <1m Normalized Density","Normalized Multivariate ENSO Index",
                                            "Normalized Pacific Decadal Oscillation","Normalized North Pacific Gyre Oscillation",
                                            "Normalized Maximum Significant Wave Height","Normalized Sea Surface Temperature")) %>%
  knitr::kable()

## ----mean density over time, echo=F,message=F----------------------------
# Species naming key for use in plotting
namekey <- tibble(dataset=c(rep("Benthic density",6),rep("Physical",5)),short=names(westend[3:13]),
                  long=c("Laminaria","Macrocystis >1m","Pterygophora","Purple urchin","Red urchin","Macrocystis <1m","Multivariate ENSO Index","Pacific Decadal Oscillation","North Pacific Gyre Oscillation","Significant Wave Height","Sea Surface Temperature"),plotting=c("L. far","M. pyr","P. cal","S. pur","M. fra","M. pyr (j)","(P) MEI","(P) PDO","(P) NPGO","(P) SWH","(P) SST"))

namekey <- bind_rows(namekey, data_frame(long="NA",dataset="NA",short="const",plotting="Constant")) %>%
  
  # define column "plotting" as a factor so we can later control plotting order on figures
  mutate(plotting=factor(plotting,levels=c("(P) MEI","(P) PDO","(P) NPGO","(P) SWH","(P) SST","M. pyr","M. pyr (j)","L. far","P. cal","S. pur","M. fra","Constant")),
         long=factor(long,levels=c("Macrocystis >1m","Macrocystis <1m","Pterygophora","Laminaria","Purple urchin","Red urchin","Multivariate ENSO Index","Pacific Decadal Oscillation","North Pacific Gyre Oscillation","Significant Wave Height","Sea Surface Temperature","Constant")))

# key relating monitoring period to its year/month
period.key <- data_frame(year=rep(1980:2014,each=12),month=rep(1:12,35),period=c(rep(NA,8),rep(1,4),rep(2:69,each=6)))
periods.first <- period.key %>% 
  distinct(period,.keep_all=T) %>%
  mutate(day=1) %>% #assume all monitoring done on first day of month
  unite(date_str,year,month,day,sep="-") %>%
  mutate(date_str=ymd(date_str))

westend.long <- westend %>%
  gather(key=spp,value=dens,-period,-site,na.rm=T) %>%
  left_join(namekey, by=c("spp"="short"))
westend.meandens <- westend.long %>%
  group_by(spp,period,long) %>% 
  summarise(meandens=mean(dens,na.rm=T),sddens=sd(dens,na.rm=T))%>%
  mutate(lower=meandens-sddens,upper=meandens+sddens)%>%
  left_join(periods.first,by="period") %>%
  ungroup()

# position adjustment for visual clarity
pd <- position_dodge(width=100)
westend.mean_dens_plot <- westend.meandens %>%
  filter(spp %in% c("mac","red","pter","ymac","purp","lam")) %>%
  mutate(type=ifelse(spp %in% c("mac","ymac","pter","lam"),"algae","urchin")) %>%
  ggplot(aes(x=date_str,y=meandens,col=long,linetype=type))+
      geom_hline(yintercept=0)+
      geom_line(lwd=2,alpha=0.8,position=pd)+
      geom_pointrange(aes(ymin=lower,ymax=upper),size=0.3,position=pd)+
      xlab("Year")+ylab("Normalized Density")+
      scale_color_manual(values=c("navyblue","gray50","darkgreen","darkcyan","mediumpurple4","darkred"),
                         guide=guide_legend(title="Species",override.aes = list(size=1)))+
      guides(linetype="none")+
      theme(text = element_text(color="black",size=16),
            axis.title=element_text(size=10),
            legend.position = c(0.15,0.7),
            legend.text = element_text(color="black",size=10),
            panel.border = element_blank())
rm(pd,westend.long)

## ----physical variables timeseries,fig.height=9,fig.width=7--------------
mei.pdo.npgo.ts.plot <-westend %>%
  ungroup()%>%
  slice(1:63)%>%
  select(period,mei,pdo,npgo)%>%
  gather("driver","value",mei:npgo,-period)%>%
  left_join(namekey, by=c("driver"="short"))%>%
  left_join(periods.first,by="period") %>%
  ggplot(aes(date_str,value,col=driver))+
  geom_line(lwd=1.5)+
  geom_hline(yintercept=0,linetype=2,col="black")+
  labs(title="Physical Drivers: Climate Modes",x="Year",y="Normalized Index Value")+
  theme(axis.title=element_text(size=10))
sst.waves.ts.plot <-westend %>%
  ungroup()%>%
  slice(1:63)%>%
  select(period,sst,waves)%>%
  gather("driver","value",sst,waves,-period)%>%
  left_join(namekey, by=c("driver"="short"))%>%
  left_join(periods.first,by="period") %>%
  ggplot(aes(date_str,value,col=driver))+
  geom_line(lwd=1.5)+
  geom_hline(yintercept=0,linetype=2,col="black")+
  labs(title="Physical Drivers: Temperature and Disturbance",x="Year",y="Normalized Index Value")+
  theme(axis.title=element_text(size=10))
lay <- rbind(1,2,3,3)
grid.arrange(mei.pdo.npgo.ts.plot,sst.waves.ts.plot,westend.mean_dens_plot,layout_matrix=lay)
rm(lay,westend.mean_dens_plot,sst.waves.ts.plot,mei.pdo.npgo.ts.plot)

## ----time series segments, include=F-------------------------------------
# segments of the time series (which rows of the data are individual transects?)
westend_segs_sites <- westend %>%
  ungroup() %>%
  mutate(ind = row_number()) %>% 
  group_by(site) %>%
  summarise(first=first(ind),last=last(ind))

# Time series segments without site identifier. We'll need this later
westend.segs <- select(westend_segs_sites,-site)

# species in the study
study_spp <- c("red","purp","lam","pter","mac","ymac")
# physical variables in the study
phys.vars <- c("mei","pdo","npgo","waves","sst")

## ----univariate simplex westend,fig.width=7,echo=F-----------------------
# List to store output of simplex projection
westend.simp.list <- list()

# Run simplex projection for each species at each site, and plot output
for(i in 1:length(study_spp)) {
  spp <- study_spp[i]
  dat <- westend %>% select(matches(spp)) %>% as.data.frame()
  out <- simplex(as.numeric(dat[,1]),lib=as.matrix(westend.segs),E=2:15,silent=T,stats_only = F) %>%
    mutate(spp=namekey$long[match(spp,namekey$short)])
  westend.simp.list[[spp]] <- out
}

# plot
bind_rows(westend.simp.list) %>%
  ggplot(aes(E,rho,color=spp))+
  geom_line(size=2)+
  facet_wrap(~spp,nrow=2,scales="free_y")+
  labs(x="Embedding Dimension (E)",y=expression(paste("Skill, ",rho)))+
  scale_x_continuous(breaks=seq(0,12,by=2))+
  scale_color_manual(values=c('Macrocystis >1m'="navyblue",'Macrocystis <1m'="gray50",'Pterygophora'="darkgreen",'Laminaria'="darkcyan",'Purple urchin'="mediumpurple4",'Red urchin'="darkred"),name="")+
  guides(color=F)

# Save best embedding dimensions
westend.bestE <- sapply(westend.simp.list,function(x) {
  temp <- x %>% filter(!is.na(rho))
  temp$E[temp$rho==max(temp$rho)]
})

rm(spp,dat,out)

# Run simplex for each physical variable
phys.simp.list <- list()
for(i in 1:length(phys.vars)) {
  ind <- phys.vars[i]
  dat <- westend %>% select(matches(ind)) %>% as.data.frame()
  out <- simplex(as.numeric(dat[,1]),lib=c(1,63),E=2:10,silent=T) %>%
    mutate(spp=namekey$long[match(ind,namekey$short)])
  phys.simp.list[[ind]] <- out
  westend.simp.list[[ind]] <- out
}

# plot
bind_rows(phys.simp.list) %>%
  ggplot(aes(E,rho))+
  geom_line(size=2)+
  facet_wrap(~spp,nrow=2,scales="free_y")+
  labs(x="Embedding Dimension (E)",y=expression(paste("Skill, ",rho)))+
  scale_x_continuous(breaks=seq(0,10,by=2))+
  guides(color=F)+
  theme(strip.text = element_text(size=8))

rm(ind,dat,out)
phys.bestE <- sapply(phys.simp.list,function(x) {
  temp <- x %>% filter(!is.na(rho))
  temp$E[temp$rho==max(temp$rho)]
})

westend.bestE <- c(westend.bestE,phys.bestE)

## ----prediction horizon test,message=F,echo=F,fig.width=7,echo=F---------
# List to hold prediction horizon results
westend.tp.list <- list()

par(mfrow=c(3,2))
for(i in 1:length(study_spp)) {
  spp <- study_spp[i]
  dat <- westend %>% select(matches(spp)) %>% as.data.frame()
  out <- simplex(as.numeric(dat[,1]),lib=as.matrix(westend.segs),silent=T,E=westend.bestE[spp],tp=1:10) %>%
    mutate(spp=namekey$long[match(spp,namekey$short)])
  westend.tp.list[[spp]] <- out
}

#Plot time horizon vs. rho
bind_rows(westend.tp.list) %>%
  ggplot(aes(tp,rho,color=spp))+
  geom_line(size=2)+
  facet_wrap(~spp,nrow=2,scales="free_y")+
  labs(x="Time to Prediction",y=expression(paste("Skill, ",rho)))+
  scale_color_manual(values=c('Macrocystis >1m'="navyblue",'Macrocystis <1m'="gray50",'Pterygophora'="darkgreen",'Laminaria'="darkcyan",'Purple urchin'="mediumpurple4",'Red urchin'="darkred"),name="")+
  scale_x_continuous(breaks=seq(0,10,by=2))+
  guides(color=F)

rm(spp,dat,out)

## ----univariate s_maps westend,fig.width=7,echo=F------------------------
westend.smap.list <- list() # list to store output

par(mfrow=c(3,2))
for(i in 1:length(study_spp)) {
  spp <- study_spp[i]
  dat <- westend %>% select(matches(spp)) %>% as.data.frame()
  out <- s_map(as.numeric(dat[,1]),lib=as.matrix(westend.segs),E=westend.bestE[spp],silent=T) %>%
    mutate(spp=namekey$long[match(spp,namekey$short)])
  westend.smap.list[[spp]] <- out
}

bind_rows(westend.smap.list) %>%
  ggplot(aes(theta,rho,color=spp))+
  geom_line(size=2)+
  facet_wrap(~spp,nrow=2,scales="free_y")+
  labs(x=expression(paste("Nonlinearity (",theta,")")),y=expression(paste("Skill, ",rho)))+
  scale_color_manual(values=c('Macrocystis >1m'="navyblue",'Macrocystis <1m'="gray50",'Pterygophora'="darkgreen",'Laminaria'="darkcyan",'Purple urchin'="mediumpurple4",'Red urchin'="darkred"),name="")+
  scale_x_continuous(breaks=seq(0,8,by=2))+
  guides(color=F)

rm(spp,dat,out)

## ----ccm example---------------------------------------------------------
tempE <- westend.bestE['red']
temp <- suppressWarnings(ccm(westend,lib=as.matrix(westend.segs),pred=as.matrix(westend.segs),E=tempE,lib_column= 'mac',target_column = 'red',lib_sizes = c(10,25,50,75,100,125,150,200,300,400,500),num_samples=100,replace=T,silent=T,RNGseed = 41389))

mac_xmap_red <- temp %>%
  group_by(lib_size)%>%
  summarise(rhomean=mean(rho,na.rm=T),upper=quantile(rho, 0.9),lower=quantile(rho, 0.1))%>%
  mutate(lower=pmax(0,lower)) %>% 
  ungroup()%>%
  ggplot(aes(lib_size,rhomean))+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3,fill="red")+
  geom_line(color="darkorchid3")+
  labs(x="Library Size",y=expression(paste(rho, " (predictive skill)")),title="Macrocystis xmap Red Urchin")

mac_xmap_red

## ----ccm all species and physical vars, echo=F,warning=F,fig.width=7-----
# Rows hold predicted variables, columns hold forcing variables. There are more columns than rows since the physical variables are included as potential forcing factors.
n_col <- dim(westend)[2]-2
n_row <- dim(westend)[2]-7
col_names <- colnames(westend)[3:(n_col+2)]
row_names <- colnames(westend)[3:(n_row+2)]
westend.xmap_mat <- array(NA,dim=c(n_row,n_col),dimnames=list(row_names,col_names))

## CCM causation criterion 1: cross-map skill greater than zero
# matrix to store a bootstrapped p-value, measuring the probability that a given xmap is greater than zero (calculated as 1 minus the number of positive results for rho divided by the number of iterations)
westend.p1.mat <- array(NA,dim=c(n_row,n_col),dimnames=list(row_names,col_names))

## CCM causation criterion 2: evidence for convergence
# similarly, matrix to store a bootstrapped p-value, this time a t-test value between library size 10 and library size 500, to see if the rho at large library is significantly greater than the rho at small library size (i.e., looking for convergence).
westend.p2.mat <- array(NA,dim=c(n_row,n_col),dimnames=list(row_names,col_names))

# if both p1 and p2 are positive, indicate overall significant causal signal
westend.ptot.mat <- array(NA,dim=c(n_row,n_col),dimnames=list(row_names,col_names))

## Run CCM for each combination of variables
for(i in 1:n_row) {
  for(j in 1:n_col) {
    if(i != j) {
      # remember, we use the best embedding dimension for the target variable (the variable we're cross-mapping to, i.e. the putative forcing variable)
      tempE=westend.bestE[col_names[j]]
      temp <- suppressWarnings(ccm(westend,lib=as.matrix(westend.segs),pred=as.matrix(westend.segs),E=tempE,lib_column= 2+i,target_column = 2+j,lib_sizes = c(10,500),num_samples=100,replace=T,silent=T,RNGseed = 41389))
      
     # mean rho at library size 500
      rhomeans <- temp %>% ccm_means()
      westend.xmap_mat[i,j] <- rhomeans$rho[rhomeans$lib_size==500]
      
      # first p-value (is cross map skill greater than zero? at library size 500)
      p1 <- temp %>% filter(lib_size==500) %>% 
        mutate(pos=ifelse(rho>0,1,0)) %>%
        summarise(p=(1-sum(pos)/n()))
      westend.p1.mat[i,j] <- as.numeric(p1)
      
      # second p-value (rho at lib-size 500 greater than rho at lib-size 10? By more than 0.1?)
      p2 <- t.test(temp$rho[temp$lib_size==10],temp$rho[temp$lib_size==500])$p.value
      if((rhomeans$rho[rhomeans$lib_size==500]-rhomeans$rho[rhomeans$lib_size==10])<0.1) p2 <- 1
      westend.p2.mat[i,j] <- as.numeric(p2)
      
      # overall significance (both p1 and p2 signficant at alpha 0.05)
      westend.ptot.mat[i,j] <- ifelse(p1<0.05 & p2<0.05,1,0)
    }
  }
}
rm(n_col,n_row,col_names,row_names,p1,p2,tempE,temp,rhomeans)

# keep only signficant cross-mappings
westend.xmap_mat <- westend.xmap_mat*westend.ptot.mat %>% as.data.frame()
westend.xmap_rast <- westend.xmap_mat %>% 
  mutate(predictee=row.names(westend.xmap_mat)) %>% 
  gather(key=predictor,value=rho,-predictee)

# If rho is zero, replace with NA (no significant causal signal)
westend.xmap_rast$rho[westend.xmap_rast$rho==0] <- NA

# Names for plotting to distinguish biological and physical variables
westend.xmap_rast <- westend.xmap_rast %>% 
  mutate(force.name=namekey$plotting[match(predictor,namekey$short)]) %>%
  mutate(pred.name=namekey$plotting[match(predictee,namekey$short)])

## plot
westend.xmap_all_plot <- ggplot(westend.xmap_rast,aes(x=force.name,y=pred.name,fill=rho)) +
  geom_raster() +
  scale_fill_gradient(low = "#9AFF9A", high = "#548B54", space = "Lab", na.value = "grey50", 
      guide = "colourbar",limits=c(0,0.8),breaks=c(0.2,0.4,0.6,0.8),name=expression(paste(rho, "(skill)"))) +
  geom_text(aes(label=round(rho,2)),family="Arial")+
  ggtitle("Kelp Forest Convergent Cross Mapping") +
  xlab("Predictor (Forcing Variable)") +
  ylab("Predicted Variable") +
  theme(text=element_text(color="black"),
        plot.background = element_rect(fill="white"),
        axis.text=element_text(color="black"),
    axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.1),
    panel.border = element_blank())

westend.xmap_all_plot

## ----network,fig.height=7,fig.width=7------------------------------------
## i graph
# connections matrix for the network based on the CCM results above
connect.mat <- t(westend.ptot.mat) %>% 
  as.data.frame() %>%
  mutate(mei=rep(0,nrow(.)),pdo=rep(0,nrow(.)),npgo=rep(0,nrow(.)),waves=rep(0,nrow(.)),sst=rep(0,nrow(.)))%>%
  as.matrix()

# edges (arrows) weighted based on CCM rho value
edge.weight.mat <- t(westend.xmap_mat)%>%
  as.data.frame() %>%
  mutate(mei=rep(0,nrow(.)),pdo=rep(0,nrow(.)),npgo=rep(0,nrow(.)),waves=rep(0,nrow(.)),sst=rep(0,nrow(.)))%>%
  as.matrix()
edge.weight.mat[is.na(edge.weight.mat)]<-0

names <- map_chr(colnames(connect.mat), function(x) as.character(namekey$plotting[match(x,namekey$short)]))
names[6] <- "M. pyr\n(j)"
names [7:11] <- c("MEI","PDO","NPGO","SWH","SST")

# initiate the igraph
g2 <- graph.adjacency(connect.mat)

# Vertex labels
V(g2)$name <-names

# vector of node colors for plotting
vert.cols <- c(rep("#238B45",3),rep("#BEBADA",2),"#238B45",rep("#8DD3C7",5))

# weighted edges, with width and color
edge.weights <- as.numeric(t(edge.weight.mat))[t(edge.weight.mat) != 0]*5
E(g2)$width <- edge.weights
col.1 <- adjustcolor("#99D8C9", alpha=0.8)
col.2 <- adjustcolor("#006D2C", alpha=0.8)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.pal <- edge.pal(100)
edge.cols<-round((edge.weights-min(edge.weights))/(max(edge.weights)-min(edge.weights))*100) +1
edge.cols <- edge.pal[edge.cols]

# some edges must be curved to look okay in the graph
curves <- rep(0,length(edge.weights))
curves[c(1,4,18,21)] <-0.4
curves[c(19)] <-0.3
curves[c(6,8,10,16,20)] <- 0.2

#plot to file if uncommented
plot(g2,asp=0.8, vertex.label.color="black",vertex.label.family="sans", vertex.label.cex=0.8, edge.arrow.size=0.8,edge.arrow.width=0.8,vertex.label.font=2, vertex.shape="circle", edge.curved=curves,margin=0,vertex.size=20,rescale=T,edge.lty=1,edge.color=edge.cols,layout=layout_in_circle,vertex.color=vert.cols,vertex.frame.color="black",main="")

## ----causal vars fxn,message=F,echo=F------------------------------------
# Pull out causal vars from CCM for each species
westend.causal.vars <- purrr::map(c("mac","purp","lam","ymac","red","pter"), function(x) {
  westend.xmap_rast %>%
    filter(predictee==x,!is.na(rho))->out
  out <- c(x,out$predictor)
  out
})

names(westend.causal.vars)=c("mac","purp","lam","ymac","red","pter")

## ----3d plot,fig.height=4,fig.width=4,dpi=300----------------------------
x <- westend$lam
y <- westend$mac
z <- westend$purp

seg.start <- as.numeric(westend.segs[1,1])
seg.end <- as.numeric(westend.segs[1,2])

#smooth
t <- 1:length(seg.start:seg.end)
tt<-seq(1,length(seg.start:seg.end),len=500)

xsmooth<-splinefun(t, x[seg.start:seg.end])(tt)
ysmooth<-splinefun(t, y[seg.start:seg.end])(tt)
zsmooth<-splinefun(t, z[seg.start:seg.end])(tt)

xtemp <- x[seg.start:seg.end]
ytemp <- y[seg.start:seg.end]
ztemp <- z[seg.start:seg.end]

# plot
par(mai=c(0.1,0.1,0.1,0.1))
scatter3D(xsmooth,ysmooth,zsmooth, col="gray30",
          type="l",phi=40,theta=55,bty="b",xlab = "Laminaria", ylab = "Macrocystis", zlab = "Purple Urchin")
points3D(xtemp,ytemp,ztemp,col="black",pch=19,cex=0.5,add=T)

for(i in 2:5) {
  seg.start <- as.numeric(westend.segs[i,1])
  seg.end <- as.numeric(westend.segs[i,2])
  
  #smooth
  t <- 1:length(seg.start:seg.end)
  tt<-seq(1,length(seg.start:seg.end),len=500)
  xsmooth<-splinefun(t, x[seg.start:seg.end])(tt)
  ysmooth<-splinefun(t, y[seg.start:seg.end])(tt)
  zsmooth<-splinefun(t, z[seg.start:seg.end])(tt)
  
  xtemp <- x[seg.start:seg.end]
  ytemp <- y[seg.start:seg.end]
  ztemp <- z[seg.start:seg.end]
  
  # plot
  
  scatter3D(xsmooth,ysmooth,zsmooth, col="gray30",
            type="l",bty="b",add=T)
  points3D(xtemp,ytemp,ztemp,col="black",pch=19,cex=0.5,add=T)
}

scatter3D(x[46],y[46],z[46],col="red",pch=19,cex=2,add=T)

#manually find nearest neighbors to point 46
temp <- westend %>% select(lam,mac,purp) %>% slice(1:315) %>% mutate(id=row_number()) %>%
  filter(!is.na(lam))
tempnn <-nn2(temp[,1:3],k=4)$nn.idx %>% as.data.frame()
names(tempnn) <- c("nn1","nn2","nn3","nn4")
temp <- bind_cols(temp,tempnn)

scatter3D(x[109],y[109],z[109],col="darkgreen",pch=19,cex=2,add=T)
scatter3D(x[235],y[235],z[235],col="darkgreen",pch=19,cex=2,add=T)
scatter3D(x[169],y[169],z[169],col="darkgreen",pch=19,cex=2,add=T)

scatter3D(x[110],y[110],z[110],col="blue",pch=19,cex=2,add=T)
scatter3D(x[236],y[236],z[236],col="blue",pch=19,cex=2,add=T)
scatter3D(x[170],y[170],z[170],col="blue",pch=19,cex=2,add=T)
#

rm(xsmooth,ysmooth,zsmooth,xtemp,ytemp,ztemp,t,tt,seg.start,seg.end,x,y,z,temp)

## ----general smap model fxn, include=F-----------------------------------
# Generalized function to build a multivariate S-map model. Function takes the data, a data frame of the segments (i.e., the denoting the different within-site swaths denoting breaks in the time series), the target species (the one for which we are building a model), and a character vector of the other, causal variables. The function returns the fitted model and the coefficients (i.e., the species interaction partials), the model statistics, and a box plot of the interactions coefficients for that species across all the data.

smap_multi <- function(sitedat,sitesegs,species,causalvars) {
  
  full.mod <- suppressWarnings(block_lnlp(sitedat,lib=as.matrix(sitesegs),pred=as.matrix(sitesegs),columns=causalvars,target_column = species,theta = c(0, 1e-04, 3e-04, 0.001,0.003, 0.01, 0.1, 0.5, 1, 2, 4, 6,8), num_neighbors=0,method="s-map",silent=T))

  opttheta <- full.mod$theta[full.mod$rho==max(full.mod$rho)]
  
  best.mod <- suppressWarnings(block_lnlp(sitedat,lib=as.matrix(sitesegs),pred=as.matrix(sitesegs),columns=causalvars,target_column = species,theta = opttheta, num_neighbors=0,method="s-map",save_smap_coefficients=T,silent=T))
  
  # gather and organize model output
  t <- best.mod$model_output[[1]]$time
  
  coeff.names <- c(causalvars,"const")
  # Save interaction coefficients
  coeffs <- best.mod$smap_coefficients[[1]]
  names(coeffs) <- coeff.names
  
  # # Model output-- time, observations, predictions. match time indicator with observations (one period ahead)
  coeffs <- coeffs %>% mutate(time=t)%>%complete(time=full_seq(x=c(1,630),1)) %>% mutate(period=rep(1:63,10))
  
  # Save coefficient variances (sensu Deyle 2016; we have a nxn variance-covariance matrix for each estimated point, where n is the total number of causal variables for whom coefficients are estimated, plus the constant term)
  # Variances for individual coefficients are the diagonals of each covariance matrix
  coeffs.variances <- best.mod$smap_coefficient_covariances[[1]] %>%
    purrr::map_dfr(function(x) {
      # if no estimate, no variance associated with it (NULL)
      if(is.null(x)) x <- matrix(NaN,nrow=length(coeff.names),ncol=length(coeff.names))
      # pull out diagonal elements, add variable names
      out <- t(diag(x)) %>% as.data.frame()
      names(out) <- coeff.names
      out
      })%>%
    
    # add period identifier and overall time stamp
    # mutate(time=t) %>% complete(time=full_seq(c(1,630),1)) %>% mutate(period=rep(1:63,10))
    mutate(time=t)%>%complete(time=full_seq(x=c(1,630),1)) %>% mutate(period=rep(1:63,10))
  
  # long form for plotting, including variance
  coeffs.long <- coeffs %>%
    gather(key=effect,value=value,-period,-time,na.rm=T) %>%
    left_join(namekey, by=c("effect"="short"))
  
  # long form variances
  coeffs.variances.long <- coeffs.variances %>%
    gather(key=effect,value=variance,-period,-time,na.rm=T) %>%
    left_join(namekey,by=c("effect"="short"))
  
  coeffs.long <- coeffs.long %>%
    left_join(coeffs.variances.long,by=c("time","period","effect","dataset","long","plotting")) %>%
    # add standard deviation and +/- 95% CI (1.96*SD)
    mutate(sd=sqrt(variance),upper=value+1.96*sd,lower=value-1.96*sd)
  
  # box plot of interactions
  sppname <- species
  int.box <- coeffs.long %>% filter(effect!='const')%>%
    ggplot(aes(x=plotting,y=value))+
    geom_boxplot() +
    geom_hline(yintercept=0,linetype=2,col="gray40")+
    xlab("Interacting Variable") +
    ylab("Coefficient") +
    ggtitle(paste("Interactive Effects on Dynamics of ",sppname))
  
  ## Outputs (to go in a list)

  outlist <- list(model=best.mod,coeffs=coeffs,coeffs.variances=coeffs.variances,coeffs.long= coeffs.long,plot=int.box)
  return(outlist)
}

## ----apply multivariate models, include=F,warnings=F---------------------
# Apply to the West End data
westend.multi.mods <- purrr::map(list("mac","purp","lam","ymac","red","pter"),function(x) {
  smap_multi(westend,westend.segs,x,westend.causal.vars[[x]])
})
names(westend.multi.mods) <- c("mac","purp","lam","ymac","red","pter")

## ----model perf,echo=F,results='asis'------------------------------------
# collect model performance statistics
westend.model.perf <- tibble(Species=c("Macrocystis","Purple Urchin","Laminaria","Young Macrocystis","Red Urchin","Pterygophora"),
                             Predictors=purrr::map(westend.causal.vars,.f=function(x){
                               y=namekey$long[match(x,namekey$short)];paste(y,sep=",",collapse=", ")}),
                             Theta = map_dbl(westend.multi.mods,function(x) x$model$theta),
                             Rho = map_dbl(westend.multi.mods,function(x) x$model$rho),
                             MAE=map_dbl(westend.multi.mods,function(x) x$model$mae),
                             mod=c("mac","purp","lam","ymac","red","pter"))

westend.model.perf %>% 
  select(-mod) %>%
  kable("html",col.names = c("Modeled Species","Predictors","$\\theta$","$\\rho$","MAE"),digits=3,escape=F) %>%
  kable_styling(full_width = F) %>%
  column_spec(1,bold=T)%>%
  column_spec(2, width = "30em")

## ----all model effects join, echo=F,warning=F----------------------------
# Join all coefficients
all.mod.coeffs.long <- purrr::map(names(westend.multi.mods),function(x){
  westend.multi.mods[[x]][["coeffs.long"]] %>% mutate(mod=x)
  }
  )%>%
  bind_rows() %>%
  mutate(mod_plotting_name = namekey$long[match(mod,namekey$short)])


# assign interaction types: competition, herbivory (effect of urchins on algae), consumption (effect of algae on urchin), or growth
algaes<-c("lam","mac","pter","ymac")
urchins <- c("red","purp")
combinations <- expand.grid(study_spp,study_spp) %>% 
  rename(mod=Var1,effect=Var2) %>%
  select(mod,effect)%>%
  mutate(type = case_when(
    mod ==   effect ~ "intraspecies",
    mod == "mac" & effect == "ymac" ~ "intraspecies",
    mod == "ymac" & effect == "mac" ~ "intraspecies",
    mod %in% algaes & effect %in% algaes ~ "algal competition",
    mod %in% urchins & effect %in% urchins ~ "urchin competition",
    mod %in% algaes & effect %in% urchins ~ "herbivory",
    mod %in% urchins & effect %in% algaes ~ "consumption",
    TRUE ~ "NA"
    ),
    type=factor(type,levels=c("algal competition","urchin competition","herbivory","consumption","intraspecies")),
    spp_type = case_when(
      mod %in% algaes ~ "algae",
      mod %in% urchins ~ "urchin",
      TRUE ~ "none"
    )
  ) %>%
  mutate_at(vars(mod,effect),as.character)

# add rho label for facet plotting
mod_rho_labels <- westend.model.perf %>%
  mutate(label1="rho==",Rho_round=round(Rho,2))%>%
  unite(label,label1,Rho_round,sep="",remove=F)%>%
  select(mod,label)

all.mod.coeffs.long %>%
  #add the type of interaction
  left_join(combinations,by=c("mod","effect"))%>%
  left_join(mod_rho_labels,by="mod") -> all.mod.coeffs.long

## ----boxplot of species interactions, echo=F,warning=F,message=F,fig.height=8,fig.width=7.5----
#
effect_categories <- all.mod.coeffs.long %>%
  #filter for just species (no physical vars yet)
  filter(effect %in% study_spp) %>%
  select(effect,plotting) %>%
  distinct(effect,plotting) %>% .$plotting %>% as.character()

# boxplot of all species interactions
all.mod.coeffs.plot.box <- all.mod.coeffs.long %>%
  #filter for just species (no physical vars yet)
  filter(effect %in% study_spp) %>%
  #plot
  ggplot(aes(x=plotting,y=value,fill=type))+
  geom_boxplot(outlier.alpha=0.4,outlier.size=1,width=0.8,position = "identity")+
  stat_summary(fun.y=mean, colour="white", geom="point", 
               shape=17, size=2,show.legend = FALSE)+
  scale_fill_manual(values=c("darkorange1","orangered2","darkorchid1","palegreen3","skyblue1"),name="")+
  scale_x_discrete(breaks=effect_categories,limits=effect_categories,labels=effect_categories,drop=FALSE)+
  geom_hline(yintercept=0,linetype=2)+
  geom_text(aes(x=1.2,y=1.8,label=label),size=4, parse=T,check_overlap=TRUE,show.legend = FALSE)+
  ylim(-2,2)+
  coord_flip()+
  labs(y="Estimated Interaction Coefficient",x="Interacting Species",title="Estimated Species Interactions in Multivariate Models")+
  facet_grid(mod_plotting_name~.,scales="free",space="free")

all.mod.coeffs.plot.box

## ----phys coeffs,echo=F,message=F,warning=F,fig.height=8,fig.width=8-----
# with boxplot
all.phys.coeffs.plot <- all.mod.coeffs.long %>%
  filter(effect %in% phys.vars)%>%
  mutate(label = case_when(
    mod == "mac" ~ "Macrocystis >1m",
    mod == "ymac" ~ "Macrocystis <1m",
    mod == "lam" ~ "Laminaria",
    mod == "pter" ~ "Pterygophora",
    mod == "red" ~ "Red urchin",
    mod == "purp" ~ "Purple urchin")
  ) %>%
  ggplot(aes(effect,value))+
  scale_x_discrete(labels=c(
    "mei"="MEI",
    "npgo"="NPGO",
    "pdo"="PDO",
    "waves"="SWH",
    "sst"="SST")
    )+
  geom_boxplot(fill="dodgerblue2",outlier.alpha=0.4,position = "identity")+
  stat_summary(fun.y=mean, colour="white", geom="point", 
               shape=17, size=2,show.legend = FALSE)+
  geom_hline(yintercept=0,linetype=2)+
  ylim(-2,2)+
  coord_flip()+
  labs(y="Estimated Interaction Coefficient",x="Physical Forcing",title="Estimated Effect of Physical Variables in Multivariate Models")+
  facet_grid(mod_plotting_name~.)+
  geom_text(aes(x=1.5,y=2,label=label),hjust=1,family="Arial",fontface="italic",check_overlap = TRUE,size=4,show.legend = FALSE)+
  theme(strip.background = element_blank(),
        strip.text = element_blank())
all.phys.coeffs.plot

## ----model spp_only,warning=F--------------------------------------------
# remove physical variables from predictors then run models
causal.vars.spp <- purrr::map(westend.causal.vars,function(x) setdiff(x,phys.vars))
multi.mods.spp.only <- purrr::map(list("mac","purp","lam","ymac","red","pter"),function(x) {
  smap_multi(westend,westend.segs,x,causal.vars.spp[[x]])
})
names(multi.mods.spp.only) <- c("mac","purp","lam","ymac","red","pter")

## ----fig.height=8,fig.width=7.5------------------------------------------
# collect model performance statistics
spp.only.model.perf <- tibble(Species=c("Macrocystis","Purple Urchin","Laminaria","Young Macrocystis","Red Urchin","Pterygophora"),
                             Predictors=purrr::map(causal.vars.spp,.f=function(x){
                               y=namekey$long[match(x,namekey$short)];paste(y,sep=",",collapse=", ")}),
                             Rho = map_dbl(multi.mods.spp.only,function(x) x$model$rho),
                             MAE=map_dbl(multi.mods.spp.only,function(x) x$model$mae),
                             mod=c("mac","purp","lam","ymac","red","pter"))
spp.only.model.perf %>% select(-mod) %>% 
  kable("html",col.names = c("Modeled Species","Predictors","$\\rho$","MAE"),digits=3,escape=F) %>%
  kable_styling(full_width = F) %>%
  column_spec(1,bold=T)%>%
  column_spec(2, width = "30em")

# Join all coefficients
spp.only.mod.coeffs <- purrr::map(names(multi.mods.spp.only),function(x){
  multi.mods.spp.only[[x]][["coeffs.long"]] %>% mutate(mod=x)
  }
  )%>%
  bind_rows() %>%
  mutate(mod_plotting_name = namekey$long[match(mod,namekey$short)])


# assign interaction types: competition, herbivory (effect of urchins on algae), consumption (effect of algae on urchin), or growth
algaes<-c("lam","mac","pter","ymac")
urchins <- c("red","purp")
combinations <- expand.grid(study_spp,study_spp) %>% 
  rename(mod=Var1,effect=Var2) %>%
  select(mod,effect)%>%
  mutate(type = case_when(
    mod ==   effect ~ "intraspecies",
    mod == "mac" & effect == "ymac" ~ "intraspecies",
    mod == "ymac" & effect == "mac" ~ "intraspecies",
    mod %in% algaes & effect %in% algaes ~ "algal competition",
    mod %in% urchins & effect %in% urchins ~ "urchin competition",
    mod %in% algaes & effect %in% urchins ~ "herbivory",
    mod %in% urchins & effect %in% algaes ~ "consumption",
    TRUE ~ "NA"
    ),
    type=factor(type,levels=c("algal competition","urchin competition","herbivory","consumption","intraspecies")),
    spp_type = case_when(
      mod %in% algaes ~ "algae",
      mod %in% urchins ~ "urchin",
      TRUE ~ "none"
    )
  ) %>%
  mutate_at(vars(mod,effect),as.character)

# add rho label for facet plotting
mod_rho_labels <- spp.only.model.perf %>%
  mutate(label1="rho==",Rho_round=round(Rho,2))%>%
  unite(label,label1,Rho_round,sep="",remove=F)%>%
  select(mod,label)

spp.only.mod.coeffs %>%
  #add the type of interaction
  left_join(combinations,by=c("mod","effect"))%>%
  left_join(mod_rho_labels,by="mod") -> spp.only.mod.coeffs

# boxplot of all species interactions
spp.only.mod.coeffs.plot.box <- spp.only.mod.coeffs %>%
   filter(effect %in% study_spp) %>%
  #plot
  ggplot(aes(x=plotting,y=value,fill=type))+
  geom_boxplot(outlier.alpha=0.4,outlier.size=1,width=0.8,position = "identity")+
  stat_summary(fun.y=mean, colour="white", geom="point", 
               shape=17, size=2,show.legend = FALSE)+
  scale_fill_manual(values=c("darkorange1","orangered2","darkorchid1","palegreen3","skyblue1"),name="")+
  scale_x_discrete(breaks=effect_categories,limits=effect_categories,labels=effect_categories,drop=FALSE)+
  geom_hline(yintercept=0,linetype=2)+
  geom_text(aes(x=1.2,y=1.8,label=label),size=4, parse=T,check_overlap=TRUE,show.legend = FALSE)+
  ylim(-2,2)+
  coord_flip()+
  labs(y="Estimated Interaction Coefficient",x="Interacting Species",title="Estimated Species Interactions in Multivariate Models\n(Physical variables removed)")+
  facet_grid(mod_plotting_name~.,scales="free_y",space="free")

spp.only.mod.coeffs.plot.box

## ----mean interactions, fig.height=8,fig.width=7.5-----------------------
# join all original density variables and physical variables to interaction data

spp.only.coeffs.with.phys<- westend %>%
  select(one_of(c(phys.vars,"period"))) %>%
  distinct() %>%
  right_join(filter(spp.only.mod.coeffs,!is.na(type)),by="period")%>%
  gather("external_forcing","ext_value",one_of(phys.vars))


## ----interactions density plot by type,fig.height=7,fig.width=6----------
# Density plot by interaction type

interaction_distribution_by_type <- spp.only.mod.coeffs %>%
  filter(!is.na(type)) %>%
  group_by(type) %>%
  mutate(mean_int=mean(value,na.rm=T))%>%
  ungroup()%>%
  ggplot(aes(value,..density..,fill=type))+
  geom_density(alpha=0.8,color=NA)+
  # geom_density_ridges(alpha=0.8,scale=1.4)+
  # stat_density_ridges(geom="density_ridges_gradient",calc_ecdf='true',quantiles=4,alpha=0.8,scale=1)+
  # scale_fill_viridis_d(name="quantiles",alpha=0.8)+
  geom_vline(xintercept=0,linetype=2)+
  geom_vline(aes(xintercept=mean_int),color='black')+
  scale_fill_manual(values=c("darkorange1","orangered2","darkorchid1","palegreen3","skyblue1"),name="")+
  scale_color_manual(values=c("darkorange1","orangered2","darkorchid1","palegreen3","skyblue1"),name="")+
  guides(fill=F,color=F)+
  xlim(-1.5,1.5)+
  labs(x="Interaction strength",y="Probability density",title="Distribution of Interaction Strengths by Type")+
  facet_grid(rows=vars(type))

interaction_distribution_by_type

## ----mac pter,fig.width=6,fig.height=6-----------------------------------
mac_effect_on_pter <-spp.only.coeffs.with.phys %>%
  ungroup()%>%
  filter(mod=="pter",effect=="mac")%>%
  mutate(phys_high_low=case_when(
    ext_value > 1 ~ "high",
    ext_value < -1 ~ "low")) %>%
  filter(!is.na(phys_high_low))%>%
  left_join(select(westend,period,one_of(phys.vars),red),by="period")%>%
  group_by(type,external_forcing,phys_high_low)%>%
  mutate(mean_int=mean(value,na.rm=T))%>%
  distinct()%>%
  ungroup()

# physical forcings labeller
phys.labeller <- c(
  mei= "MEI",
  npgo= "NPGO",
  pdo= "PDO",
  sst="SST",
  waves="SWH"
)

mac_effect_on_pter_plot <-mac_effect_on_pter %>%
  select(-red) %>% distinct() %>%
ggplot(aes(value,..density..,fill=phys_high_low))+
  geom_density(col=NA,alpha=0.6)+
  geom_vline(xintercept=0,linetype=2)+
  geom_vline(aes(xintercept=mean_int,color=phys_high_low),size=1.5)+
  scale_fill_manual(values=c("darkorange","dodgerblue4"),name="Physical\nForcing Value",labels=c("high","low"))+
  scale_color_manual(values=c("darkorange","dodgerblue4"))+
  guides(color="none")+
  labs(x="Effect of Macrocystis on Pterygophora",y="Probability Density", title="Varying Interaction Strength across Environmental Gradients\n(n=520 interactions for each box)")+
  facet_grid(rows=vars(external_forcing),labeller=labeller(external_forcing=phys.labeller))+
  theme(strip.text=element_text(size=10))

mac_effect_on_pter_plot

## ----algal comp context,fig.width=8,fig.height=6-------------------------
# pull out all estimate algal competition coefficients and establish whether they're neutral
alg_comp <- spp.only.mod.coeffs %>%
  filter(type %in% c("algal competition","herbivory"))%>%
  mutate(neutral=ifelse(lower<0 & upper>0,1,0))%>%
  left_join(select(westend,period,one_of(phys.vars)),by="period") %>% distinct()%>%
  select(mod,effect,long,plotting,value,neutral,one_of(phys.vars))

alg_comp_contexts_simple <- alg_comp%>%
  mutate_at(vars(one_of(phys.vars)),function(x) {case_when(
    x < -1 ~ "Neg",
    x < 1 & x > -1 ~ "Norm",
    x > 1 ~ "Pos"
  )
    })%>%
  gather("driver","context",mei:sst)%>%
  mutate(context=factor(context,levels=c("Neg",'Norm','Pos')))%>%
  distinct()%>%
  group_by(long,plotting,effect,driver,context)%>%
  summarise(n=n(),perc_negative=sum(neutral==0 & value<0)/n()*-100,perc_positive=sum(neutral==0 & value>0)/n()*100,
            perc_neutral=sum(neutral)/n()*100, perc_nonneutral=sum(neutral==0)/n()*100)%>%
  ungroup()%>%
  filter(perc_neutral<100,!is.na(context),effect!='ymac')%>%
  gather("comp_or_mutualism","prevalance",perc_negative:perc_positive)%>%
  ungroup()%>%
  arrange(effect)

competition_herbivory_context_plot <- alg_comp_contexts_simple%>%
  filter(effect != "pter")%>%
  ggplot(aes(context,prevalance,fill=comp_or_mutualism))+
  geom_col(color=NA)+
  scale_y_continuous(breaks=seq(-60,60,by=10),labels=c(60,50,40,30,20,10,0,10,20,30,40,50,60))+
  scale_fill_npg(labels=c("Negative","Positive"),name="Interaction Direction")+
  labs(x="Level of Physical Driver",y="Percent of Estimated Interactions")+
  facet_grid(long~driver,labeller = labeller(driver=phys.labeller),scales="free")
competition_herbivory_context_plot

