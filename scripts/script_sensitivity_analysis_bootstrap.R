#Loading required packages

library(here)
library(igraph)
library(dplyr)
library(ggpubr)
library(rptR)
library(piecewiseSEM)
library(lme4)
library(nlme)
library(multcompView)
library(DHARMa)
library(glmmTMB)
library(data.table)
library(ggnetwork)
library(patchwork)
library(expm)
library(visreg)
library(jtools)

# Loading functions
source(here("src", "species_metrics.R"))
source(here("src", "katz_function.R"))
source(here("src", "criteria.R"))
source(here("src", "summarySE.R"))

net_dry<-read.table(here("src", "network_dry.txt"), header=TRUE)
net_rainy<-read.table(here("src", "network_rainy.txt"), header=TRUE)
data<-read.table(here("src", "plant_dataset.txt"), header=TRUE)

inc_dry<-as.matrix(net_dry[,3:20])
inc_rainy<-as.matrix(net_rainy[,3:20])

rownames(inc_dry)<-net_dry$plant
rownames(inc_rainy)<-net_rainy$plant

#----------------------------------------------------#
#-Bootstrap confidence intervals and network metrics-#
#----------------------------------------------------#

#----------------#
#-Preparing data-#
#----------------#

incidence_dry<-inc_dry[, which(colSums(inc_dry)>0)]
incidence_rainy<-inc_rainy[, which(colSums(inc_rainy)>0)]

types_dry<-data.frame(node=c(rownames(incidence_dry), colnames(incidence_dry)), 
                      type=c(rep(FALSE, nrow(incidence_dry)), rep(TRUE, ncol(incidence_dry))))

types_rainy<-data.frame(node=c(rownames(incidence_rainy), colnames(incidence_rainy)), 
                        type=c(rep(FALSE, nrow(incidence_rainy)), rep(TRUE, ncol(incidence_rainy))))

g_dry<-graph_from_incidence_matrix(incidence_dry, weighted=TRUE)
edges_dry<-data.table(igraph::as_data_frame(g_dry)); edges_dry<-edges_dry[, .SD[rep(.I, weight)]]; edges_dry<-edges_dry[,1:2]

g_rainy<-graph_from_incidence_matrix(incidence_rainy, weighted=TRUE)
edges_rainy<-data.table(igraph::as_data_frame(g_rainy)); edges_rainy<-edges_rainy[, .SD[rep(.I, weight)]]; edges_rainy<-edges_rainy[,1:2]

#-------------------------#
#-Giant component metrics-#
#-------------------------#

net_bootstrap_sensit<-function(M, n_samples, sp_types){
  
  B<-M[sample(n_samples, replace=TRUE), ]
  g<-graph_from_data_frame(B, directed=FALSE)
  g<-simplify(g)
  V(g)$type=sp_types$type[match(V(g)$name, sp_types$node)]
  A<-get.adjacency(g, sparse=FALSE)
  A2<-get.incidence(g, sparse=FALSE)
  
  crit<-criteria(A2)
  ncomp<-components(g)$no
  giantc_prop<-max(components(g)$csize)/sum(components(g)$csize)
  
  r<-data.frame(ncomp, giantc_prop, crit, n_samples)
  
  return(r)
}

n_samples_steps<-seq(29, 244, 5)

boot_sensit_dry<-lapply(n_samples_steps, FUN=function(x){do.call(rbind, replicate(n=1000, expr=net_bootstrap_sensit(M=edges_dry, sp_types=types_dry, n_samples=x), simplify=FALSE))})
boot_sensit_rainy<-lapply(n_samples_steps, FUN=function(x){do.call(rbind, replicate(n=1000, expr=net_bootstrap_sensit(M=edges_rainy, sp_types=types_rainy, n_samples=x), simplify=FALSE))})

boot_sensit_dry<-do.call(rbind, boot_sensit_dry); boot_sensit_dry$season<-"Dry"
boot_sensit_rainy<-do.call(rbind, boot_sensit_rainy); boot_sensit_rainy$season<-"Rainy"
boot_df<-rbind(boot_sensit_rainy, boot_sensit_dry)
boot_df$n_samples<-as.numeric(boot_df$n_samples)

#--------------------------#
#-Giant component plotting-#
#--------------------------#

p2<-ggline(data=boot_df, x="n_samples", y="ncomp", add="mean_sd", color="season", plot_type="p", numeric.x.axis=TRUE)+
  xlab("Number of samples")+ylab("Number of components")+
  scale_color_manual(values=c("#DC3220", "#005AB5"))+
  theme_pubr()

p3<-ggline(data=boot_df, x="n_samples", y="giantc_prop", add="mean_sd", color="season", plot_type="p", numeric.x.axis=TRUE)+
  xlab("Number of samples")+ylab("Fraction of nodes in largest component")+
  scale_color_manual(values=c("#DC3220", "#005AB5"))+
  theme_pubr()

p4<-ggline(data=boot_df, x="n_samples", y="crit", add="mean_sd", color="season", plot_type="p", numeric.x.axis=TRUE)+
  xlab("Number of samples")+ylab("Connectivity parameter")+
  scale_color_manual(values=c("#DC3220", "#005AB5"))+
  theme_pubr()

p2+p3+p4

summarySE(data=boot_df, measurevar="crit", groupvars="season")