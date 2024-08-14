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

dry<-species_metrics(network=inc_dry, normalized=FALSE)[1:47,]
rainy<-species_metrics(network=inc_rainy, normalized=FALSE)[1:20,]

d_dry<-data.frame(plant=dry$sample, prob=dry$degree_bin/sum(dry$degree_bin), w_prob=dry$degree_w/sum(dry$degree_w), abund=dry$degree_w, total_abund=sum(dry$degree_w), total_rich=sum(dry$degree_bin))
d_rainy<-data.frame(plant=rainy$sample, prob=rainy$degree_bin/sum(rainy$degree_bin), w_prob=rainy$degree_w/sum(rainy$degree_w), abund=rainy$degree_w, total_abund=sum(rainy$degree_w), total_rich=sum(rainy$degree_bin))

d_dry$season<-"dry"
d_rainy$season<-"rainy"

df<-rbind(d_dry, d_rainy)

df_results<-left_join(df, data)

df_results$season<-ifelse(df_results$season=="dry", 2, 1)
df_results$herbivory<-log((df_results$herbivory/(1-df_results$herbivory)))

#---------------------------#
#-Fitting individual models-#
#---------------------------#

h<-lm(herbivory~compositional_pd+structural_pd, data=df_results)
h2<-lm(herbivory~season, data=df_results)
pw<-glm(w_prob~structural_pd+compositional_pd+season, family="binomial", weights=total_abund, data=df_results)
s_pd1<-lm(structural_pd~season, data=df_results)
s_pd2<-lm(structural_pd~herbivory, data=df_results)
c_pd1<-lm(compositional_pd~season, data=df_results)
c_pd2<-lm(compositional_pd~herbivory, data=df_results)

#-------------------------------------------------#
#-Using the package piecewiseSEM to fit SEM model-#
#-------------------------------------------------#

#---------------#
#-Main text SEM-#
#---------------#

sem1<-psem(
  
  h, pw, s_pd1, c_pd1,
  structural_pd %~~% compositional_pd,
  
  data=df_results
  
  
)

summary(sem1, standardize="scale")
plot(sem1)

#----------------------#
#-Alternative SEM (SI)-#
#----------------------#

sem2<-psem(
  
  h2, pw, s_pd2, c_pd2,
  structural_pd %~~% compositional_pd,
  
  data=df_results
  
)

summary(sem2)
plot(sem2)

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

net_bootstrap_gc<-function(M, sp_types, nboot){
  
  B<-M[sample(nrow(M), replace=TRUE), ]
  g<-graph_from_data_frame(B, directed=FALSE)
  g<-simplify(g)
  V(g)$type=sp_types$type[match(V(g)$name, sp_types$node)]
  A<-get.adjacency(g, sparse=FALSE)
  A2<-get.incidence(g, sparse=FALSE)
  
  crit<-criteria(A2)
  ncomp<-components(g)$no
  giantc_prop<-max(components(g)$csize)/sum(components(g)$csize)
  
  r<-data.frame(ncomp, giantc_prop, crit, nboot)
  
  return(r)
}

boot_dry<-lapply(1:10000, net_bootstrap_gc, M=edges_dry, sp_types=types_dry)
boot_rainy<-lapply(1:10000, net_bootstrap_gc, M=edges_rainy, sp_types=types_rainy)

boot_dry<-do.call(rbind, boot_dry); boot_dry$season<-"Dry"
boot_rainy<-do.call(rbind, boot_rainy); boot_rainy$season<-"Rainy"

boot_df<-rbind(boot_rainy, boot_dry)

#--------------------------#
#-Giant component plotting-#
#--------------------------#

g_r<-graph_from_incidence_matrix(incidence_rainy, weighted=TRUE)
g_d<-graph_from_incidence_matrix(incidence_dry, weighted=TRUE)

gr_df<-ggnetwork(g_r)
gr_df$season<-"Rainy season"
gd_df<-ggnetwork(g_d)
gd_df$season<-"Dry season"

net_df<-rbind(gd_df, gr_df)

data$season<-ifelse(data$season=="dry", "Dry season", "Rainy season")

data_struct<-summarySE(data=data, measurevar="structural_pd", groupvars="season")
data_comp<-summarySE(data=data, measurevar="compositional_pd", groupvars="season")

p1<-ggplot(data=net_df, aes(x=x,y=y, xend=xend, yend=yend))+
  geom_edges(aes(size=weight, alpha=weight), color="gray40", curvature=0.1, show.legend=FALSE)+
  geom_nodes(aes(fill=type), size=2, shape=21, show.legend=FALSE)+
  scale_size_continuous(range=c(0.6, 1.5))+
  scale_alpha_continuous(range=c(0.5, 0.7))+
  scale_fill_manual(values=c("palegreen3", "royalblue3"))+
  facet_wrap(~season)+
  theme_facet()

p2<-ggline(data=boot_df, x="season", y="ncomp", add="mean_sd", color="season", plot_type="p")+xlab("Season")+ylab("Number of components")+
  scale_color_manual(values=c("#DC3220", "#005AB5"))+
  theme_pubr()
p3<-ggline(data=boot_df, x="season", y="giantc_prop", add="mean_sd", color="season", plot_type="p")+xlab("Season")+ylab("Fraction of nodes in largest component")+
  scale_color_manual(values=c("#DC3220", "#005AB5"))+
  theme_pubr()
p4<-ggline(data=boot_df, x="season", y="crit", add="mean_sd", color="season", plot_type="p")+xlab("Season")+ylab("Connectivity parameter")+
  scale_color_manual(values=c("#DC3220", "#005AB5"))+
  theme_pubr()

p1 / (p2+p3+p4) + plot_annotation(tag_levels="a")

summarySE(data=boot_df, measurevar="crit", groupvars="season")

#-------------------------#
#-Network pathway metrics-#
#-------------------------#

net_bootstrap_path<-function(M, nboot){
  
  B<-M[sample(nrow(M), replace=TRUE), ]
  g<-graph_from_data_frame(B, directed=FALSE)
  g<-simplify(g)
  A<-get.adjacency(g, sparse=FALSE)
  
  A.eigen<-eigen(A)$values #Calculating eigenvalues
  A.maxeigen<-ifelse(is.complex(A.eigen), max(Re(A.eigen[which(Im(A.eigen)==0)])), max(A.eigen)) # Largest eigenvalue
  
  alfa <- 1/(A.maxeigen+0.1) #Calculating alfa
  
  Q<-alfa*A #Multiplying A matrix with alfa
  
  I<-diag(1, nrow=nrow(Q), ncol=ncol(Q)) #Getting the identity matrix of I
  
  T_mat<-solve(I-Q) #Obtaining analytical solution to Katz centrality
  
  diag(T_mat)<-0
  
  F_mat<-T_mat;F_mat[A==1]<-0
  
  ind_contrib<-sum(F_mat)/sum(T_mat)
  
  return(data.frame(total_pathways=sum(T_mat), ind_pathways=sum(F_mat), ind_contrib, nboot))
}

path_dry<-lapply(1:10000, net_bootstrap_path, M=edges_dry)
path_rainy<-lapply(1:10000, net_bootstrap_path, M=edges_rainy)

path_dry<-do.call(rbind, path_dry); path_dry$season<-"Dry"
path_rainy<-do.call(rbind, path_rainy); path_rainy$season<-"Rainy"

path_df<-rbind(path_rainy, path_dry)
summarySE(data=path_df, measurevar="ind_contrib", groupvars="season")

phist<-ggplot(data=path_df, aes(group=season))+
  geom_histogram(aes(x=ind_contrib, y = stat(count) / sum(count), color=season, fill=season), position="identity")+
  scale_fill_manual(values=c("#DC3220", "#005AB5"))+
  scale_color_manual(values=c("firebrick3", "royalblue4"))+
  labs(color="Season", fill="Season")+
  xlab("Relative contribution of indirect pathways")+ylab("Frequency")+
  theme_pubr()

phist

#--------------------------#
#-Network pathway plotting-#
#--------------------------#

A_rainy<-graph_from_incidence_matrix(incidence_rainy, directed=FALSE)
A_rainy<-get.adjacency(A_rainy, sparse=FALSE)

A_dry<-graph_from_incidence_matrix(incidence_dry, directed=FALSE)
A_dry<-get.adjacency(A_dry, sparse=FALSE)

max(shortest.paths(g_rainy))

D_1<-graph_from_adjacency_matrix(A_dry, weighted=TRUE, mode=c("undirected")); layd<-layout_nicely(D_1)
V(D_1)$type<-types_dry$type[match(V(D_1)$name, types_dry$node)]
D_1<- ggnetwork(D_1, layout=layd); D_1$season<-"Dry season"; D_1$p_length<-"Direct pathways"
D_3<- A_dry %^% 3; diag(D_3)<-0; D_3<-graph_from_adjacency_matrix(D_3, weighted=TRUE, mode=c("undirected"))
V(D_3)$type<-types_dry$type[match(V(D_3)$name, types_dry$node)]
D_3<- ggnetwork(D_3, layout=layd); D_3$season<-"Dry season"; D_3$p_length<-"Length 3 pathways"
D_5<- A_dry %^% 5; diag(D_5)<-0; D_5<-graph_from_adjacency_matrix(D_5, weighted=TRUE, mode=c("undirected"))
V(D_5)$type<-types_dry$type[match(V(D_5)$name, types_dry$node)]
D_5<- ggnetwork(D_5, layout=layd); D_5$season<-"Dry season"; D_5$p_length<-"Length 5 pathways"

R_1<-graph_from_adjacency_matrix(A_rainy, weighted=TRUE, mode=c("undirected")); layr<-layout_nicely(R_1)
V(R_1)$type<-types_rainy$type[match(V(R_1)$name, types_rainy$node)]
R_1<- ggnetwork(R_1, layout=layr); R_1$season<-"Rainy season"; R_1$p_length<-"Direct pathways"
R_3<- A_rainy %^% 3; diag(R_3)<-0; R_3<-graph_from_adjacency_matrix(R_3, weighted=TRUE, mode=c("undirected"))
V(R_3)$type<-types_rainy$type[match(V(R_3)$name, types_rainy$node)]
R_3<- ggnetwork(R_3, layout=layr); R_3$season<-"Rainy season"; R_3$p_length<-"Length 3 pathways"
R_5<- A_rainy %^% 5; diag(R_5)<-0; R_5<-graph_from_adjacency_matrix(R_5, weighted=TRUE, mode=c("undirected"))
V(R_5)$type<-types_rainy$type[match(V(R_5)$name, types_rainy$node)]
R_5<- ggnetwork(R_5, layout=layr); R_5$season<-"Rainy season"; R_5$p_length<-"Length 5 pathways"

path_networks<-rbind(R_1, R_3, R_5, D_1, D_3, D_5)

plot_paths<-ggplot(data=path_networks, aes(x=x,y=y, xend=xend, yend=yend))+
  geom_edges(aes(size=log(weight), alpha=weight), curvature=0.1, show.legend=TRUE, color="gray10")+
  geom_nodes(aes(fill=type), size=2, shape=21, show.legend=FALSE)+
  scale_size_continuous(range=c(0.5, 0.8))+
  scale_alpha_continuous(range=c(0.4, 0.5))+
  scale_fill_manual(values=c("palegreen3", "royalblue3"))+
  facet_grid(season~p_length)+
  theme_facet()

phist/plot_paths+plot_annotation(tag_levels="a")
