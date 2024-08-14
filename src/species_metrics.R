species_metrics<-function(network, normalized=FALSE){

#---------------------------------------------------------------------#  
#-Creating an igraph network and getting the network adjacency matrix-# 
#---------------------------------------------------------------------#

# Weighted network  
net_weighted<-graph_from_incidence_matrix(network, weighted=TRUE)
# Binary network
net_binary<-graph_from_incidence_matrix(network, weighted=NULL)

# Retrieving weighted adjacency matrix
adj.weight<-get.adjacency(net_weighted, attr="weight", sparse=FALSE)
# Retrieving binary adjacency matrix
adj.binary<-get.adjacency(net_binary, sparse=FALSE)

#-------------#
#-Node degree-#
#-------------#

# Node degree correspond to the summation of the number of interactions/weights of a given species
# It can be normalized by dividing a given species degree by the maximum possible degree in the network
# For bipartite networks in which there is no intraspecific interactions, the maximum degree is the number of species in the other bipartite set.

deg_bin<-rowSums(adj.binary) #Non-normalized
deg_w<-rowSums(adj.weight) #Non-normalized

if(normalized==TRUE){deg_bin<-deg_bin/(nrow(network)*ncol(network)); deg_w<-deg_w/sum(deg_w)}

#---------------------------#
#-Node closeness centrality-#
#---------------------------#

# Closeness centrality corresponds to the sum of the reciprocal of the shortest paths lengths (sum of 1/d)
# The function distances() from igraph calculates the shortest paths lengths between all pairs of nodes
# Optionally the function distances() accepts a vector of link weights and use these as the cost of traveling through a given path
# When higher link weights corresponds to an easier path to follow (e.g. abundance, interaction strength, etc)
# we need to use the inverse of the weights as costs such that pathways with larger weights are easier to follow
# Closeness centrality can be normalized by dividing by the number of nodes minus one (N-1)

spaths_bin<-1/distances(net_binary)
diag(spaths_bin)<-0
spaths_w<-1/distances(net_weighted, weights=(1/E(net_weighted)$weight))
diag(spaths_w)<-0

closeness_bin<-rowSums(spaths_bin)
closeness_w<-rowSums(spaths_w)

if(normalized==TRUE){closeness_bin<-closeness_bin/sum(closeness_bin); closeness_w<-closeness_w/sum(closeness_w)}

#-----------------#
#-Katz centrality-#
#-----------------#

# Katz centrality corresponds to the sum of pathways of all lenghts that connects to a given node, not only the shortest ones
# It can be calculated analytically by inverting the matrix K = (I-c*A), in which c is a constant that is lower than
# 1/maximum eigenvalue of A. In an undirected network we can calculate Katz centrality through both the sum of the rows
# or the sum of the columns of matrix K. However, in a directed network the row or columns sums can mean the outgoing
# pathways or the ingoing pathways depending on how the adjacency matrix is structured. For instance, if the direction of links
# in the adjacency matrix is an effect that goes out of the node on the column to the node on the row, then columns sums
# of matrix K corresponds the the outgoing pathways from each node, and row sums to the ingoing pathways.

# Katz centrality can be normalized by dividing the row sums of the matrix K by the sum of all entries of matrix K

# Here we calculate Katz centrality through a custom function that returns the row sums of the K matrix



katz_bin<-katz(adj.binary)
katz_w<-katz(adj.weight)

if(normalized==TRUE){katz_bin<-katz_bin/sum(katz_bin); katz_w<-katz_w/sum(katz_w)}

#-----------------------#
#-Data Frame of metrics-#
#-----------------------#

results<-data.frame(sample=c(rownames(network), colnames(network)), degree_bin=deg_bin, closeness_bin, katz_bin,
                    degree_w=deg_w, closeness_w, katz_w)

#df_bin<-data.frame(sample=c(rownames(network), colnames(network)), type="binary", degree=deg_bin, closeness=closeness_bin, katz=katz_bin)
#df_w<-data.frame(sample=c(rownames(network), colnames(network)), type="weighted", degree=deg_w, closeness=closeness_w, katz=katz_w)
#results<-rbind(df_bin, df_w)

return(results)

}