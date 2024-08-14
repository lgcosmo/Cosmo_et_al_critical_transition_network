#Function that calculates Katz centrality for a given adjacency matrix

katz<-function(A){

#-----------------------------#
#-Calculating Katz centrality-#
#-----------------------------#

A.eigen<-eigen(A)$values #Calculating eigenvalues
A.maxeigen<-ifelse(is.complex(A.eigen), max(Re(A.eigen[which(Im(A.eigen)==0)])), max(A.eigen)) # Largest eigenvalue

alfa <- 1/(A.maxeigen+0.1) #Calculating alfa

T<-alfa*A #Multiplying A matrix with alfa

I<-diag(1, nrow=nrow(T), ncol=ncol(T)) #Getting the identity matrix of I

katz<-solve(I-T) #Obtaining analytical solution to Katz centrality

katz_cent<-rowSums(katz) #Calculating Katz centrality for each node

return(katz_cent)

}
