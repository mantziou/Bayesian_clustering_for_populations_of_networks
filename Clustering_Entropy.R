#Calculation of Clustering Entropy index for posterior samples of z labels FOR 3 CLUSTERS and 90 NETWORKS IN TOTAL

Clustering_Entr_alt<-function(z_draw,z_true){
  entr<-array(0,dim = c(nrow(z_draw),1))
  for (i in 1:nrow(z_draw)){
    e_1<-e_2<-e_3<-0
    z_draw_c<-z_draw[i, ]
    if(length(which(z_draw_c==1))!=0){ #Calculation of entropy of cluster 1
      p_1<-(1/length(which(z_draw_c==1)))*c(length(which(z_draw_c[which(z_true==1)]==1)),length(which(z_draw_c[which(z_true==2)]==1)),length(which(z_draw_c[which(z_true==3)]==1)))
      #calculates the proportion of times that subject belonging to cluster 1 found in cluster 1, belonging to cluster 2 found in cluster 1, belonging to cluster 3 found in cluster 1
      if(p_1[1]!=0){
        e_1<-e_1+p_1[1]*log2(p_1[1])
      }
      if(p_1[2]!=0){
        e_1<-e_1+p_1[2]*log2(p_1[2])
      }
      if(p_1[3]!=0){
        e_1<-e_1+p_1[3]*log2(p_1[3])
      }
    }
    if(length(which(z_draw_c==2))!=0){ #Calculation of entropy of cluster 2
      p_2<-(1/length(which(z_draw_c==2)))*c(length(which(z_draw_c[which(z_true==1)]==2)),length(which(z_draw_c[which(z_true==2)]==2)),length(which(z_draw_c[which(z_true==3)]==2)))
      #calculates the proportion of times that subject belonging to cluster 1 found in cluster 1, belonging to cluster 2 found in cluster 1, belonging to cluster 3 found in cluster 1
      if(p_2[1]!=0){
        e_2<-e_2+p_2[1]*log2(p_2[1])
      }
      if(p_2[2]!=0){
        e_2<-e_2+p_2[2]*log2(p_2[2])
      }
      if(p_2[3]!=0){
        e_2<-e_2+p_2[3]*log2(p_2[3])
      }
    }
    if(length(which(z_draw_c==3))!=0){ #Calculation of entropy of cluster 3
      p_3<-(1/length(which(z_draw_c==3)))*c(length(which(z_draw_c[which(z_true==1)]==3)),length(which(z_draw_c[which(z_true==2)]==3)),length(which(z_draw_c[which(z_true==3)]==3)))
      #calculates the proportion of times that subject belonging to cluster 1 found in cluster 1, belonging to cluster 2 found in cluster 1, belonging to cluster 3 found in cluster 1
      if(p_3[1]!=0){
        e_3<-e_3+p_3[1]*log2(p_3[1])
      }
      if(p_3[2]!=0){
        e_3<-e_3+p_3[2]*log2(p_3[2])
      }
      if(p_3[3]!=0){
        e_3<-e_3+p_3[3]*log2(p_3[3])
      }
    }
    entr[i]<-(1/length(z_draw_c))*(length(which(z_draw_c==1))*(-e_1)+length(which(z_draw_c==2))*(-e_2)+length(which(z_draw_c==3))*(-e_3))
  }
  return(entr)
}
