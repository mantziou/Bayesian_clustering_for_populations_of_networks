#Calculation of Clustering Purity index for posterior samples of z labels FOR 3 CLUSTERS and 90 NETWORKS IN TOTAL


Clustering_Purity_alt<-function(z_draw,z_true){
  pur<-array(0,dim = c(nrow(z_draw),1))
  for (i in 1:nrow(z_draw)){
    pur_1<-pur_2<-pur_3<-0
    z_draw_c<-z_draw[i, ]
    #purity of cluster 1
    if (length(which(z_draw_c==1))!=0){
      pur_1<-max((1/length(which(z_draw_c==1)))*c(length(which(z_draw_c[which(z_true==1)]==1)),length(which(z_draw_c[which(z_true==2)]==1)),length(which(z_draw_c[which(z_true==3)]==1))))
    }
    #purity of cluster 2
    if (length(which(z_draw_c==2))!=0){
      pur_2<-max((1/length(which(z_draw_c==2)))*c(length(which(z_draw_c[which(z_true==1)]==2)),length(which(z_draw_c[which(z_true==2)]==2)),length(which(z_draw_c[which(z_true==3)]==2))))
    }
    #purity of cluster 3
    if (length(which(z_draw_c==3))!=0){
      pur_3<-max((1/length(which(z_draw_c==3)))*c(length(which(z_draw_c[which(z_true==1)]==3)),length(which(z_draw_c[which(z_true==2)]==3)),length(which(z_draw_c[which(z_true==3)]==3))))
    }
    #total purity
    pur[i]<-(1/length(z_draw_c))*(length(which(z_draw_c==1))*pur_1+length(which(z_draw_c==2))*pur_2+length(which(z_draw_c==3))*pur_3)
  }
  return(pur)
}
