univalboot<-function(y,FP,relip){

  k<-size(FP)[2]

  dife<-matrix(0,1,k)

  for (i in 1:k){
    tmp<-cor(FP[,i],y)
    tmp<-tmp/sqrt(relip[i])
    dife[i]<-tmp
    remove(tmp)
  }
  dife=t(dife)

}
