univalboot2<-function(y,FP,fg,relip){

  k<-size(FP)[2]

  dife<-matrix(0,1,k)
  difeg<-matrix(0,1,k)

  for (i in 1:k){

    tmp1<-cor(FP[,i],y)
    tmp1<-tmp1/sqrt(relip[i])
    dife[i]<-tmp1
    remove(tmp1)
    tmp2<-cor(FP[,i],fg)
    tmp2<-tmp2/sqrt(relip[i])
    difeg[i]<-tmp2
    remove(tmp2)

  }
  dife=t(dife)
  difeg=t(difeg)

  OUT<-list('dife'=dife,'difeg'=difeg)
  return(OUT)

}
