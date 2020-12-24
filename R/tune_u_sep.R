
tune_u_sep = function(m, u_candi, K, X, C=1, iter.max=500, stop=1e-3, trueY=NULL){
  Gm = rep(0,length(u_candi))
  dimen = dim(X)[-length(dim(X))]
  n = dim(X)[length(dim(X))]
  M = length(dim(X))-1
  u = dimen
  #dimen_m = dimenrod(dimen)/dimen[m]
  
  for(i in 1:length(u_candi)) {
    u[m] = u_candi[i] 
    #df[i] = (K-1)*dimenrod(u) + sum((dimen-u)*u + u*(u+1)/2 + (dimen-u)*(dimen-u+1)/2)
    env = TEMM(Xn=X, u=u, K=K, initial="kmeans", iter.max=iter.max, stop=stop, trueY=trueY)
    Gamma = env$Gamma.est
    Mm = env$Mm
    Nm = env$Nm
    
    Gm[i] = log(det(t(Gamma[[m]])%*%Mm[[m]]%*%Gamma[[m]]))+
      log(det(t(Gamma[[m]])%*%solve(Nm[[m]])%*%Gamma[[m]]))
  }
  
  bic = Gm + C*u_candi*log(n)/n
  opt.u = u_candi[which.min(bic)]
  return(list(opt.u=opt.u, bic=bic))
}