### set parameters in simulation###
#' Generate simulation parameters
#'
#' @param model model number
#' @param seed random seed
#'
#' @return List of parameters
#' @export
#' @import stats
#'
#' @examples
#' para_gen(model=1)
para_gen = function(model,seed=999){
  set.seed(seed)

  if(model==1){
    dimen = c(20,20)
    M = length(dimen)
    p = prod(dimen)
    u = c(2,2)
    pi0 = c(0.5,0.5)
    K = length(pi0)

    Omega = Omega0 = list()
    Gamma = Gamma0 = list()
    SIG = SIGsqrt = list()

    for (i in 1:M) {
      Gamma[[i]] = matrix(runif(dimen[i]*u[i]), dimen[i], u[i])
      Gamma[[i]] = qr.Q(qr(Gamma[[i]]))
      O = matrix(runif(u[i]^2), u[i], u[i])
      O = qr.Q(qr(O))
      D = 5*diag(1:u[i])
      Omega[[i]] = O %*% D %*% t(O)

      if(u[i]==dimen[i]){
        Gamma0[[i]] = matrix(0,dimen[i],dimen[i])
        Omega0[[i]] = matrix(0,dimen[i],dimen[i])
      }
      else{
        Gamma0[[i]] = qr.Q(qr(Gamma[[i]]), complete=TRUE)[ ,(u[i]+1):dimen[i]]
        Op = matrix(runif((dimen[i]-u[i])^2), (dimen[i]-u[i]), (dimen[i]-u[i]))
        Op = qr.Q(qr(Op))
        D0 = diag(exp(seq(i, -10, length.out=dimen[i]-u[i])), nrow=dimen[i]-u[i])
        Omega0[[i]] = Op %*% D0 %*% t(Op)
      }

      SIG[[i]] = Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
        Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
      SIG[[i]] = 0.55*SIG[[i]]/norm(SIG[[i]], type="F")
      SIGsqrt[[i]] = pracma::sqrtm(SIG[[i]])$B
    }
    SIGinv = lapply(SIG, FUN=MASS::ginv)

    Alpha_temp = array(runif(prod(u)*K), c(u,K))
    Alpha = asplit(Alpha_temp,M+1)
    Muk = list()  #cluster mean, a list of length K
    for(k in 1:K) {
      Muk[[k]] = rTensor::ttl(rTensor::as.tensor(Alpha[[k]]), Gamma, ms=c(1:M))@data
    }
  }


  if(model==2.1){
    dimen = c(50,50)
    M = length(dimen)
    p = prod(dimen)
    u = c(5,5)
    pi0 = rep(0.5,2)
    K = length(pi0)

    Omega = Omega0 = list()
    Gamma = Gamma0 = list()
    SIG = SIGsqrt = list()

    for (i in 1:M) {
      Gamma[[i]] = matrix(runif(dimen[i]*u[i]), dimen[i], u[i])
      Gamma[[i]] = qr.Q(qr(Gamma[[i]]))
      O = matrix(runif(u[i]^2), u[i], u[i])
      O = qr.Q(qr(O))
      D = 5*diag(1:u[i])
      Omega[[i]] = O %*% D %*% t(O)

      if(u[i]==dimen[i]){
        Gamma0[[i]] = matrix(0,dimen[i],dimen[i])
        Omega0[[i]] = matrix(0,dimen[i],dimen[i])
      }
      else{
        Gamma0[[i]] = qr.Q(qr(Gamma[[i]]), complete=TRUE)[ ,(u[i]+1):dimen[i]]
        Op = matrix(runif((dimen[i]-u[i])^2), (dimen[i]-u[i]), (dimen[i]-u[i]))
        Op = qr.Q(qr(Op))
        D0 = diag(exp(seq(i, -10, length.out=dimen[i]-u[i])), nrow=dimen[i]-u[i])
        Omega0[[i]] = Op %*% D0 %*% t(Op)
      }

      SIG[[i]] = Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
        Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
      SIG[[i]] = 2*SIG[[i]]/norm(SIG[[i]], type="F")
      SIGsqrt[[i]] = pracma::sqrtm(SIG[[i]])$B
    }
    SIGinv = lapply(SIG, FUN=MASS::ginv)

    Alpha_temp = array(runif(prod(u)*K), c(u,K))
    Alpha = asplit(Alpha_temp,M+1)
    Muk = list()  #cluster mean, a list of length K
    for(k in 1:K) {
      Muk[[k]] = rTensor::ttl(rTensor::as.tensor(Alpha[[k]]), Gamma, ms=c(1:M))@data
    }
  }


  # based on model 2.1, change number of clusters to 3
  if(model==2.2){
    dimen = c(50,50)
    M = length(dimen)
    p = prod(dimen)
    u = c(5,5)
    pi0 = rep(1/3,3)
    K = length(pi0)

    Omega = Omega0 = list()
    Gamma = Gamma0 = list()
    SIG = SIGsqrt = list()

    for (i in 1:M) {
      Gamma[[i]] = matrix(runif(dimen[i]*u[i]), dimen[i], u[i])
      Gamma[[i]] = qr.Q(qr(Gamma[[i]]))
      O = matrix(runif(u[i]^2), u[i], u[i])
      O = qr.Q(qr(O))
      D = 5*diag(1:u[i])
      Omega[[i]] = O %*% D %*% t(O)

      if(u[i]==dimen[i]){
        Gamma0[[i]] = matrix(0,dimen[i],dimen[i])
        Omega0[[i]] = matrix(0,dimen[i],dimen[i])
      }
      else{
        Gamma0[[i]] = qr.Q(qr(Gamma[[i]]), complete=TRUE)[ ,(u[i]+1):dimen[i]]
        Op = matrix(runif((dimen[i]-u[i])^2), (dimen[i]-u[i]), (dimen[i]-u[i]))
        Op = qr.Q(qr(Op))
        D0 = diag(exp(seq(i, -10, length.out=dimen[i]-u[i])), nrow=dimen[i]-u[i])
        Omega0[[i]] = Op %*% D0 %*% t(Op)
      }

      SIG[[i]] = Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
        Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
      SIG[[i]] = 2*SIG[[i]]/norm(SIG[[i]], type="F")
      SIGsqrt[[i]] = pracma::sqrtm(SIG[[i]])$B
    }
    SIGinv = lapply(SIG, FUN=MASS::ginv)

    Alpha_temp = array(runif(prod(u)*K), c(u,K))
    Alpha = asplit(Alpha_temp,M+1)
    Muk = list()  #cluster mean, a list of length K
    for(k in 1:K) {
      Muk[[k]] = rTensor::ttl(rTensor::as.tensor(Alpha[[k]]), Gamma, ms=c(1:M))@data
    }
  }


  # based on model 2.1, change number of clusters to 4
  if(model==2.3){
    dimen = c(50,50)
    M = length(dimen)
    p = prod(dimen)
    u = c(5,5)
    pi0 = rep(0.25,4)
    K = length(pi0)

    Omega = Omega0 = list()
    Gamma = Gamma0 = list()
    SIG = SIGsqrt = list()

    for (i in 1:M) {
      Gamma[[i]] = matrix(runif(dimen[i]*u[i]), dimen[i], u[i])
      Gamma[[i]] = qr.Q(qr(Gamma[[i]]))
      O = matrix(runif(u[i]^2), u[i], u[i])
      O = qr.Q(qr(O))
      D = 5*diag(1:u[i])
      Omega[[i]] = O %*% D %*% t(O)

      if(u[i]==dimen[i]){
        Gamma0[[i]] = matrix(0,dimen[i],dimen[i])
        Omega0[[i]] = matrix(0,dimen[i],dimen[i])
      }
      else{
        Gamma0[[i]] = qr.Q(qr(Gamma[[i]]), complete=TRUE)[ ,(u[i]+1):dimen[i]]
        Op = matrix(runif((dimen[i]-u[i])^2), (dimen[i]-u[i]), (dimen[i]-u[i]))
        Op = qr.Q(qr(Op))
        D0 = diag(exp(seq(i, -10, length.out=dimen[i]-u[i])), nrow=dimen[i]-u[i])
        Omega0[[i]] = Op %*% D0 %*% t(Op)
      }

      SIG[[i]] = Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
        Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
      SIG[[i]] = 2*SIG[[i]]/norm(SIG[[i]], type="F")
      SIGsqrt[[i]] = pracma::sqrtm(SIG[[i]])$B
    }
    SIGinv = lapply(SIG, FUN=MASS::ginv)

    Alpha_temp = array(runif(prod(u)*K), c(u,K))
    Alpha = asplit(Alpha_temp,M+1)
    Muk = list()  #cluster mean, a list of length K
    for(k in 1:K) {
      Muk[[k]] = rTensor::ttl(rTensor::as.tensor(Alpha[[k]]), Gamma, ms=c(1:M))@data
    }
  }


  if(model==3.1){
    dimen = c(10,10,10)
    M = length(dimen)
    p = prod(dimen)
    u = c(1,2,3)
    pi0 = rep(0.5,2)
    K = length(pi0)

    Omega = Omega0 = list()
    Gamma = Gamma0 = list()
    SIG = SIGsqrt = list()

    for (i in 1:M) {
      Gamma[[i]] = matrix(runif(dimen[i]*u[i]), dimen[i], u[i])
      Gamma[[i]] = qr.Q(qr(Gamma[[i]]))
      O = matrix(runif(u[i]^2), u[i], u[i])
      O = qr.Q(qr(O))
      D = 5*diag(1:u[i])
      Omega[[i]] = O %*% D %*% t(O)

      if(u[i]==dimen[i]){
        Gamma0[[i]] = matrix(0,dimen[i],dimen[i])
        Omega0[[i]] = matrix(0,dimen[i],dimen[i])
      }
      else{
        Gamma0[[i]] = qr.Q(qr(Gamma[[i]]), complete=TRUE)[ ,(u[i]+1):dimen[i]]
        Op = matrix(runif((dimen[i]-u[i])^2), (dimen[i]-u[i]), (dimen[i]-u[i]))
        Op = qr.Q(qr(Op))
        D0 = diag(exp(seq(i, -10, length.out=dimen[i]-u[i])), nrow=dimen[i]-u[i])
        Omega0[[i]] = Op %*% D0 %*% t(Op)
      }

      SIG[[i]] = Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
        Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
      SIG[[i]] = 1.1*SIG[[i]]/norm(SIG[[i]], type="F")
      SIGsqrt[[i]] = pracma::sqrtm(SIG[[i]])$B
    }
    SIGinv = lapply(SIG, FUN=MASS::ginv)

    Alpha_temp = array(runif(prod(u)*K), c(u,K))
    Alpha = asplit(Alpha_temp,M+1)
    Muk = list()  #cluster mean, a list of length K
    for(k in 1:K) {
      Muk[[k]] = rTensor::ttl(rTensor::as.tensor(Alpha[[k]]), Gamma, ms=c(1:M))@data
    }
  }


  # based on model 3.1, change the structure of covariance matrices
  if(model==3.2){
    dimen = c(10,10,10)
    M = length(dimen)
    p = prod(dimen)
    u = c(1,2,3)
    pi0 = rep(0.5,2)
    K = length(pi0)

    Omega = Omega0 = list()
    Gamma = Gamma0 = list()
    SIG = SIGsqrt = list()

    for (i in 1:M) {
      Gamma[[i]] = matrix(runif(dimen[i]*u[i]), dimen[i], u[i])
      Gamma[[i]] = qr.Q(qr(Gamma[[i]]))
      O = matrix(runif(u[i]^2), u[i], u[i])
      O = qr.Q(qr(O))
      D = 5*diag(1:u[i])
      Omega[[i]] = O %*% D %*% t(O)

      if(u[i]==dimen[i]){
        Gamma0[[i]] = matrix(0,dimen[i],dimen[i])
        Omega0[[i]] = matrix(0,dimen[i],dimen[i])
      }
      else{
        Gamma0[[i]] = qr.Q(qr(Gamma[[i]]), complete=TRUE)[ ,(u[i]+1):dimen[i]]
        Op = matrix(runif((dimen[i]-u[i])^2), (dimen[i]-u[i]), (dimen[i]-u[i]))
        Op = qr.Q(qr(Op))
        D0 = diag(exp(seq(i, -10, length.out=dimen[i]-u[i])), nrow=dimen[i]-u[i])
        Omega0[[i]] = Op %*% D0 %*% t(Op)
      }

      SIG[[i]] = Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
        Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
      SIG[[i]] = SIG[[i]] + matrix(0.2,dimen[i],dimen[i]) + diag(1-0.2,dimen[i],dimen[i])
      SIG[[i]] = 1.1*SIG[[i]]/norm(SIG[[i]], type="F")
      SIGsqrt[[i]] = pracma::sqrtm(SIG[[i]])$B
    }
    SIGinv = lapply(SIG, FUN=MASS::ginv)

    Alpha_temp = array(runif(prod(u)*K), c(u,K))
    Alpha = asplit(Alpha_temp,M+1)
    Muk = list()  #cluster mean, a list of length K
    for(k in 1:K) {
      Muk[[k]] = rTensor::ttl(rTensor::as.tensor(Alpha[[k]]), Gamma, ms=c(1:M))@data
    }
  }


  # model 4a, based on model 3.1, add a fourth tensor mode
  if(model==4.1){
    dimen = c(10,10,10,2)
    M = length(dimen)
    p = prod(dimen)
    u = c(1,2,3,2)
    pi0 = rep(0.5,2)
    K = length(pi0)

    Omega = Omega0 = list()
    Gamma = Gamma0 = list()
    SIG = SIGsqrt = list()

    for (i in 1:M) {
      Gamma[[i]] = matrix(runif(dimen[i]*u[i]), dimen[i], u[i])
      Gamma[[i]] = qr.Q(qr(Gamma[[i]]))
      O = matrix(runif(u[i]^2), u[i], u[i])
      O = qr.Q(qr(O))
      D = 5*diag(1:u[i])
      Omega[[i]] = O %*% D %*% t(O)

      if(u[i]==dimen[i]){
        Gamma0[[i]] = matrix(0,dimen[i],dimen[i])
        Omega0[[i]] = matrix(0,dimen[i],dimen[i])
      }
      else{
        Gamma0[[i]] = qr.Q(qr(Gamma[[i]]), complete=TRUE)[ ,(u[i]+1):dimen[i]]
        Op = matrix(runif((dimen[i]-u[i])^2), (dimen[i]-u[i]), (dimen[i]-u[i]))
        Op = qr.Q(qr(Op))
        D0 = diag(exp(seq(i, -10, length.out=dimen[i]-u[i])), nrow=dimen[i]-u[i])
        Omega0[[i]] = Op %*% D0 %*% t(Op)
      }

      SIG[[i]] = Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
        Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
      if(i==4){
        SIG[[i]] = diag(dimen[4])
      }
      SIG[[i]] = 1.2*SIG[[i]]/norm(SIG[[i]], type="F")
      SIGsqrt[[i]] = pracma::sqrtm(SIG[[i]])$B
    }
    SIGinv = lapply(SIG, FUN=MASS::ginv)

    Alpha_temp = array(runif(prod(u)*K), c(u,K))
    Alpha = asplit(Alpha_temp,M+1)
    Muk = list()  #cluster mean, a list of length K
    for(k in 1:K) {
      Muk[[k]] = rTensor::ttl(rTensor::as.tensor(Alpha[[k]]), Gamma, ms=c(1:M))@data
    }
  }


  # model 4b
  if(model==4.2){
    dimen = c(10,10,10,2)
    M = length(dimen)
    p = prod(dimen)
    u = c(1,2,3,2)
    pi0 = rep(0.5,2)
    K = length(pi0)

    Omega = Omega0 = list()
    Gamma = Gamma0 = list()
    SIG = SIGsqrt = list()

    for (i in 1:M) {
      Gamma[[i]] = matrix(runif(dimen[i]*u[i]), dimen[i], u[i])
      Gamma[[i]] = qr.Q(qr(Gamma[[i]]))
      O = matrix(runif(u[i]^2), u[i], u[i])
      O = qr.Q(qr(O))
      D = 5*diag(1:u[i])
      Omega[[i]] = O %*% D %*% t(O)

      if(u[i]==dimen[i]){
        Gamma0[[i]] = matrix(0,dimen[i],dimen[i])
        Omega0[[i]] = matrix(0,dimen[i],dimen[i])
      }
      else{
        Gamma0[[i]] = qr.Q(qr(Gamma[[i]]), complete=TRUE)[ ,(u[i]+1):dimen[i]]
        Op = matrix(runif((dimen[i]-u[i])^2), (dimen[i]-u[i]), (dimen[i]-u[i]))
        Op = qr.Q(qr(Op))
        D0 = diag(exp(seq(i, -10, length.out=dimen[i]-u[i])), nrow=dimen[i]-u[i])
        Omega0[[i]] = Op %*% D0 %*% t(Op)
      }

      SIG[[i]] = Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
        Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
      if(i==4){
        SIG[[i]] = diag(dimen[4])
      }
      SIG[[i]] = 2*SIG[[i]]/norm(SIG[[i]], type="F")
      SIGsqrt[[i]] = pracma::sqrtm(SIG[[i]])$B
    }
    SIGinv = lapply(SIG, FUN=MASS::ginv)

    Alpha_temp = array(runif(prod(u)*K), c(u,K))
    Alpha = asplit(Alpha_temp,M+1)
    Muk = list()  #cluster mean, a list of length K
    for(k in 1:K) {
      Muk[[k]] = rTensor::ttl(rTensor::as.tensor(Alpha[[k]]), Gamma, ms=c(1:M))@data
    }
  }


  # similar to model 4.2, change fourth dimension to 4
  if(model==4.3){
    dimen = c(10,10,10,4)
    M = length(dimen)
    p = prod(dimen)
    u = c(1,2,3,4)
    pi0 = rep(0.5,2)
    K = length(pi0)

    Omega = Omega0 = list()
    Gamma = Gamma0 = list()
    SIG = SIGsqrt = list()

    for (i in 1:M) {
      Gamma[[i]] = matrix(runif(dimen[i]*u[i]), dimen[i], u[i])
      Gamma[[i]] = qr.Q(qr(Gamma[[i]]))
      O = matrix(runif(u[i]^2), u[i], u[i])
      O = qr.Q(qr(O))
      D = 5*diag(1:u[i])
      Omega[[i]] = O %*% D %*% t(O)

      if(u[i]==dimen[i]){
        Gamma0[[i]] = matrix(0,dimen[i],dimen[i])
        Omega0[[i]] = matrix(0,dimen[i],dimen[i])
      }
      else{
        Gamma0[[i]] = qr.Q(qr(Gamma[[i]]), complete=TRUE)[ ,(u[i]+1):dimen[i]]
        Op = matrix(runif((dimen[i]-u[i])^2), (dimen[i]-u[i]), (dimen[i]-u[i]))
        Op = qr.Q(qr(Op))
        D0 = diag(exp(seq(i, -10, length.out=dimen[i]-u[i])), nrow=dimen[i]-u[i])
        Omega0[[i]] = Op %*% D0 %*% t(Op)
      }

      SIG[[i]] = Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
        Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
      if(i==4){
        SIG[[i]] = diag(dimen[4])
      }
      SIG[[i]] = 2*SIG[[i]]/norm(SIG[[i]], type="F")
      SIGsqrt[[i]] = pracma::sqrtm(SIG[[i]])$B
    }
    SIGinv = lapply(SIG, FUN=MASS::ginv)

    Alpha_temp = array(runif(prod(u)*K), c(u,K))
    Alpha = asplit(Alpha_temp,M+1)
    Muk = list()  #cluster mean, a list of length K
    for(k in 1:K) {
      Muk[[k]] = rTensor::ttl(rTensor::as.tensor(Alpha[[k]]), Gamma, ms=c(1:M))@data
    }
  }


  # similar to model 4.2, change fourth dimension to 6
  if(model==4.4){
    dimen = c(10,10,10,6)
    M = length(dimen)
    p = prod(dimen)
    u = c(1,2,3,6)
    pi0 = rep(0.5,2)
    K = length(pi0)

    Omega = Omega0 = list()
    Gamma = Gamma0 = list()
    SIG = SIGsqrt = list()

    for (i in 1:M) {
      Gamma[[i]] = matrix(runif(dimen[i]*u[i]), dimen[i], u[i])
      Gamma[[i]] = qr.Q(qr(Gamma[[i]]))
      O = matrix(runif(u[i]^2), u[i], u[i])
      O = qr.Q(qr(O))
      D = 5*diag(1:u[i])
      Omega[[i]] = O %*% D %*% t(O)

      if(u[i]==dimen[i]){
        Gamma0[[i]] = matrix(0,dimen[i],dimen[i])
        Omega0[[i]] = matrix(0,dimen[i],dimen[i])
      }
      else{
        Gamma0[[i]] = qr.Q(qr(Gamma[[i]]), complete=TRUE)[ ,(u[i]+1):dimen[i]]
        Op = matrix(runif((dimen[i]-u[i])^2), (dimen[i]-u[i]), (dimen[i]-u[i]))
        Op = qr.Q(qr(Op))
        D0 = diag(exp(seq(i, -10, length.out=dimen[i]-u[i])), nrow=dimen[i]-u[i])
        Omega0[[i]] = Op %*% D0 %*% t(Op)
      }

      SIG[[i]] = Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
        Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
      if(i==4){
        SIG[[i]] = diag(dimen[4])
      }
      SIG[[i]] = 2*SIG[[i]]/norm(SIG[[i]], type="F")
      SIGsqrt[[i]] = pracma::sqrtm(SIG[[i]])$B
    }
    SIGinv = lapply(SIG, FUN=MASS::ginv)

    Alpha_temp = array(runif(prod(u)*K), c(u,K))
    Alpha = asplit(Alpha_temp,M+1)
    Muk = list()  #cluster mean, a list of length K
    for(k in 1:K) {
      Muk[[k]] = rTensor::ttl(rTensor::as.tensor(Alpha[[k]]), Gamma, ms=c(1:M))@data
    }
  }


  # similar to model 4.2, change fourth dimension to 8
  if(model==4.5){
    dimen = c(10,10,10,8)
    M = length(dimen)
    p = prod(dimen)
    u = c(1,2,3,8)
    pi0 = rep(0.5,2)
    K = length(pi0)

    Omega = Omega0 = list()
    Gamma = Gamma0 = list()
    SIG = SIGsqrt = list()

    for (i in 1:M) {
      Gamma[[i]] = matrix(runif(dimen[i]*u[i]), dimen[i], u[i])
      Gamma[[i]] = qr.Q(qr(Gamma[[i]]))
      O = matrix(runif(u[i]^2), u[i], u[i])
      O = qr.Q(qr(O))
      D = 5*diag(1:u[i])
      Omega[[i]] = O %*% D %*% t(O)

      if(u[i]==dimen[i]){
        Gamma0[[i]] = matrix(0,dimen[i],dimen[i])
        Omega0[[i]] = matrix(0,dimen[i],dimen[i])
      }
      else{
        Gamma0[[i]] = qr.Q(qr(Gamma[[i]]), complete=TRUE)[ ,(u[i]+1):dimen[i]]
        Op = matrix(runif((dimen[i]-u[i])^2), (dimen[i]-u[i]), (dimen[i]-u[i]))
        Op = qr.Q(qr(Op))
        D0 = diag(exp(seq(i, -10, length.out=dimen[i]-u[i])), nrow=dimen[i]-u[i])
        Omega0[[i]] = Op %*% D0 %*% t(Op)
      }

      SIG[[i]] = Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
        Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
      if(i==4){
        SIG[[i]] = diag(dimen[4])
      }
      SIG[[i]] = 2*SIG[[i]]/norm(SIG[[i]], type="F")
      SIGsqrt[[i]] = pracma::sqrtm(SIG[[i]])$B
    }
    SIGinv = lapply(SIG, FUN=MASS::ginv)

    Alpha_temp = array(runif(prod(u)*K), c(u,K))
    Alpha = asplit(Alpha_temp,M+1)
    Muk = list()  #cluster mean, a list of length K
    for(k in 1:K) {
      Muk[[k]] = rTensor::ttl(rTensor::as.tensor(Alpha[[k]]), Gamma, ms=c(1:M))@data
    }
  }


  # similar to model 4.2, change fourth dimension to 10
  if(model==4.6){
    dimen = c(10,10,10,10)
    M = length(dimen)
    p = prod(dimen)
    u = c(1,2,3,10)
    pi0 = rep(0.5,2)
    K = length(pi0)

    Omega = Omega0 = list()
    Gamma = Gamma0 = list()
    SIG = SIGsqrt = list()

    for (i in 1:M) {
      Gamma[[i]] = matrix(runif(dimen[i]*u[i]), dimen[i], u[i])
      Gamma[[i]] = qr.Q(qr(Gamma[[i]]))
      O = matrix(runif(u[i]^2), u[i], u[i])
      O = qr.Q(qr(O))
      D = 5*diag(1:u[i])
      Omega[[i]] = O %*% D %*% t(O)

      if(u[i]==dimen[i]){
        Gamma0[[i]] = matrix(0,dimen[i],dimen[i])
        Omega0[[i]] = matrix(0,dimen[i],dimen[i])
      }
      else{
        Gamma0[[i]] = qr.Q(qr(Gamma[[i]]), complete=TRUE)[ ,(u[i]+1):dimen[i]]
        Op = matrix(runif((dimen[i]-u[i])^2), (dimen[i]-u[i]), (dimen[i]-u[i]))
        Op = qr.Q(qr(Op))
        D0 = diag(exp(seq(i, -10, length.out=dimen[i]-u[i])), nrow=dimen[i]-u[i])
        Omega0[[i]] = Op %*% D0 %*% t(Op)
      }

      SIG[[i]] = Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
        Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
      if(i==4){
        SIG[[i]] = diag(dimen[4])
      }
      SIG[[i]] = 2*SIG[[i]]/norm(SIG[[i]], type="F")
      SIGsqrt[[i]] = pracma::sqrtm(SIG[[i]])$B
    }
    SIGinv = lapply(SIG, FUN=MASS::ginv)

    Alpha_temp = array(runif(prod(u)*K), c(u,K))
    Alpha = asplit(Alpha_temp,M+1)
    Muk = list()  #cluster mean, a list of length K
    for(k in 1:K) {
      Muk[[k]] = rTensor::ttl(rTensor::as.tensor(Alpha[[k]]), Gamma, ms=c(1:M))@data
    }
  }


  return(list(dimen=dimen, M=M, p=p, u=u, pi0=pi0, K=K, Muk=Muk,
              Gamma=Gamma, SIG=SIG, SIGsqrt=SIGsqrt, SIGinv=SIGinv))
}













