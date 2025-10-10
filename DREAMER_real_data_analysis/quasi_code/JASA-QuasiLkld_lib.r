library(dplyr)
library(MASS)
library(glmnet)
library(scalreg)
library(tictoc)


Psi.def<- function(eta,q, cov='diag'){ ##generating different Psi
  if (cov=='diag') {psi = diag(eta,q)}
  if (cov=='sym') {psi=toeplitz(eta^(0:(q-1))); d=1}
  if (cov=='semi') {psi = diag(c(rep(eta[1],q/2),rep(0,q/2)));d=q}
  psi
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


CV_quasi_fixed_effects = function(X, y, Z, grp, a.seq=seq(1, 10, 0.5) ) {
  if (!is.element("matrix", class(X) )) {
    X = matrix(X, ncol=1);
  }

  if (!is.element("matrix", class(Z) )) {
    Z = matrix(Z, ncol=1);
  }

  stopifnot(dim(X)[1] == length(y) );
  stopifnot(dim(Z)[1] == length(y) );
  stopifnot(length(grp) == length(y) );

  order_ixs = order(grp);
  X = matrix(X[order_ixs,], ncol=ncol(X) );
  y = y[order_ixs];
  Z = matrix(Z[order_ixs,], ncol=ncol(Z) );
  grp = grp[order_ixs];

  grp_uniq = unique(grp);

  N = length(y);
  n = length(grp_uniq);

  group_start_stop_ixs = (
    data.frame(id=grp) %>%
    mutate(ix=1:N) %>%
    group_by(id) %>%
    summarize(ax=min(ix), bx=max(ix), .groups="drop") );

  for (gx in 1:n) {
    ax = group_start_stop_ixs$ax[gx];
    bx = group_start_stop_ixs$bx[gx];
    stopifnot(all(group_start_stop_ixs$id[gx] == grp[ax:bx]) );
  }

  group_start_stop_ixs = t(as.matrix(group_start_stop_ixs[,c("ax", "bx")]) );

  n_folds = 10;
  grps_random_order = sample(1:n);
  n_leftover = n %% n_folds;

  if (n_leftover > 0) {
    grps_random_order = c(grps_random_order, rep(NA, n_folds - n_leftover) );
  }

  test_group_ixs = t(matrix(grps_random_order, n_folds) );
  test_mses = matrix(NA, length(a.seq), n_folds);

  for (tx in 1:n_folds) {
    cat(sprintf("Beginning fold %d ...\n", tx) );

    test_group_tx_ixs = sort(as.numeric(na.omit(test_group_ixs[,tx]) ));
    test_ixs = as.vector(
      unlist(
        sapply(
          test_group_tx_ixs,
          \(gx) group_start_stop_ixs[1,gx]:group_start_stop_ixs[2,gx]) ));

    train_ixs = setdiff(1:N, test_ixs);

    X_train = matrix(X[train_ixs,], ncol=ncol(X) );
    Z_train = matrix(Z[train_ixs,], ncol=ncol(Z) );
    y_train = y[train_ixs];
    grp_train = grp[train_ixs];

    X_test = matrix(X[test_ixs,], ncol=ncol(X) );
    y_test = y[test_ixs];

    for (lx in 1:length(a.seq) ) {
      quasi_mod = Fixed_effects_estimation(
        X_train,
        y_train,
        Z_train,
        grp_train,
        a.seq[lx]);

      test_mses[lx,tx] = mean((y_test - (X_test %*% quasi_mod$beta.hat) )^2);
    }
  }

  test_mses = cbind(a.seq, as.data.frame(test_mses) );
  colnames(test_mses) = c("a", paste0("fold_", 1:n_folds) );

  test_mses
};


compute_blup_quasi = function(resids, Z, sigma_Z, grp) {
  grp_uniq = unique(grp);

  n_dim = length(grp_uniq);
  m_dim = length(grp);

  stopifnot(nrow(Z) == m_dim);

  Q_dim = ncol(Z);

  blup = matrix(NA, Q_dim, n_dim);

  for (gx in 1:n_dim) {
    ixs_grp_gx = which(grp == grp_uniq[gx]);
    m_dim_gx = length(ixs_grp_gx);
    resids_gx = resids[ixs_grp_gx];
    Z_gx = matrix(Z[ixs_grp_gx,], ncol=ncol(Z) );
    sigma_Y_gx = tcrossprod(Z_gx %*% sigma_Z, Z_gx) + diag(1, m_dim_gx);
    blup[,gx] = sigma_Z %*% crossprod(Z_gx, solve(sigma_Y_gx) ) %*% resids_gx;
  }

  list(group_ids=grp_uniq, blup=blup)
};


predict_quasi = function(X_test, Z_test, b_hat, blup_map, grp_test) {
  all_factor_levels = unique(
    c(levels(blup_map$group_ids), levels(grp_test) ));

  grp_test = factor(grp_test, levels=all_factor_levels);
  blup_map$group_ids = factor(blup_map$group_ids, levels=all_factor_levels);

  grp_uniq = unique(grp_test);

  n_dim = length(grp_uniq);
  m_dim = length(grp_test);

  stopifnot(nrow(X_test) == m_dim);
  stopifnot(nrow(Z_test) == m_dim);

  Q_dim = ncol(Z_test);

  y_hat = X_test %*% b_hat;

  for (gx in 1:n_dim) {
    if (is.element(grp_uniq[gx], blup_map$group_ids) ) {
      ixs_grp_gx = which(grp_test == grp_uniq[gx]);
      ix_blup = which(blup_map$group_ids == grp_uniq[gx]);
      Z_gx = matrix(Z_test[ixs_grp_gx,], ncol=ncol(Z_test) );
      y_hat[ixs_grp_gx] = y_hat[ixs_grp_gx] + (Z_gx %*% blup_map$blup[ix_blup]);
    }
  }

  as.numeric(y_hat)
};


#main function for fixed effect estimation & inference
#Inputs-- X: design for fixed effects; y: response; z: design for random effects;
##grp: a vector of length N indicating group membership; a : tuning parameter in $\Sig_a$
##lm: linear model fitting; inf.coord: the coordinates of beta for inference.
#Outputs-- beta.hat: the Lasso estimator based on psedo-likelihood; 
##beta.db:debiased Lasso estimators for the fixed effects in inf.coord;
##beta.db.sd: standard deviation of the debiased Lasso estimators for the fixed effects in inf.coord.
Fixed_effects_estimation <- function(X, y, z, grp, a=8, lm=F, inf.coord=NULL){
  grp_uniq = unique(grp);

  N=length(y)
  n=length(grp_uniq);
  q=ncol(z)
  p=ncol(X)
  ###preprocessing
  #X <- scale (X)
  #X.sd <- attributes(X)$`scaled:scale`
  #y <- y - mean(y)

  if(lm){#linear model fitting
    sig.init<-scalreg(X,y)$hsigma
    beta.hat <- as.numeric(glmnet(X, y, lambda = sig.init*sqrt(2*log(p)/N))$beta)
    return(list(beta.hat=beta.hat, Sig.a.inv.half=NULL))
  }else{#mixed effect model fitting
    X.a<-X
    y.a<-y
    tr.a<-0
    for (i in 1:n){
      cur.mem=which(grp==grp_uniq[i])
      mi=length(cur.mem)
      zi = as.matrix(z[cur.mem,])
      sigmai = a*tcrossprod(zi)+diag(rep(1,mi))
      sigmai.svd=svd(sigmai)
      Sig.a.inv.half <- sigmai.svd$u %*% diag(1/sqrt(sigmai.svd$d)) %*% t(sigmai.svd$u)
      X.a[cur.mem,] =  Sig.a.inv.half %*% as.matrix(X[cur.mem,]) 
      y.a[cur.mem] =  Sig.a.inv.half %*% as.matrix(y[cur.mem])
      tr.a<-tr.a+ sum(1/sigmai.svd$d) #compute the effective sample size
    }
    sig.init<-scalreg(X.a,y.a)$hsigma #scaled-lasso to get tuning parameter
    beta.hat = glmnet(X.a, y.a, lambda = sig.init*sqrt(2*log(p)/N))$beta
    #cv.init<-cv.glmnet(X.a, y.a, lam=seq(1.2,0.2,-0.1)*sqrt(2*log(p)/N))
    #beta.hat<-coef(cv.init, s=cv.init$lambda.min)[-1]

    
    ###debiased Lasso for mlm
    beta.db.sd.mlm<-rep(NA,length=length(inf.coord))
    beta.db.mlm<-rep(NA,length=length(inf.coord))
    if(is.null(inf.coord)){
      return( list(beta.hat = beta.hat, beta.db=beta.db.mlm, beta.db.sd=beta.db.sd.mlm, tr.a=tr.a))
    }
    res<-y.a-X.a %*% beta.hat
    lam.seq<-seq(1.2,0.1,-0.1)*sqrt(2*log(p)/N)
    for(j in 1:length(inf.coord)){
      col.j<-inf.coord[j]
      #cv.x<-cv.glmnet(X.a[,-col.j], X.a[,col.j], lambda=lam.seq)
      #gam.j<- coef(cv.x, s=cv.x$lambda.min)[-1]
      sig.x<-scalreg(X.a[,-col.j], X.a[,col.j])$hsigma
      kappa.hat.mlm <- glmnet(X.a[,-col.j], X.a[,col.j], lambda=sig.x*sqrt(2*log(p)/N))
      gam.j <- as.numeric(kappa.hat.mlm$beta)
      wj.mlm <- X.a[,col.j] - X.a[,-col.j]%*%gam.j
      beta.db.mlm[j] = beta.hat[col.j] + sum( wj.mlm * res)/sum(wj.mlm * X.a[,col.j])
      num=0
      for(i in 1:n){
        cur.mem=((i-1)*m+1):(i*m)
        num <-num+sum(wj.mlm[cur.mem]*res[cur.mem])^2
      }
      beta.db.sd.mlm[j] = sqrt(num)/sum(wj.mlm*X.a[,col.j])
    }
    return( list(beta.hat = beta.hat, beta.db=beta.db.mlm, beta.db.sd=beta.db.sd.mlm, tr.a=tr.a))
  }
}

##################variance components estimation ####################
#main function for variance components estimation
#Inputs-- r.hat: empirical residuals computed from y-X%*%beta.hat; z: design for random effects;
##G.list: a list of basis functions for Psi_eta; grp: a vector of length N indicating group membership; a : tuning parameter in $\Sig_a$
#Outputs-- eta.hat: estimator of eta; sig.e.sq.hat: estimator of sig_e^2;
##is.spd: is Psi_{eta.hat} positive definite?
Varcomp.est<- function(r.hat,z, G.list, a, grp){
  grp_uniq = unique(grp);

  q=ncol(z)
  n = length(grp_uniq)
  ###estimate sig^2.e##
  denom<-0
  num<-0
  for (i in 1:n){
    cur.mem=which(grp==grp_uniq[i])
    if(length(cur.mem) <= q) next
    denom = denom + length(cur.mem)-q
    num <- num + sum((lm(r.hat[cur.mem]~z[cur.mem,]-1)$res)^2)
  }
  sig.e.sq.hat<- num/denom
  sig.e.sq.hat = max(sig.e.sq.hat, 10^(-4))
  ######estimate eta
  G=G.list
  d<-length(G)
  B.a = matrix(0,nrow=d,ncol=d)
  omega.hat<-rep(0,d)
  for (i in 1:n){
    cur.mem=which(grp==grp_uniq[i])
    zi = as.matrix(z[cur.mem,])
    ri<-r.hat[cur.mem]
    mi=length(cur.mem)
    Sig.ai.inv = solve(a*zi%*%t(zi)+ diag(1,mi)) #(Sig_a^i)^{-1}
    for (j in 1:d){
      for (k in 1:d){
        B.a[j,k] = B.a[j,k]+sum(diag(as.matrix(Sig.ai.inv%*%zi%*%G[[j]]%*%t(zi)%*%Sig.ai.inv%*%zi%*%G[[k]]%*%t(zi))))
      }
      
      omega.hat[j] = omega.hat[j] + t(ri)%*%Sig.ai.inv%*%zi%*%G[[j]]%*%t(zi)%*%Sig.ai.inv%*%ri-
        sig.e.sq.hat*sum(diag(Sig.ai.inv%*%zi%*%G[[j]]%*%t(zi)%*%Sig.ai.inv))
    }
  }
  eta.hat = solve(B.a)%*%omega.hat
  psi.hat<- sum(sapply(1:d, function(j) eta.hat[j]*G[[j]]))
  if(length(G)==1){
    eta.hat<-max(eta.hat,0)
  }
  is.spd<- (min(psi.hat)>=0)
  return(list(eta.hat = eta.hat, sig.e.sq.hat = sig.e.sq.hat, is.spd=is.spd))
}
