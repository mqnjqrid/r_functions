library(glmnet)
library(nnet)
library(ranger)
library(gam)
library(mlogit)

source("C:/Users/manja/Dropbox/capture_recapture/codes/indep_cov_Tilling_simulation.R")

psihats_p = function(List_matrix_cov, K, i, j, simul){
  #      p1 = profvis({
  N = nrow(List_matrix_cov)
  eps = 0.005
  psiq = numeric(0)
  qnphi = numeric(0)
  for (s in 1:simul){
    set1 = sample(N, N/2, replace = FALSE)
    List1 = as.data.frame(List_matrix_cov[set1,])
    List2 = as.data.frame(List_matrix_cov[-set1,])
    
    fiti = glm(formula(paste("L", i, " ~.", sep = '')), family = binomial(link = "logit"), data = List1[,-c(1:K)[-i]])
    fitj = glm(formula(paste("L", j, " ~.", sep = '')), family = binomial(link = "logit"), data = List1[,-c(1:K)[-j]])
    fitij = glm(formula(paste("L", i, "*L", j, " ~.", sep = '')), family = binomial(link = "logit"), data = List1[,c(i, j, (K + 1):ncol(List1))])
    q1dot = pmax(predict(fiti, newdata = List2, type = "response"), eps)
    qdot1 = pmax(predict(fitj, newdata = List2, type = "response"), eps)
    q11 = pmax(predict(fitij, newdata = List2, type = "response"), eps)
    psihat = 1/mean(q1dot*qdot1/q11)
    psiq = c(psiq, psihat)
    yi = List2[,paste("L", i, sep = '')]
    yj = List2[,paste("L", j, sep = '')]
    Qnphihat = mean(-psihat^2 * ((yj - qdot1)*q1dot/q11
                                 + (yi - q1dot)*qdot1/q11 
                                 - (yi*yj - q11)*q1dot*qdot1/q11^2 + q1dot*qdot1/q11) + psihat)
    Qnphihat
    qnphi = c(qnphi, Qnphihat)
  }
 #    })
  return(list(psiq = mean(psiq), qnphi = mean(qnphi)))
}

psihats_mlogit = function(List_matrix_cov, K, i, j, simul){
  #      p1 = profvis({
  N = nrow(List_matrix_cov)
  eps = 0.005
  psiq = numeric(0)
  qnphi = numeric(0)
  for (s in 1:simul){
    set1 = sample(N, N/2, replace = FALSE)
    List1 = as.data.frame(List_matrix_cov[set1,])
    List2 = as.data.frame(List_matrix_cov[-set1,])
    Listy = cbind(List1[,"L1"] + 2*List1[,"L2"], List1[,-(1:K)])
    colnames(Listy) = c("choice", colnames(List_matrix_cov)[-c(1:K)])
    Listy = as.data.frame(Listy)
    mml_train = mlogit.data(Listy, choice = "choice", shape = "wide", alt.levels = 1:3)
    Listy = cbind(List2[,"L1"] + 2*List2[,"L2"], List2[,-(1:K)])
    colnames(Listy) = c("choice", colnames(List_matrix_cov)[-c(1:K)])
    Listy = as.data.frame(Listy)
    mml_test = mlogit.data(Listy, choice = "choice", shape = "wide", alt.levels = 1:3)
    
    mfit = try(mlogit(choice ~ 1 | x1, data = mml_train, outcome = FALSE))
    if(class(mfit) != "try-error"){
      coef = mfit$coefficients
      xlist = List2[,"x1"]
      l1 = rep(1, length(xlist));
      l2 = coef[1] + coef[3]*xlist;
      l3 = coef[2] - coef[4]*xlist
      ##    pred = predict(mfit, newdata = mml_test)
      q1dot = pmin(pmax(l1/(l1 + l2 + l3), eps), 1 - eps)
      qdot1 = pmin(pmax(l2/(l1 + l2 + l3), eps), 1 - eps)
      q11 = pmin(pmax(l3/(l1 + l2 + l3), eps), 1 - eps)
      psihat = 1/mean(q1dot*qdot1/q11)
      psiq = c(psiq, psihat)
      yi = List2[,paste("L", i, sep = '')]
      yj = List2[,paste("L", j, sep = '')]
      Qnphihat = mean(-psihat^2 * ((yj - qdot1)*q1dot/q11
                                   + (yi - q1dot)*qdot1/q11 
                                   - (yi*yj - q11)*q1dot*qdot1/q11^2 + q1dot*qdot1/q11) + psihat)
      qnphi = c(qnphi, Qnphihat)
    }
   
  }
  return(list(psiq = mean(psiq), qnphi = mean(qnphi)))
}

psihats_np = function(List_matrix_cov, K, i, j, simul){
  N = nrow(List_matrix_cov)
  eps = 0.005
  #      p1 = profvis({
  psiq = numeric(0)
  qnphi = numeric(0)
  for(s in simul){
    set1 = sample(N, ceiling(N/2), replace = FALSE)
    List1 = as.data.frame(List_matrix_cov[set1,])
    List2 = as.data.frame(List_matrix_cov[-set1,])
    #sapply((K + 1):ncol(List1), function(l){ paste(rep(paste('x', l, sep = ' + '), 3), c(' + ', '^2 + ', '^3 + '), sep = '') })
    fiti = ranger(formula(paste("factor(L", i, ") ~.", sep = '')), data = List1[,-c(1:K)[-i]], probability = TRUE)
    fitj = ranger(formula(paste("factor(L", j, ") ~.", sep = '')), data = List1[,-c(1:K)[-j]], probability = TRUE)
    fitij = ranger(formula(paste("factor(L", i, "*L", j, ") ~.", sep = '')), data = List1[,c(i, j, (K + 1):ncol(List1))], probability = TRUE)
    qpred10 = predict(fiti, data = List2)$predictions
    qpred01 = predict(fitj, data = List2)$predictions
    qpred11 = predict(fitij, data = List2)$predictions
    if('1' %in% colnames(qpred10)){
      q1dot = qpred10[,'1']
    }else {
      q1dot = eps
    }
    if('1' %in% colnames(qpred01)){
      qdot1 = qpred01[,'1']
    }else {
      qdot1 = eps
    }
    if('1' %in% colnames(qpred11)){
      q11 = qpred11[,'1']
    }else {
      q11 = eps
    }
    psihat = 1/mean(q1dot*qdot1/q11)
    psiq = c(psiq, psihat)
    yi = List2[,paste("L", i, sep = '')]
    yj = List2[,paste("L", j, sep = '')]
    Qnphihat = mean(-psihat^2 * ((yj - qdot1)*q1dot/q11
                                 + (yi - q1dot)*qdot1/q11 
                                 - (yi*yj - q11)*q1dot*qdot1/q11^2 + q1dot*qdot1/q11) + psihat)
    qnphi = c(qnphi, Qnphihat)
  }
  #    })
  return(list(psiq = mean(psiq), qnphi = mean(qnphi)))
}

psihats_gam = function(List_matrix_cov, K, i, j, simul){
  N = nrow(List_matrix_cov)
  List_matrix_cov = as.data.frame(List_matrix_cov)
  eps = 0.005
  #      p1 = profvis({
  psiq = numeric(0)
  qnphi = numeric(0)
  for(s in simul){
    set1 = sample(N, ceiling(N/2), replace = FALSE)
    List1 = as.data.frame(List_matrix_cov[set1,])
    List2 = as.data.frame(List_matrix_cov[-set1,])
    
    #sapply((K + 1):ncol(List1), function(l){ paste(rep(paste('x', l, sep = ' + '), 3), c(' + ', '^2 + ', '^3 + '), sep = '') })
    fiti = gam(formula(paste("L", i, " ~", paste("s(x", 1:l, ",df = 2)", sep = '', collapse = ' + '), sep = '')), data = List1, family = gaussian())
    fitj = gam(formula(paste("L", j, " ~", paste("s(x", 1:l, ",df = 2)", sep = '', collapse = ' + '), sep = '')), data = List1, family = gaussian())
    fitij = gam(formula(paste("L", i, "*L", j, " ~", paste("s(x", 1:l, ",df = 2)", sep = '', collapse = ' + '), sep = '')), data = List1, family = gaussian())
    q1dot = pmax(pmin(predict(fiti, newdata = List2), 0.90), 0.10)
    qdot1 = pmax(pmin(predict(fitj, newdata = List2), 0.90), 0.10)
    q11 = pmax(pmin(predict(fitij, newdata = List2), 0.90), 0.10)
    
    psihat = 1/mean(q1dot*qdot1/q11)
    psiq = c(psiq, psihat)
    yi = List2[,paste("L", i, sep = '')]
    yj = List2[,paste("L", j, sep = '')]
    Qnphihat = mean(-psihat^2 * ((yj - qdot1)*q1dot/q11
                                 + (yi - q1dot)*qdot1/q11 
                                 - (yi*yj - q11)*q1dot*qdot1/q11^2 + q1dot*qdot1/q11) + psihat)
    qnphi = c(qnphi, Qnphihat)
  }
  #    })
  return(list(psiq = mean(psiq), qnphi = mean(qnphi)))
}

estim_pnp = function(List_matrix, n, K, whichfunc){
  if(missing(n)){
    n = nrow(List_matrix)
  }
  if(missing(K)){
    K = ncol(List_matrix) - 1
  }
  
  #removing all rows with only 0's
  List_matrix_cov = List_matrix[which(rowSums(List_matrix[,1:K]) > 0),]
  colnames(List_matrix_cov) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List_matrix_cov) - K), sep = ''))
  #N = number of observed or captured units
  N = nrow(List_matrix_cov)
  simul = 100
  
  output_summary = matrix(NA, nrow = K*(K - 1)/2, ncol = 2*length(whichfunc))
  rownames(output_summary) = unlist(sapply(1:(K - 1), function(k) {
    sapply((k + 1):K, function(s) {
      return(paste(k, ", ", s, sep = ''))
    })}))                                                          
  #colnames(output_summary) =  c("phiq_p", "phihat_p", "phiq_np", "phihat_np", "phiq_gam", "phihat_gam")
  
  for(i in 1:(K - 1)){
    for(j in (i + 1):K){
      for(fi in 1:length(whichfunc)){
        estims = whichfunc[[fi]](List_matrix_cov, K, i, j, simul)
        psiq = unlist(estims$psiq)
        qnphi = unlist(estims$qnphi)
        output_summary[paste(i, ", ", j, sep = ''),(2*fi - 1):(2*fi)] = c(mean(psiq), mean(psiq + qnphi))
      }

    }
  }
  return(list(output_summary = output_summary, psi = N/n))
}

psi_bias_var = function(nsize, K, l, simuldraw, whichfunc){
  psi_all = matrix(NA, nrow = simuldraw, ncol = 4*length(whichfunc))

  
  for(s in 1:simuldraw){
    print(s)
    datap = dat_p(nsize, l)
    List_matrix = datap$List_matrix
    List_matrix_xstar = datap$List_matrix_xstar
    est_val = estim_pnp(List_matrix, K = 2, whichfunc = whichfunc)
    psi_all[s,1:(2*length(whichfunc))] = est_val$output_summary - est_val$psi
    est_val = estim_pnp(List_matrix_xstar, K = 2, whichfunc = whichfunc)
    psi_all[s,(2*length(whichfunc)) + 1:(2*length(whichfunc))] = est_val$output_summary - est_val$psi
    psi_all[1,]
  }
  psi_12 = colMeans(abs(psi_all), na.rm = TRUE)
  psivar_12 = colMeans(psi_all*psi_all, na.rm = TRUE) -
              colMeans(psi_all, na.rm = TRUE)^2
  return(list(psi = psi_12, psivar = psivar_12))
}
