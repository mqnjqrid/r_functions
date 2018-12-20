source("C:/Users/manja/Dropbox/capture_recapture/codes/indep_cov_Tilling_simulation.R")

simuldraw = 3
nums = c(1000, 5000, 10000, 15000, 20000)
cntb = list()
cntv = list()
whichfunc = c(psihats_p, psihats_mlogit, psihats_np, psihats_gam)
funcname = c("P", "Mlogit", "NP", "Gam")
for (iter in 1:length(nums)) {
  print(iter)
  n = nums[iter]
  psi_bv = psi_bias_var(n, 2, 4, simuldraw, whichfunc)
  cntbias = psi_bv$psi
  cntvar = psi_bv$psivar
  cntbmat <- matrix(cntbias, ncol = 2*length(whichfunc), byrow = T)
  cntb[[iter]] = do.call(rbind,(lapply(1:length(whichfunc), function(fi) {return(cntbmat[,(2*fi):(2*fi - 1)])})))

  cntvmat <- matrix(cntvar, ncol = 2*length(whichfunc), byrow = T)
  cntv[[iter]] = do.call(rbind,(lapply(1:length(whichfunc), function(fi) {return(cntvmat[,(2*fi):(2*fi - 1)])})))

}
par(mfrow = c(2, length(nums)), las = 3)
ucap = 3
#colvec = c("cornflowerblue", "coral", "blue", "lightsalmon")
#colvec = c("cornflowerblue", "darkgoldenrod1", "blue", "darkgoldenrod")
#colvec = c("cadetblue", "aquamarine", "darkslategray", "darkolivegreen4")
#colvec = c("aquamarine", "seagreen2", "darkolivegreen4", "darkslategray")
#colvec = c("darkseagreen1", "springgreen3", "darkolivegreen4", "darkslategray")
colvec = c("darkolivegreen4", "darkseagreen1")#, "cadetblue", "aquamarine")
for(i in 1:length(nums)){
  counts = cntb[[i]]
  counts = pmin(counts, ucap)
  #counts[,4] = 0
  rownames(counts) = paste(rep(c("Cor", "Mis"), length(whichfunc)), rep(funcname, each = 2), sep =' ')
  colnames(counts) = c("bias proposed", "bias plugin")
  barplot((t(counts)), ylab = "Bias of psi_hat", main = paste("n = ", nums[i], sep = ''),
          xlab="Data", col = colvec,
          legend = colnames(counts), args.legend = list(x = "topright", bg = "transparent"),
          beside=TRUE)
}
for(i in 1:length(nums)){
  counts = sqrt(cntv[[i]] + (cntb[[i]])^2)
  counts = pmin(counts, ucap)
  rownames(counts) = paste(rep(c("Cor", "Mis"), length(whichfunc)), rep(funcname, each = 2), sep =' ')
  colnames(counts) = c("sqmse proposed", "sqmse plugin")
  barplot((t(counts[1:8,])), ylab = "RMSE of psi_hat", main = paste("n = ", nums[i], sep = ''),
          xlab="Data", col = colvec,
          legend = colnames(counts), args.legend = list(x = "topright", bg = "transparent"), beside=TRUE)
}


save(nums, simuldraw, cntb, cntv, file = "C:/Users/manja/Dropbox/capture_recapture/codes/bias_sqmse_2.Rdata")
#matbias = rbind(cntb[[1]], cntb[[2]], cntb[[3]], cntb[[4]], cntv[[1]], cntv[[2]], cntv[[3]], cntv[[4]])
#write.table(matbias, file = "C:/Users/manja/Dropbox/capture_recapture/data/psi_bias_sqmse_300_mlogit_nonls.txt", row.names = FALSE, col.names = FALSE)
