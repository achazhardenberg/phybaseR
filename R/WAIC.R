WAIC<-function(model,n.iter = 2000,n.burnin = 500,n.thin = 10){
  samples <- rjags::jags.samples(model$model,
                              c("WAIC","deviance"),
                              type = "mean",
                              n.iter = n.iter,
                              n.burnin = n.burnin,
                              n.thin = n.thin)

  samples$p_waic <- samples$WAIC
  samples$waic <- samples$deviance + samples$p_waic
  tmp <- sapply(samples, sum)
  waic <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
  return(waic)
}
