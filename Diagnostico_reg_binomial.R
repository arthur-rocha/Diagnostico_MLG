#Pacotes (todos necessários para o funcionamento das funções)
library(ggplot2)
library(dplyr)
library(cowplot)
library(gridExtra)
library(hnp)

#------------Função de diagnóstico geral (binomial)--------------------

diag_binom=function(fit.model){
  X <- model.matrix(fit.model)
  n <- nrow(X)
  p <- ncol(X)
  w <- fit.model$weights
  W <- diag(w)
  H <- solve(t(X)%*%W%*%X)
  H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h <- diag(H)
  ts <- resid(fit.model,type="pearson")/sqrt(1-h)
  td <- resid(fit.model,type="deviance")/sqrt(1-h)
  di <- (h/(1-h))*(ts^2)
  a <- max(td)
  b <- min(td)
  x11(title = paste("Diagnosticos do",deparse(substitute(fit.model))))
  par(mfrow=c(2,2))
  plot(fitted(fit.model),h,xlab="Valores Ajustados", ylab="Medida h",
       main="Pontos de Alavanca", pch=16,panel.first=grid(col = rgb(.2,.2,.2,.1),lty=1),
       frame.plot = F,col.axis="gray70")
  axis(1,col = "gray70",col.ticks ="gray70")
  axis(2,col = "gray70",col.ticks ="gray70")
  # identify(fitted(fit1.model), h, n=1)
  #
  plot(di,xlab="Índice", ylab="Distância de Cook",
       main="Pontos Influentes",pch=16,panel.first=grid(col = rgb(.2,.2,.2,.1),lty=1),
       frame.plot = F,col.axis="gray70")
  axis(1,col = "gray70",col.ticks ="gray70")
  axis(2,col = "gray70",col.ticks ="gray70")
  # identify(di, n=1)
  #
  plot(td,xlab="Índice", ylab="Resíduo Componente do Desvio",
       main="Pontos Aberrantes", ylim=c(b-1,a+1), pch=16,
       panel.first=grid(col = rgb(.2,.2,.2,.1),lty=1),
       frame.plot = F,col.axis="gray70")
  axis(1,col = "gray70",col.ticks ="gray70")
  axis(2,col = "gray70",col.ticks ="gray70")
  abline(2,0,lty=2,col="orange")
  abline(-2,0,lty=2,col="orange")
  # identify(td, n=1)
  #
  plot(predict(fit.model),td,xlab="Preditor Linear", 
       ylab="Residuo Componente do Desvio",
       main="Função de Ligação", ylim=c(b-1,a+1), pch=16,
       panel.first=grid(col = rgb(.2,.2,.2,.1),lty=1),
       frame.plot = F,col.axis="gray70")
  axis(1,col = "gray70",col.ticks ="gray70")
  axis(2,col = "gray70",col.ticks ="gray70")
  abline(2,0,lty=2,col="orange")
  abline(-2,0,lty=2,col="orange")
  
  #envelope
  x11(title = paste("Envelope do",deparse(substitute(fit.model))))
  par(mfrow=c(1,1))
  fit.model %>% hnp(halfnormal = F,xlab="Quantis Teóricos",
                    ylab='Componente de Desvio',pch=16)
  grid(lty=1)
  
  #resultados
  cat("-------------------Resultados do modelo-------------------")
  sumario=fit.model %>% summary()
  x11(title = paste("Resultados do",deparse(substitute(fit.model))))
  grid.arrange(bottom=paste(sumario$call[2],""),top=paste("Resultados do",deparse(substitute(fit.model))),
               tableGrob(rbind(sumario$coefficients,
                               deviance=c(sumario$deviance,NA,NA,NA),
                               AIC=c(sumario$aic,NA,NA,NA))%>% round(digits = 3)))
  
  
}

#---------------- Fim função diagnóstico ------------------------------
