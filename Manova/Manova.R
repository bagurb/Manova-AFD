
library(rstatix)
MANOVA_class <- function(dfw)
{
  #-> détection automatique des variables numériques
  id_num <- which(sapply(dfw, is.numeric))
  X <- dfw[,id_num]
  
  #-> .. on cherche la variable NON numérique ..
  id_string  <- which(!sapply(dfw, is.numeric))
  
  #-> pour faciliter et uniformiser les recherches des classes, on renomme cette variable
  names(dfw)[id_string] <- 'fac' 
  
  # ici code à compléter en utilisant ce qui a été développé dans la première partie
  
  #  création de la liste
  
  L <- lapply(levels(dfw$fac), function(x){ return(list('data' = subset(dfw,fac==x),
                                                         mean = (1/nrow(subset(dfw,fac==x)))*matrix(1:1,ncol=nrow(subset(dfw,fac==x)),nrow=1,byrow=FALSE)%*%as.matrix(cbind(subset(dfw,fac==x)[,1:ncol(subset(dfw,fac==x))-1])),
                                                         s = (1/(nrow(subset(dfw,fac==x))-1))*(t(as.matrix(cbind(subset(dfw,fac==x)[,1:ncol(subset(dfw,fac==x))-1])))%*%as.matrix(cbind(subset(dfw,fac==x)[,1:ncol(subset(dfw,fac==x))-1])) - t(as.matrix(cbind(subset(dfw,fac==x)[,1:ncol(subset(dfw,fac==x))-1])))%*%((1/(nrow(subset(dfw,fac==x))))*matrix(1:1,ncol=nrow(subset(dfw,fac==x)),nrow=nrow(subset(dfw,fac==x)),byrow=FALSE))%*%as.matrix(cbind(subset(dfw,fac==x)[,1:ncol(subset(dfw,fac==x))-1]))),                                                      
                                                         n = nrow(subset(dfw,fac==x)) ))  })
  L <- setNames(L,c("G1","G2"))
  # calcul de la variance commune
  VC = box_m(data=dfw[,1:ncol(dfw)-1], group=dfw$fac)
    
    Sd2 <- as.matrix(((L$G1$n - 1)*(L$G1$s) +(L$G2$n - 1)*(L$G2$s))/((L$G1$n+L$G2$n) - 2))
    
    # matrice des différences des moyennes
    Y <- matrix(L$G1$mean - L$G2$mean)
    Yt <- t(Y)
    
    # calcul du T2 de Hotteling
    T2 <- Yt%*%solve(((1/L$G1$n)+(1/L$G2$n))*Sd2)%*%(Y)
    
    # calcul de la p value
    Fobs <- ((L$G1$n + L$G2$n - (ncol(dfw)-1) - 1)/((L$G1$n + L$G2$n - 2)*(ncol(dfw)-1)))*T2
    pvalue <- df(Fobs,(ncol(dfw)-1),L$G1$n + L$G2$n - (ncol(dfw)-1) - 1)
    
  # les résultats sont stockés dans une liste
  
  test_result        <- data.frame('T² Hotteling' = T2, 'Fobs' = Fobs, 'p value' = pvalue )
  variance           <- data.frame(Sd2) 
  colnames(variance) <- names(df)[id_num] ; rownames(variance) <- names(df)[id_num]
  class_mean         <- data.frame(t(sapply(L, function(x){return(x$mean)})))
  
  ret   <- list(test = test_result, variance = variance, class_mean = class_mean)
  
  return(ret)
}

df2           <- data.frame(subset(iris, iris$Species != 'virginica'))
df2$Species  <- factor(df2$Species )

result <- MANOVA_class(df2)
result
