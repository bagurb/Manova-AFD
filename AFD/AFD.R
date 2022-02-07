rm(list = ls())
library(formattable) # formatage des résultats sous forme de tableau widget
df  <- read.table('VIN_QUALITE.txt', header = T)
df$Quality <- factor(df$Quality)

#-> .. on cherche la variable NON numérique ..
id_string  <- which(!sapply(df, is.numeric))

#-> pour faciliter et uniformiser les recherches des classes, on renomme cette variable
names(df)[id_string] <- 'fac' 

#LIST FOR General
mG <- list(TP = mean(df$TP),Sun = mean(df$Sun),Heat = mean(df$Heat),Rain = mean(df$Rain))

#LIST FOR Class 1
data <- subset(df,df$fac=='bad')
n1 <- nrow(data)
mC1 <- list(mean = list(TP= mean(data$TP),Sun= mean(data$Sun),Heat= mean(data$Heat),Rain= mean(data$Rain)), var = var(data[,1:4])*((nrow(data[,1:4])-1)/nrow(data[,1:4])))

#LIST FOR Class 2
data <- subset(df,df$fac=='medium')
n2 <- nrow(data)
mC2 <- list(mean = list(TP= mean(data$TP),Sun= mean(data$Sun),Heat= mean(data$Heat),Rain= mean(data$Rain)), var = var(data[,1:4])*((nrow(data[,1:4])-1)/nrow(data[,1:4])))

#LIST FOR Class 3
data <- subset(df,df$fac=='good')
n3 <- nrow(data)
mC3 <- list(mean = list(TP= mean(data$TP),Sun= mean(data$Sun),Heat= mean(data$Heat),Rain= mean(data$Rain)), var = var(data[,1:4])*((nrow(data[,1:4])-1)/nrow(data[,1:4])))

#LIST FOR ALL STATS

L <- list(G = mG,C1 = mC1,C2 = mC2,C3 = mC3)

L$C1$var

#INERTIE INTRA#

N <- nrow(df)
W <- (1/N)*(n1*L$C1$var+n2*L$C2$var+n3*L$C3$var)
W <- as.matrix(W)
W
#INERTIE INTER#

B <- (1/N)*((n1*(unlist(L$C1$mean)-unlist(L$G))%*%t(unlist(L$C1$mean)- unlist(L$G)))+(n2*(unlist(L$C2$mean)- unlist(L$G))%*%t(unlist(L$C2$mean)-unlist(L$G)))+(n3*(unlist(L$C3$mean)-unlist(L$G))%*%t(unlist(L$C3$mean)- unlist(L$G))))
B

#INERTIE TOTALE#

Tot2 <- W+B
Tot2


#DIAGONALISATION#

#Methode de Fisher#
k <- 3

B2 <- B * (N/(k-1))
B2

W2 <- W * (N/(N-k))
W2


rapportBW <- as.matrix(B2%*%solve(W2))
rapportBW

#VALEUR PROPRE#
u <- eigen(rapportBW)
uReal <- Re(u$vectors)

reduxUvectors <- uReal[,1:2]
reduxUvalues <- Re(u$values[1:2])
reduxUvalues
reduxUvectors

normeVal <- sqrt(diag(t(reduxUvectors)%*%W%*%reduxUvectors))
normeVal
normeVec <- reduxUvectors/normeVal

normeVal
normeVec

VecP <- sweep(reduxUvectors,2,normeVal,'/') 
VecP

#INFERENCE STATISTIQUE#
p <- 4
#AXE 1#
Iw1 = reduxUvalues[1]/sum(reduxUvalues)
Iw1

Ncorr1 = sqrt(reduxUvalues[1]/(1 + reduxUvalues[1]))
Ncorr1

A1 = prod(1 - Ncorr1**2)
A1

E1 = -(N-((p+k)/2)-1)*log(A1)
E1
#Calcul du test de Chi2#
#AXE 2#
Iw2 = reduxUvalues[2]/sum(reduxUvalues)
Iw2

Ncorr2 = sqrt(reduxUvalues[2]/(1 + reduxUvalues[2]))
Ncorr2

A2 = prod(1 - Ncorr2**2)
A2

E2 = -(N-((p+k)/2)-1)*log(A2)
E2

L = list(Axe1 = list(Iw = Iw1, Ncorr = Ncorr1),Axe2 = list(Iw = Iw2,Ncorr = Ncorr2))


## PART 2##
N = nrow(df)
df_matrix = unname(as.matrix(sapply(df[1:4],as.numeric)))
mg_matrix = matrix(unlist(mG), nrow=N, ncol = 4, byrow = TRUE)

Z <- as.matrix(df_matrix - mg_matrix)
Z

pre_score = Z%*%VecP
score = cbind(pre_score,df[5])

names(score)<-c("Axe_1","Axe_2","fac")
score

names(score)[3] <- 'class' <- 'class'
names(score)
library(ggplot2)
gr_ind <- ggplot() +  geom_point(data = score  ,aes( x = Axe_1, y = Axe_2, colour = class, shape = class))                                # trace les points
gr_ind <- gr_ind   +  geom_hline(yintercept = 0 , size = 0.1, colour = '#CCCCCC') + geom_vline(xintercept  = 0, size = 0.1,colour = '#CCCCCC') # trace les lignes
gr_ind <- gr_ind   +  xlim(c(-5,4)) + ylim(c(-3,3))  # uniformités des échelles
gr_ind

#3.1
Axe_1 = as.numeric(score[,1])
Axe_2 = as.numeric(score[,2])

G <- cbind(pre_score,df[5])

G_fac <- aggregate(cbind(Axe_1,Axe_2),list(df$fac),mean)

G_fac

gr_ind <- ggplot() +  geom_point(data = G_fac  ,aes( x = Axe_1, y = Axe_2, colour = Group.1, shape =  Group.1), size = 3)  +  # trace les centres de gravité
  geom_hline(yintercept = 0 , size = 0.1, colour = '#CCCCCC') + geom_vline(xintercept  = 0, size = 0.1,colour = '#CCCCCC') + # trace les lignes
  xlim(c(-5,4)) + ylim(c(-3,3)) 
gr_ind


#3.2
G <- as.matrix(subset(G_fac,select=-Group.1))
G

distance <- NULL

for (i in 1:nrow(G))
{ delta <- score[1:2] - matrix(G[i,],ncol = ncol(G), nrow = N, byrow = T)
distance <- cbind(distance,apply(delta,1,function(x){return(sqrt(sum(x*x)))}))
}

df_distance <- as.data.frame(distance) ; names(df_distance) <- G_fac[,1]
df_distance

Quality_pred <- factor(apply(df_distance,1,function(x){return(names(x)[which.min(x)])})) 
Quality <- head(df$fac, "ind")

Qtab <- as.data.frame(cbind(Quality= as.character(Quality),Quality_pred = as.character(Quality_pred)))
Qtab
#3.3
Reference <- table(Quality_pred,Quality)
Reference
confusionMatrix(Reference)

G <- cbind(pre_score,df[5])

as.data.frame(G)
Scores <- cbind(score, good_class = ifelse(Qtab$Quality_pred ==  Qtab$Quality, 'good','bad') )

gr_ind <- ggplot() + geom_point(data = Scores, aes(x = Axe_1, y = Axe_2, colour = class), size = 2)    +
  geom_text (data = Scores, aes(x = Axe_1, y = Axe_2, label = good_class, colour = good_class), hjust = 1, vjust = 0 ) +
  geom_point(data = G_fac,  aes(x = Axe_1, y = Axe_2, colour = Group.1, shape =  Group.1), size = 3)  +  # trace les centres de gravité
  geom_hline(yintercept = 0 , size = 0.1, colour = '#CCCCCC') + geom_vline(xintercept  = 0, size = 0.1,colour = '#CCCCCC') 
gr_ind  

