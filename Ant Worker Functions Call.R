# Function Call

# Plot J throughout to see behavior of ants
# What starting density does red aphid species have to be at in order to succeed the black species?
# What growth rate?
# Start thinking of figures for poster.

# Changes/differences in alpha = f(handling time, sugar quality)
# Changes in ant nutritional requirements. Nutritional requirement equivalent to colony size.


Q <- matrix(data = NA, nrow = 3, ncol = 100, byrow = FALSE, dimnames = NULL)
B <- matrix(data = NA, nrow = 1, ncol = 100, byrow = FALSE, dimnames = NULL)
A <- matrix(data = NA, nrow = 8, ncol = 100, byrow = FALSE, dimnames = NULL)
R <- matrix(data = NA, nrow = 2, ncol = 100, byrow = FALSE, dimnames = NULL)
U <- matrix(data = NA, nrow = 2, ncol = 100, byrow = FALSE, dimnames = NULL)
D <- matrix(data = NA, nrow = 1, ncol = 100, byrow = FALSE, dimnames = NULL)
######################

#m and n change
c1 <- 30
c2 <-10

m <- 30
n <- 10
Ssurp <- 0
Psurp <- 0
D[1,] <- "30:10"

F <- 0
G <- 0

ws <- .01
wp <- .01

dp1 <- .2
dp2 <- .8
dp3 <- 1
dp <- (dp1+dp2+dp3)/3
B[1,] <- ".6"

ds1 <- dp1
ds2 <- dp2
ds3 <- 1
ds <- (ds1+ds2+ds3)/3

s1 <- 1
s2 <- 1
s3 <- 1
s <- (s1+s2+s3)/3


beta <-5
cs3 <- 5
Cs <- 5
  
cp3 <- 1
Cp <- 1

#######################

#L1, L2, w1, w2, a, and k are fixed values
L1 <-2
L2 <-2
w1 <-1
w2 <-1
a <- .01
k <- 1

#Initial host 1
Q[1,1] <- 100
#Initial host 2
Q[2,1] <-100
#Initial parasitoid
Q[3,1] <- 10

for(i in 1:99){
  
  J <- matrix(data = 0, nrow = 8, ncol = 3, byrow = FALSE, dimnames = NULL)  
  
  
  H1 <- Q[1,i]
  H2 <- Q[2,i]
  Avg <- (H1 + H2 + 10)/3
  P <- Q[3,i]
  
  
  Minimum.workers_2(Avg, s, ds, dp, ws, wp, G, F, m, n, Cs)
  tbl1 <- read.table("F:/data1.txt")
  J[,1] <- tbl1[,1]
  file.remove("F:/data1.txt")
  F <- J[1,1]
  G <- J[2,1]

  
  
  Maximum.nutrients(H1,H2, F, G, beta, ws, wp, s1, s2, s3, dp1, dp2, dp3, ds1, ds2, ds3)
  tbl2 <- read.table("F:/data2.txt")
  J[,2] <- tbl2[,1]
  file.remove("F:/data2.txt")

  
  Nutrient.adjust(c1, c2, m, n, J[7,2], J[8,2], Ssurp, Psurp)
  tbl3 <- read.table("F:/data3.txt")
  J[,3] <- tbl3[,1]
  file.remove("F:/data3.txt")
  m <- J[2,3]
  Ssurp <- J[1,3]
  n <- J[4,3]
  Psurp <- J[3,3]
 
  
  A[,i] <- J[,2]
  #########################################################################
  
  Q[1,i+1] <- L1*J[1,2]*beta + L1*((H1-J[1,2]*beta) - (J[4,2]*(1-dp1+wp)))*((1+ (a*P/k))^-k)
  Q[2,i+1] <- L2*J[2,2]*beta + L2*((H2-J[2,2]*beta) - (J[5,2]*(1-dp2+wp)))*((1+ (a*P/k))^-k)
  Q[3,i+1] <- w1*L1*((H1-J[1,2]*beta) - (J[4,2]*(1-dp1+wp)))*(1-((1+ (a*P/k))^-k)) + w2*L2*((H2-J[2,2]*beta) - (J[5,2]*(1-dp2+wp)))*(1-((1+ (a*P/k))^-k)) + 1
  
  R[1,i+1] <- (L1*((H1-J[1,2]*beta) - (J[4,2]*(1-dp1+wp)))*(1-((1+ (a*P/k))^-k)) +1)/((L1*J[1,2]*beta + L1*((H1-J[1,2]*beta) - (J[4,2]*(1-dp1+wp)))*((1+ (a*P/k))^-k)) + L1*((H1-J[1,2]*beta) - (J[4,2]*(1-dp1+wp)))*(1-((1+ (a*P/k))^-k)) +1)
  R[2,i+1] <- (L2*((H2-J[2,2]*beta) - (J[5,2]*(1-dp2+wp)))*(1-((1+ (a*P/k))^-k))+1)/((L2*J[2,2]*beta + L2*((H2-J[2,2]*beta) - (J[5,2]*(1-dp2+wp)))*((1+ (a*P/k))^-k)) + L2*((H2-J[2,2]*beta) - (J[5,2]*(1-dp2+wp)))*(1-((1+ (a*P/k))^-k))+1)
  U[1,i+1] <-(J[1,2] - J[4,2])/(J[1,2]+J[2,2]+J[3,2]+J[4,2]+J[5,2]+J[6,2] + 1e-15)
  U[2,i+1] <-(J[2,2] - J[5,2])/(J[1,2]+J[2,2]+J[3,2]+J[4,2]+J[5,2]+J[6,2] + 1e-15) 
}
# df6 <- data.frame(D[1,])
# df5 <- data.frame(B[1,])
# df1 <- data.frame(R[1,])
# df2 <- data.frame(R[2,])
# df3 <- data.frame(U[1,])
# df4 <- data.frame(U[2,])
# df6D <- cbind(df6,df5,df1,df2,df3,df4)
# colnames(df6D) <- c("Nutrient_Reqs","Distance","H1_Par", "H2_Par", "H1_Ants", "H2_Ants")
df1 <- data.frame(Q[1,])
id1 <- rownames(df1)
df1 <- cbind(id = id1, df1)
df2 <- data.frame(Q[2,])
df3 <- data.frame(Q[3,])

gg <- ggplot(data=df1, aes(x = id1)) +
  geom_line(aes(y = log10(Q[1,]), group=1, color="Host 1")) +
  geom_line(aes(y = log10(Q[2,]), group=1, color="Host 2")) +
  geom_line(aes(y = log10(Q[3,]), group=1, color="Parasitoid")) 
gg + ggtitle("Sample Population Trajectories") +
  theme_set(theme_classic(base_size = 24)) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(colour = "Species", x = "Time (Generations)", y = "Log Population Density")

# plot( log10(Q[1,]), main= "d",xlab="Initial Frequency",ylab="Frequency", ylim = c(-1,4), type ='l')
# points (log10(Q[2,]),main= "d",xlab="Initial Frequency",ylab="Frequency",col = "red", type = 'l')
# points (log10(Q[3,]),main= "d",xlab="Initial Frequency",ylab="Frequency",col = "blue", type = 'l')

# plot(mean(U[1,]), mean(R[1,]), main= "d",xlab="Initial Frequency",xlim = c(-1,1), ylim = c(0,1), ylab="Frequency")
# arrows(mean(U[1,]), mean(R[1,])-sd(R[1,]), mean(U[1,]), mean(R[1,])+sd(R[1,]), length=0.05, angle=90, code=3)
# arrows(mean(U[1,])-sd(U[1,]), mean(R[1,]),mean(U[1,])+sd(U[1,]) ,mean(R[1,]) , length=0.05, angle=90, code=3)

 #boxplot( U[1,], U[2,])
# points(mean(R), col="red")
# boxplot(R[2,], ylim=c(0,1))
# points(mean(R[2,]), col="red")
# boxplot(U[1,], ylim=c(-1,1))
# points(mean(U[1,]), col="red")
# boxplot(U[2,], ylim=c(-1,1))
# points(mean(U[2,]), col="red")

# layout(t(1:4))
# par(oma=c(2, 4, 4, 0), mar=rep(1, 4), cex=1)
# 
# boxplot(df1A$H2_Par~df1A$Distance, ylim=c(-1,1))
# points(mean(df1A$H1_Par), col="red")
# boxplot(df1A$H2_Par~df1A$Distance, ylim=c(-1,1))
# points(mean(df1A$H1_Par), col="red")

# points(mean(U[2,]), mean(R[2,]), main= "d",xlab="Initial Frequency",xlim = c(-1,1), ylim = c(0,1), ylab="Frequency", col = "red")
# arrows(mean(U[2,]), mean(R[2,])-sd(R[2,]), mean(U[2,]), mean(R[2,])+sd(R[2,]), length=0.05, angle=90, code=3)
# arrows(mean(U[2,])-sd(U[2,]), mean(R[2,]),mean(U[2,])+sd(U[2,]) ,mean(R[2,]) , length=0.05, angle=90, code=3)

# plot( (R[1,]), main= "d",xlab="Initial Frequency",ylab="Frequency", ylim = c(0,1), type ='l')
# points ((R[2,]),main= "d",xlab="Initial Frequency",ylab="Frequency",col = "red", type = 'l')
