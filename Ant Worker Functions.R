




#Minimum ant workers function

Minimum.workers <- function(N, alpha, beta, c1, c2, c3, m, n){
  
R <- matrix(data = 0, nrow = 4, ncol = 1, byrow = FALSE, dimnames = NULL)
  #F
  R[1,1] <- sqrt((m+c1)/(alpha*beta))
  
  h <- (N-beta*(sqrt((c1 +m)/(alpha*beta))))
  #G
  R[2,1] <- (h*(n+c3-c2))/(N-beta*(sqrt((c1 +m)/(alpha*beta))))
  #Minimum Protein
  R[3,1] <- R[2,1]*(N-beta*R[1,1])/h +c2-c3
  #Minimum Sugar
  R[4,1] <- alpha*R[1,1]^2*beta -c1
 
  write.table(R, "F:/data1.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
  
}

#Maximum sugar function
Maximum.sugar <- function(H1, H2, F, alpha1, alpha2, beta, c1){
  
S <- matrix(data = 0, nrow = 4, ncol = 1, byrow = FALSE, dimnames = NULL)  

  if(alpha1 > alpha2){
    if(H1 < F*beta){
      #f1
      S[1,1] <- H1/beta
      #f2
      S[2,1] <- F-H1/beta
    }
    else{
      #f1
      S[1,1] <- F
      #f2
      S[2,1] <- 0
    }
  }
  else {
    if(H2 < F*beta){
      #f2
      S[2,1] <- H2/beta
      #f1
      S[1,1] <- F-H2/beta
    }
    else{
      #f2
      S[2,1] <- F
      #f1
      S[1,1] <- 0
    }

  }
  #Maximum Sugar
  S[3,1] <- alpha1*S[1,1]^2*beta + alpha2*S[2,1]^2*beta - c1

  write.table(S, "F:/data2.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
}

#Maximum protein function
Maximum.protein <- function(H1, H2, f1, f2, G, beta, h1, h2, c2, c3){

T <- matrix(data = 0, nrow = 4, ncol = 1, byrow = FALSE, dimnames = NULL) 

  if(((H1 - f1*beta)) > ((H2 - f2*beta))){
    if(H1 < G){
      #g1
      T[1,1] <- H1
      #g2
      T[2,1] <- G-H1 
    }
    else{
      #g1
      T[1,1] <- G
      #g2
      T[2,1] <- 0 
    }
  }
  else{
    if(H2 < G){
      #g2
      T[2,1] <- H2
      #g1
      T[1,1] <- G-H2 
    }
    else{
      #g2
      T[2,1] <- G
      #g1
      T[1,1] <- 0 
    }
  }
  
  #proteinMax
  T[3,1] <- T[1,1]*((H1-f1*beta)/h1) + T[2,1]*((H2-f2*beta)/h2) + c2 - c3
  
  write.table(T, "F:/data3.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
}

