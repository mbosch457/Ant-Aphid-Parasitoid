#Minimum ant workers function

Minimum.workers_2 <- function(N, s, ds, dp, ws, wp, G, F, m, n, Cs, Cp){
  
  R <- matrix(data = 0, nrow = 4, ncol = 1, byrow = FALSE, dimnames = NULL)
  #F
  R[1,1] <- sqrt(((1-dp+wp)^(-1/2)*G + m)/(s*Cs*(1-ds+ws)^2))
  
  #G
  R[2,1] <- n/((1-dp+wp)*(N-(R[1,1]*beta)/3)) 
  
  #Minimum Sugar
  R[4,1] <- s*(1-ds+ws)^2*R[1,1]^2*Cs - ((1-dp+wp)^(-1/2))*R[2,1]
  
  #Minimum Protein
  R[3,1] <- (1-dp+wp)*R[2,1]*(N-(R[1,1]*beta)/3)
  
  write.table(R, "F:/data1.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
  
}




Maximum.nutrients <- function(H1, H2, F, G, beta, ws, wp, s1, s2, s3, dp1, dp2, dp3, ds1, ds2, ds3){
  
  S <- matrix(data = 0, nrow = 8, ncol = 1, byrow = FALSE, dimnames = NULL)  
  
  if(s1*(1-ds1+ws)^2 > s2*(1-ds2+ws)^2 && s2*(1-ds2+ws)^2 > s3*(ws)^2){
    if(H1 < F*beta){
      #f1
      S[1,1] <- H1/beta
      if(H2/beta < F-H1/beta){
        #f2
        S[2,1] <- H2/beta
        #f3
        S[3,1] <- F - H1/beta - H2/beta
      }
      else{
        
        S[2,1] <- F-H1/beta
        S[3,1] <- 0
      }
    }
    else{
  
      S[1,1] <- F
      S[2,1] <- 0
      S[3,1] <- 0
    }
  }
    
  else if(s2*(1-ds2+ws)^2 > s1*(1-ds1+ws)^2 && s1*(1-ds1+ws)^2 > s3*(ws)^2){
    if(H2 < F*beta){
      
      S[2,1] <- H2/beta
      
      if(H1/beta < F-H2/beta){
        
        S[1,1] <- H1/beta
        
        S[3,1] <- F - H2/beta - H1/beta
      }
      else{
        
        S[1,1] <- F-H2/beta
        S[3,1] <- 0
      }
    }
    else{
      
      S[1,1] <- 0
      S[2,1] <- F
      S[3,1] <- 0
    }

  }
  else if(s1*(1-ds1+ws)^2 > s3*(ws)^2  && s3*(ws)^2 > s2*(1-ds2+ws)^2){
    if(H1 < F*beta){
      
      S[1,1] <- H1/beta
      S[2,1] <- 0
      S[3,1] <- F -  H1/beta
      
    }
    else{
      
      S[1,1] <- F
      S[2,1] <- 0
      S[3,1] <- 0
    }
  }
  else if(s3*(ws)^2 > s1*(1-ds1+ws)^2  && s1*(1-ds1+ws)^2 > s2*(1-ds2+ws)^2){
      
      S[1,1] <- 0
      S[2,1] <- 0
      S[3,1] <- F
    
  }
  else if(s3*(ws)^2 > s2*(1-ds2+ws)^2  && s2*(1-ds2+ws)^2 > s1*(1-ds1+ws)^2){
    
    S[1,1] <- 0
    S[2,1] <- 0
    S[3,1] <- F
    
  }
  else if(s2*(1-ds2+ws)^2 > s3*(ws)^2  && s3*(ws)^2 > s1*(1-ds1+ws)^2){
    if(H2 < F*beta){
      
      S[1,1] <- 0
      S[2,1] <- H2/beta
      S[3,1] <- F -  H2/beta
      
    }
    else{
      
      S[1,1] <- 0
      S[2,1] <- F
      S[3,1] <- 0
    }
  }
  else if(s1*(1-ds1+ws)^2 == s2*(1-ds2+ws)^2 && s2*(1-ds2+ws)^2 > s3*(ws)^2){
    if(H1 + H2 < F*beta){
      #f1
      S[1,1] <- H1/beta
      S[2,1] <- H2/beta
      S[3,1] <- F - H1/beta - H2/beta
    }
    else{
      S[1,1] <- F/2
      S[2,1] <- F/2
      S[3,1] <- 0
    }
  }
  else if(s2*(1-ds2+ws)^2 == s3*(ws)^2 && s3*(ws)^2 > s1*(1-ds1+ws)^2){
    if(H2 < F*beta){
      
      S[1,1] <- 0
      S[2,1] <- H2/beta
      S[3,1] <- F - H2/beta
    }
    else{
      S[1,1] <- 0
      S[2,1] <- F/2
      S[3,1] <- F/2
    }
  }
  else if(s1*(1-ds1+ws)^2 == s3*(ws)^2 && s3*(ws)^2 > s2*(1-ds2+ws)^2 ){
    if(H1 < F*beta){
      
      S[1,1] <- H1/beta
      S[2,1] <- 0
      S[3,1] <- F - H1/beta
    }
    else{
      S[1,1] <- F/2
      S[2,1] <- 0
      S[3,1] <- F/2
    }
  }
  else if(s1*(1-ds1+ws)^2 > s2*(1-ds2+ws)^2 && s2*(1-ds2+ws)^2 == s3*(ws)^2 ){
    if(H1 < F*beta){
      #f1
      S[1,1] <- H1/beta
      if(H2/beta < F-H1/beta){
        #f2
        S[2,1] <- H2/beta
        #f3
        S[3,1] <- F - H1/beta - H2/beta
      }
      else{
        
        S[2,1] <- (F-H1/beta)/2
        S[3,1] <- (F-H1/beta)/2
      }
    }
    else{
      
      S[1,1] <- F
      S[2,1] <- 0
      S[3,1] <- 0
    }
  }
  else if(s2*(1-ds2+ws)^2 > s1*(1-ds1+ws)^2 && s1*(1-ds1+ws)^2 == s3*(ws)^2 ){
    if(H2 < F*beta){
      #f1
      S[2,1] <- H2/beta
      if(H1/beta < F-H2/beta){
        #f2
        S[1,1] <- H1/beta
        #f3
        S[3,1] <- F - H1/beta - H2/beta
      }
      else{
        
        S[1,1] <- (F-H2/beta)/2
        S[3,1] <- (F-H2/beta)/2
      }
    }
    else{
      
      S[1,1] <- 0
      S[2,1] <- F
      S[3,1] <- 0
    }
  }
  else if(s3*(ws)^2 > s2*(1-ds2+ws)^2 && s2*(1-ds2+ws)^2 == s1*(1-ds1+ws)^2 ){
   
      S[1,1] <- 0
      S[2,1] <- 0
      S[3,1] <- F
  }
  
  if((1-dp1+wp)*(H1-(S[1,1]*beta)) > (1-dp2+wp)*(H2-(S[2,1]*beta)) && (1-dp2+wp)*(H2-(S[2,1]*beta)) > wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 )){
    
    if((H1-(S[1,1]*beta)) < G){
      #g1
      S[4,1] <- (H1-(S[1,1]*beta))
      if((H2-(S[2,1]*beta)) < G - (H1-(S[1,1]*beta))){
        #g2
        S[5,1] <- (H2-(S[2,1]*beta))
        #g3
        S[6,1] <- G - (H1-(S[1,1]*beta))- (H2-(S[2,1]*beta))
      }
      else{
        
        S[5,1] <- G - (H1-(S[1,1]*beta))
        S[6,1] <- 0
      }
    }
    else{
      
      S[4,1] <- G
      S[5,1] <- 0
      S[6,1] <- 0
    }
  }
    
  else if((1-dp2+wp)*(H2-(S[2,1]*beta)) > (1-dp1+wp)*(H1-(S[1,1]*beta)) && (1-dp1+wp)*(H1-(S[1,1]*beta)) >  wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 )){
    if((H2-(S[2,1]*beta)) < G){
      #g1
      S[5,1] <- (H2-(S[2,1]*beta))
      if((H1-(S[1,1]*beta)) < G - (H2-(S[2,1]*beta))){
        #g2
        S[4,1] <- (H1-(S[1,1]*beta))
        #g3
        S[6,1] <- G - (H1-(S[1,1]*beta))- (H2-(S[2,1]*beta))
      }
      else{
        
        S[1,1] <- G - (H2-(S[2,1]*beta))
        S[3,1] <- 0
      }
    }
    else{
      
      S[4,1] <- 0
      S[5,1] <- G
      S[6,1] <- 0
    }
    
  }
  
  else if((1-dp1+wp)*(H1-(S[1,1]*beta)) > wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 )  && wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 ) > (1-dp2+wp)*(H2-(S[2,1]*beta))){
    if((H1-(S[1,1]*beta)) < G){
      
      S[4,1] <- (H1-(S[1,1]*beta))
      S[5,1] <- 0
      S[6,1] <- G - (H1-(S[1,1]*beta))
      
    }
    else{
      
      S[4,1] <- G
      S[5,1] <- 0
      S[6,1] <- 0
    }
  }
  else if(wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 ) > (1-dp1+wp)*(H1-(S[1,1]*beta))  && (1-dp1+wp)*(H1-(S[1,1]*beta)) > (1-dp2+wp)*(H2-(S[2,1]*beta))){
    
    S[4,1] <- 0
    S[5,1] <- 0
    S[6,1] <- G
    
  }
  else if(wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 ) > (1-dp2+wp)*(H2-(S[2,1]*beta))  && (1-dp2+wp)*(H2-(S[2,1]*beta)) > (1-dp1+wp)*(H1-(S[1,1]*beta))){
    
    S[4,1] <- 0
    S[5,1] <- 0
    S[6,1] <- G
    
  }
  else if((1-dp2+wp)*(H2-(S[2,1]*beta)) > wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 )  && wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 ) > (1-dp1+wp)*(H1-(S[1,1]*beta))){
    if((H2-(S[2,1]*beta)) < G){
      
      S[4,1] <- 0
      S[5,1] <- (H2-(S[2,1]*beta))
      S[6,1] <- G - (H2-(S[2,1]*beta))
      
    }
    else{
      
      S[4,1] <- 0
      S[5,1] <- G
      S[6,1] <- 0
    }
  }
  else if((1-dp1+wp)*(H1-(S[1,1]*beta)) == (1-dp2+wp)*(H2-(S[2,1]*beta)) && (1-dp2+wp)*(H2-(S[2,1]*beta)) > wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 )){
    
    if((H1-(S[1,1]*beta)) + (H2-(S[2,1]*beta)) < G){
      #g1
      S[4,1] <- (H1-(S[1,1]*beta))
      S[5,1] <- (H2-(S[2,1]*beta))
      S[6,1] <- G - (H1-(S[1,1]*beta)) - (H2-(S[2,1]*beta))
      }
    else{
      
      S[4,1] <- G/2
      S[5,1] <- G/2
      S[6,1] <- 0
    }
  }
  else if((1-dp2+wp)*(H2-(S[2,1]*beta)) == wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 ) && wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 ) > (1-dp1+wp)*(H1-(S[1,1]*beta))){
    
    if((H2-(S[2,1]*beta)) < G){
      #g1
      S[4,1] <- 0
      S[5,1] <- (H2-(S[2,1]*beta))
      S[6,1] <- G -  (H2-(S[2,1]*beta))
    }
    else{
      
      S[4,1] <- 0
      S[5,1] <- G/2
      S[6,1] <- G/2
    }
  }
  else if((1-dp1+wp)*(H1-(S[1,1]*beta)) == wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 ) && wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 ) > (1-dp2+wp)*(H2-(S[2,1]*beta))){
    
    if((H1-(S[1,1]*beta)) < G){
      #g1
      S[4,1] <- (H1-(S[1,1]*beta))
      S[5,1] <- 0
      S[6,1] <- G -  (H1-(S[1,1]*beta))
    }
    else{
      
      S[4,1] <- G/2
      S[5,1] <- 0
      S[6,1] <- G/2
    }
  }
  else if((1-dp1+wp)*(H1-(S[1,1]*beta)) > (1-dp2+wp)*(H2-(S[2,1]*beta)) && (1-dp2+wp)*(H2-(S[2,1]*beta)) == wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 )){
    
    if((H1-(S[1,1]*beta)) < G){
      #g1
      S[4,1] <- (H1-(S[1,1]*beta))
      if((H2-(S[2,1]*beta)) < G - (H1-(S[1,1]*beta))){
        #g2
        S[5,1] <- (H2-(S[2,1]*beta))
        #g3
        S[6,1] <- G - (H1-(S[1,1]*beta))- (H2-(S[2,1]*beta))
      }
      else{
        
        S[5,1] <- (G - (H1-(S[1,1]*beta)))/2
        S[6,1] <- (G - (H1-(S[1,1]*beta)))/2
      }
    }
    else{
      
      S[4,1] <- G
      S[5,1] <- 0
      S[6,1] <- 0
    }
  }
  else if((1-dp2+wp)*(H2-(S[2,1]*beta)) > (1-dp1+wp)*(H1-(S[1,1]*beta)) && (1-dp1+wp)*(H1-(S[1,1]*beta)) == wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 )){
    
    if((H2-(S[2,1]*beta)) < G){
      #g1
      S[5,1] <- (H2-(S[2,1]*beta))
      if((H1-(S[1,1]*beta)) < G - (H2-(S[2,1]*beta))){
        #g2
        S[4,1] <- (H1-(S[1,1]*beta))
        #g3
        S[6,1] <- G - (H1-(S[1,1]*beta))- (H2-(S[2,1]*beta))
      }
      else{
        
        S[5,1] <- (G - (H2-(S[2,1]*beta)))/2
        S[6,1] <- (G - (H2-(S[2,1]*beta)))/2
      }
    }
    else{
      
      S[4,1] <- 0
      S[5,1] <- G
      S[6,1] <- 0
    }
  }
  else if(wp*(((H1-(S[1,1]*beta))+(H2-(S[2,1]*beta)))/2 ) > (1-dp2+wp)*(H2-(S[2,1]*beta))  && (1-dp2+wp)*(H2-(S[2,1]*beta)) == (1-dp1+wp)*(H1-(S[1,1]*beta))){
    
    S[4,1] <- 0
    S[5,1] <- 0
    S[6,1] <- G
  }
  #Maximum Sugar
  S[7,1] <- s1*(1-ds1+ws)^2*S[1,1]^2*beta + s2*(1-ds2+ws)^2*S[2,1]^2*beta + s3*(ws)^2*S[3,1]^2*beta - ((1-dp1+wp)^(-1/2))*S[4,1] - ((1-dp2+wp)^(-1/2))*S[5,1] - ((wp)^(-1/2))*S[6,1]
  
  S[8,1] <- (1-dp1+wp)*S[4,1]*(H1-(S[1,1]*beta)) + (1-dp2+wp)*S[5,1]*(H2-(S[2,1]*beta)) + (wp)*S[6,1]*((H1+H2)/2)
  write.table(S, "F:/data2.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
}

###########################################################################
Nutrient.adjust <- function(c1, c2, m, n, m1, n1, u, u2){
  
  C <- matrix(data = 0, nrow = 4, ncol = 1, byrow = FALSE, dimnames = NULL) 
  C[1,1] <- u
  C[3,1] <- u2
  
  if(m1 > m){
    C[1,1] <- u+(m1-m)
    if(c1 < C[1,1]){
      C[2,1] <- c1 - c1
      C[1,1] <- C[1,1] - c1
    }
    else{
      C[2,1] <- c1 - C[1,1]
      C[1,1] <- C[1,1] - C[1,1]
    }
  }
  else if(m1 < m){
    if(C[1,1] < m-m1){
      C[2,1] <- c1-((m-m1)-C[1,1])
      C[1,1] <- u - u
    }
    else{
      C[1,1] <- u-(m-m1)
      if(c1 < C[1,1]){
        C[2,1] <- c1 - c1
        C[1,1] <- C[1,1] - c1
      }
      else{
        C[2,1] <- c1 - C[1,1]
        C[1,1] <- C[1,1] - C[1,1]
      }
    }
  }
  else if(m1==m){
    if(m >= 0 && u == 0){
      C[2,1] <- c1
      C[1,1] <- u
    }
    else if(m >= 0 && u > 0){
      if(u <= c1){
        C[2,1] <- c1 - u
        C[1,1] <- u-u
      }
      else{
       C[2,1] <- c1-c1
       C[1,1] <- u - c1
      }
    }
 
  }
 ###################################
  if(n1 > n){
    C[3,1] <- u2+(n1-n)
    if(c2 < C[3,1]){
      C[4,1] <- c2 - c2
      C[3,1] <- C[3,1] - c2
    }
    else{
      C[4,1] <- c2 - C[3,1]
      C[3,1] <- C[3,1] - C[3,1]
    }
  }
  else if(n1 < n){
    if(u2 < n-n1){
      C[4,1] <- c2-((n-n1)-u2)
      C[3,1] <- u2 - u2
    }
    else{
      C[3,1] <- u2-(n-n1)
      if(c2 < C[3,1]){
        C[4,1] <- c2 - c2
        C[3,1] <- C[3,1] - c2
      }
      else{
        C[4,1] <- c2 - C[3,1]
        C[3,1] <- C[3,1] - C[3,1]
      }
    }
  }
  else{
    if(n >= 0 && u2 == 0){
      C[4,1] <- c2
      C[3,1] <- u2
    }
    else if(n >= 0 && u2 > 0){
      if(u2 <= c2){
        C[4,1] <- c2 - u2
        C[3,1] <- u2-u2
      }
      else{
        C[4,1] <- c2-c2
        C[3,1] <- u2 - c2
      }
    }
    
  }
  write.table(C, "F:/data3.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
  
  }
  

