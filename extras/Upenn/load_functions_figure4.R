plotLocalLik <- function(plot.data){
  ggplot(plot.data, aes(x = hr, y = ll))+
    geom_line(size = 1,color = "grey46")+
    scale_x_continuous(expand = c(0, 0),position = "bottom",trans = "log2", breaks = c(0.1, 0.5,2,4))+
    scale_y_continuous(expand = c(0, 0),limits = c(-8,0), breaks = c(0,-3,-6))+
    facet_grid(site~., switch = "y")+
    theme_bw()+
    ggplot2::theme(legend.position = "top",
                   legend.title = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank(),
                   strip.text.y.left = element_text(angle = 0),
                   panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.major.x = ggplot2::element_blank(),
                   #panel.grid.minor.y = ggplot2::element_blank(), 
                   #panel.grid.minor.x = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   #axis.title.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA, size = 0.8)
    )+
    theme(text = element_text(size = 18),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 17)
    )+
    #theme(strip.placement = "outside",
    #      strip.background = ggplot2::element_blank())+
    labs(x = "Hazard Ratio", y="Relative log-likelihood")
}

plotGlobalLik.comparison <- function(plot.data){
  cols <- c("gray46","#FF9999","#00BFC4")
  ggplot(plot.data, aes(x = x, y = ll, color = method))+
    #geom_hline(yintercept  = -1.92, size = 0.3, alpha = 0.7)+
    geom_line(aes(linetype = method), size = 1)+
    scale_colour_manual(values=cols)+
    scale_x_continuous(expand = c(0, 0),position = "bottom",trans = "log2", breaks = c(0.1, 0.5,2,4))+
    scale_y_continuous(expand = c(0, 0),limits = c(-8,0), breaks = c(0,-3,-6))+
    theme_bw()+
    ggplot2::theme(legend.position = "top",
                   legend.title = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank(),
                   #panel.grid.major.y = ggplot2::element_blank(), 
                   #panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(), 
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   #axis.title.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA, size = 0.8)
    )+
    theme(text = element_text(size = 18),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 18)
    )+
    #theme(strip.placement = "outside",
    #      strip.background = ggplot2::element_blank())+
    labs(x = "", y="")
}

## cumsum from pda
rcpp_aggregate <- function(x, indices, simplify = TRUE, cumulative = FALSE, reversely = FALSE) {
  .Call('_pda_rcpp_aggregate', PACKAGE = 'pda', x, indices, simplify, cumulative, reversely)
}

## functions for local log-likelihoods and derivatives
logL_single_deriv<-function(bbar,ipdata){
  if(sum(ipdata$status==TRUE)==0){
    logL=0;logL_D1=0;logL_D2=0;logL_D3=0;logL_D4=0
  }else{
    
    T_all<-sort(unique(ipdata$time[ipdata$status==TRUE]))
    nt<-length(T_all)
    t_max <- max(ipdata$time)+1
    #ipdata0 <- ipdata
    ipdata <- ipdata[!ipdata$time<min(T_all),]
    
    #generate dataframe in format expected by ODAC
    pfdata <- cbind(T_all, 0, rep(0))
    colnames(pfdata)<-colnames(ipdata)
    pfdata <- rbind(ipdata, pfdata)
    pfdata <- pfdata[order(pfdata$time),]
    pfdata$interval <- cut(pfdata$time, breaks = c(T_all, t_max), labels = 1:nt, right=FALSE)
    pfdata <- pfdata[order(pfdata$interval),]
    pfdata$interval[is.na(pfdata$interval)]<-nt
    X <- c(pfdata$x)
    # summary stats: U, W, Z
    eXb <- c(exp(X*bbar))
    UWZ <- eXb * cbind(1, X, X^2, X^3, X^4)
    # rcpp_aggregate() is a function written in rcpp for calculating column-wise (reverse) cumsum
    # credit to Dr Wenjie Wang
    UWZ <- rcpp_aggregate(x = UWZ, indices = pfdata$interval, cumulative = T, reversely = T)
    
    # since fake X=0, cumulative W and Z will be the same, 
    # but exp(Xb)=1, so need to remove cumulated ones from each time pts
    U <- UWZ[,1] - c(nt:1)
    W <- UWZ[,2]
    Z <- UWZ[,3]
    Z3 <- UWZ[,4]
    Z4 <- UWZ[,5]
    #d <- c(table(ipdata[ipdata$status==TRUE,c(1,2)]))
    ipdataE <- ipdata[ipdata$status==TRUE,]
    d<-c(aggregate(status~time,ipdataE, FUN=sum)$status)
    #X <- as.matrix(ipdata[ipdata$status==TRUE, -c(1,2)])
    X<-c(aggregate(x~time,ipdataE, FUN=mean)$x) #breslow's method for ties
    eXb <- c(exp(X*bbar))
    logL_D1 <- sum(d*X) - sum(d * W / U,na.rm = TRUE)
    logL_D2 <- sum(d * (W^2 - U*Z) / U^2)
    logL_D3 <- -sum(d * (Z3/U - 3*Z*W / (U^2) + 2*W^3/(U^3)))
    logL_D4 <- sum(d* (-Z4/U+4*Z3*W/(U^2)+3*(Z/U)^2-12*Z*W^2/(U^3)+6*W^4/(U^4)))
    logL <- sum(d*log(eXb/U)) 
  }
  derivatives <- list(logL=logL,logL_D1=logL_D1,logL_D2=logL_D2,logL_D3=logL_D3,logL_D4=logL_D4)
  return(derivatives)
}

GetLocalDeriv <- function(InputData,InputBeta){  # get pade approximation for each site
  InputData <- cbind.data.frame(InputData$time,InputData$y,InputData$x,InputData$stratumId)
  colnames(InputData) <- c("time","status","x","stratumId")
  nStrata <- max(InputData$stratumId)
  logL.deriv <- rep(0)
  for (StrataID in 1:nStrata){
    ipdata <- InputData[InputData$stratumId==StrataID,c("time","status","x")]
    logL.deriv<-logL.deriv+unlist(logL_single_deriv(InputBeta,ipdata))
  }
  return(logL.deriv)
}
GetGlobalPadeCoef<- function(InputBeta,populations){
  Derivs<-lapply(populations, GetLocalDeriv,InputBeta=InputBeta)
  Derivs<-do.call("rbind", Derivs)
  GlobalDeriv <- colSums(Derivs)
  c0 <- GlobalDeriv[1]
  c1 <- GlobalDeriv[2]
  c2 <- GlobalDeriv[3]/2
  c3 <- GlobalDeriv[4]/6
  c4 <- GlobalDeriv[5]/24
  
  b2 <- (c3^2-c2*c4)/(c2^2-c1*c3)
  b1 <- -c3/c2-c1/c2*b2
  a0 <- c0
  a1 <- c1+c0*b1
  a2 <- c2+c1*b1+c0*b2
  P2_coef <- c(a0,a1,a2)
  Q2_coef <- c(1,b1,b2)
  return(list(num_coef = P2_coef,denom_coef=Q2_coef))
}
PadeEstCI<- function(InputBeta,PadeCoef){
  num_vec<-c(PadeCoef$num_coef)
  denom_vec<-c(PadeCoef$denom_coef)
  Pade<-function(bbar,beta,PadeCoef){  # bbar is the initial point, beta is the point you'd like to evaluate
    x <- beta-bbar
    Pm<- cbind(1,x,x^2)%*%PadeCoef$num_coef
    Qn <- cbind(1,x,x^2)%*%PadeCoef$denom_coef
    pade <- Pm/Qn
    return(pade)
  }
  # Closed form for estimates
  A<- num_vec[3]*denom_vec[2]-num_vec[2]*denom_vec[3]
  B<- 2*(num_vec[3]-num_vec[1]*denom_vec[3])
  C<- num_vec[2]-num_vec[1]*denom_vec[2]
  est1<-InputBeta+(-B-sqrt(B^2-4*A*C))/2/A
  #est2<-InputBeta+(-B+sqrt(B^2-4*A*C))/2/A
  #PadeEst<- ifelse(Pade(InputBeta,est1,PadeCoef)>Pade(InputBeta,est2,PadeCoef),est1,est2)
  PadeEst<-est1
  # Closed form for CI
  coef<-num_vec-denom_vec*c(Pade(InputBeta,PadeEst,PadeCoef)-1.92)
  PadeRoots<- c((-coef[2]+sqrt(coef[2]^2-4*coef[1]*coef[3]))/(2*coef[3]),(-coef[2]-sqrt(coef[2]^2-4*coef[1]*coef[3]))/(2*coef[3]))
  #PadeCI<-InputBeta+c(min(PadeRoots),max(PadeRoots))
  PadeCI<-InputBeta+PadeRoots
  return(list(Est=PadeEst,CI=PadeCI))
}

