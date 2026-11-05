# Direct reserve actions under the corrected insurance budgets.
# a: normal-state cost of high reserves; b: crisis-state benefit.
# R0: false-alarm probability; R1: missed-crisis probability.
crra_utility <- function(consumption,sigma) {
 stopifnot(length(sigma)==1L,is.finite(sigma),sigma>0,
           all(is.finite(consumption)),all(consumption>0))
 if(sigma==1)log(consumption)else consumption^(1-sigma)/(1-sigma)
}
action_consumption <- function(rho,th) {
 z<-insurance_parts(th)
 stopifnot(length(rho)>0L,all(is.finite(rho)),all(rho>=0),all(rho<z$A/z$k))
 list(normal=z$A-z$k*rho,crisis=z$B+(1-z$k)*rho)
}
action_components <- function(th,low,high) {
 stopifnot(length(low)==1L,length(high)==1L,low<high)
 c<-action_consumption(c(low,high),th)
 un<-crra_utility(c$normal,th$sigma);us<-crra_utility(c$crisis,th$sigma)
 list(a=un[1]-un[2],b=us[2]-us[1],u_n_low=un[1],u_n_high=un[2],
      u_s_low=us[1],u_s_high=us[2])
}
realized_utility <- function(y,rho,th) {
 stopifnot(length(y)>0L,all(y%in%c(0L,1L)),length(rho)%in%c(1L,length(y)))
 c<-action_consumption(rho,th)
 crra_utility(ifelse(y==1,c$crisis,c$normal),th$sigma)
}
expected_action_utility <- function(p,rho,th) {
 stopifnot(length(p)>0L,all(is.finite(p)),all(p>=0 & p<=1),
           length(rho)%in%c(1L,length(p)) || length(p)==1L)
 c<-action_consumption(rho,th)
 p*crra_utility(c$crisis,th$sigma)+(1-p)*crra_utility(c$normal,th$sigma)
}
ce_gain <- function(U,U_base,sigma) {
 stopifnot(length(U)==1L,length(U_base)==1L,length(sigma)==1L,
           all(is.finite(c(U,U_base,sigma))),sigma>0)
 if(sigma==1)return(expm1(U-U_base))
 stopifnot((1-sigma)*U>0,(1-sigma)*U_base>0)
 expm1(log(U/U_base)/(1-sigma))
}
