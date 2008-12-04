



schmidt <- function (T=25,gas="O2")
{

#Coefficients for the fit of the schmidth number according to:
#Sc = A- BT+CT^@-DT^3, with T in dg C

SchmidtCoeff <- data.frame(
A = c(410.14,855.1,2206.1,1638,1909.1,2205,3412.8,2039.2,2073.1,2301.1,3845.4,3501.8,3531.6,4295.8),
B = c(20.503,46.299,144.86,81.83,125.09,135.71,224.3,120.31,125.62,151.1,228.95,210.31,231.4,281.52),
C=c(0.53175,1.254,4.5413,1.483,3.9012,3.9549,6.7954,3.4209,3.6276,4.7364,6.1908,6.1851,7.2168,8.7826),
D=c(0.006011,0.01449,0.056988,0.008004,0.048953,0.047339,0.083,0.040437,0.043219,0.059431,0.06743,0.07513,0.090558,0.11025))
rownames(SchmidtCoeff) <- c("He", "Ne", "N2", "O2", "Ar",
        "Kr", "Rn", "CH4","CO2", "N2O", "CCl2F2", "CCL3F",
        "SF6", "CCl4")

Sc<-SchmidtCoeff[gas,]

schmidt <-Sc$A-Sc$B*T+Sc$C*T*T-Sc$D*T*T*T
schmidt
}

gastransfer<- function (T=25,u10=1,gas="O2",method=c("Liss","Nightingale","Wanninkhof1","Wanninkhof2"),
           Schmidt=schmidt(T=T,gas=gas))
{
 method <- match.arg(method)
 S600<-Schmidt/600
 tr <- switch(method, Liss =
 ifelse (u10<3.6,0.17*u10*Schmidt^(-2/3),ifelse (u10<=13, (u10-3.4)*2.8/sqrt(S600),
 (u10-8.4)*5.9/sqrt(S600))), Nightingale = (0.33*u10+0.222*u10*u10)/sqrt(S600),
  Wanninkhof1=0.31*u10^2*sqrt(S600),Wanninkhof2=0.0283*u10^3/sqrt(S600))
 attr(tr, "unit") = "cm/hr"
 tr
}
