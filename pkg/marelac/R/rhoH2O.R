`rhoH2O` <-
function(T, p=0, S=0, method=c("Chen", "Fofonoff")) {
   method <- match.arg(method)
   switch(method,
     Chen     = rhoH2O_Chen(T, p, S),
     Fofonoff = 1e-3 / sv_Fofonoff(T, p, S)
   )
}

