"sethelp" <-
function(opt=1) {
    set <- c(FALSE, FALSE, FALSE)
    #if (opt %in% 1:3) 
    set[opt] <- TRUE
    options(winhelp=set[1])
    options(htmlhelp=set[2])
    options(chmhelp=set[3])
}

