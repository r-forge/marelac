## =============================================================================
## From a table in "database format" (3-columns) to a "crosstab"
## =============================================================================

                                               
db2cross <- function(input, row = 1, col = 2, value = 3) {

  if (is.character(value))
    value <- which(colnames(input) == value)
  if (is.character(row))
    row <- which(colnames(input) == row)
  if (is.character(col))
    col <- which(colnames(input) == col)

  # Check if input has only 3 columns
  dim.input <- dim(IN <- input[,c(col,row,value)])
  if (ncol(IN) != 3)
    stop ("'input' should have three columns")

  cols <- sort(unique(IN[,1]))
  rows <- sort(unique(IN[,2]))

  nr = length(rows)
  nc = length(cols)
  out <- .Fortran("crosstab",t(IN),as.integer(nrow(input)),
             as.integer(1),as.integer(2), as.integer(3),
             as.double(cols), as.double(rows), nr = nr, nc = nc, 
             cross = matrix(nr = nr, nc = nc, data = as.double(0.)),
             package = "kaso")  
  list(x=rows, y= cols, z = out$cross)             
}


