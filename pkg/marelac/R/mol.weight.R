
# converts from mol to g
mol.weight2 <- function(species) {
  with(AtomicWeight, {
    trace = FALSE # for debugging, lines have to be removed in production versions
    # insert * before number (with one or more digits)
    s1 <- gsub("([0-9]+)", "\\*\\1\\+", species)
    if (trace) cat(s1, "\n")
    # insert + after capital letters
    s1 <- gsub("([a-z,A-Z])", "\\1\\+", s1, perl=TRUE)
    if (trace) cat(s1, "\n")
    # remove + before lower case letters
    s1 <- gsub("\\+([a-z])", "\\1", s1, perl=TRUE)
    if (trace) cat(s1, "\n")
    # replace +* with only *
    s1 <- gsub("\\+\\*", "\\*", s1, perl=TRUE)
    if (trace) cat(s1, "\n")
    # remove trailing +
    s1 <- gsub("\\+$", "", s1, perl=TRUE)
    if (trace) cat(s1, "\n")
    # remove trailing +)
    s1 <- gsub("\\+\\)", "\\)", s1, perl=TRUE)
    if (trace) cat(s1, "\n")
    # calculate molar mass
    eval(parse(text=s1))
  })
}





