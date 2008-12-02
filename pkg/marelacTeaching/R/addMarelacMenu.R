remMarelacMenu <- function() {
  if ("marelac" %in% winMenuNames()) winMenuDel("marelac")
}

addMarelacMenu <- function() {
  winMenuAdd("marelac")
  winMenuAddItem("marelac", "Text Help", "sethelp(1)")
  winMenuAddItem("marelac", "HTML Help", "sethelp(2)")
  winMenuAddItem("marelac", "Windows Help", "sethelp(3)")
  winMenuAddItem("marelac", "remove the marelac Menue", "remSpecialMenu()")
}



