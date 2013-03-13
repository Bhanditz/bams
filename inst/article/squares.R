squares <- function(algos,square.inches="0.08"){
  sprintf("%s \\textcolor{%s.color}{\\rule{%sin}{%sin}}",
          algos,
          algos,
          square.inches,
          square.inches)
}
