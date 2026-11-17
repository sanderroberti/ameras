transform_contract_data <- function() {
  i <- seq_len(60)
  dose <- seq(-1, 1, length.out = length(i))
  data.frame(Y = 1 + 0.5 * dose + 0.2 * sin(i * 1.7), D1 = dose, D2 = dose + 0.01)
}
