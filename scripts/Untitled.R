genera_tabella_modalita <- function(lista_var, p) {
  n <- length(lista_var)

  # 1. Invertiamo p: expand.grid vuole prima la variabile più VELOCE.
  # Nel tuo vettore p, la più veloce si trova alla fine, quindi usiamo rev(p).
  ordine_frequenza <- rev(p)

  # 2. Generiamo il prodotto cartesiano con le variabili ordinate per frequenza
  griglia_ordinata <- do.call(expand.grid, lista_var[ordine_frequenza])

  # 3. Ripristiniamo l'ordine originale delle colonne (Colore, Taglia, Marca, Attivo)
  # senza alterare l'ordine delle righe appena generato
  nomi_originali <- names(lista_var)
  df_finale <- griglia_ordinata[, nomi_originali, drop = FALSE]

  # Pulizia estetica dell'output
  rownames(df_finale) <- NULL

  return(df_finale)
}

mie_variabili <- list(
  Colore = c("Rosso", "Verde"),
  Taglia = c("S", "M", "L"),
  Marca  = c("X", "Y"),
  Attivo = c(1, 0)
)

permutazione <- c(2, 1, 3, 4)
matrice_output <- genera_tabella_modalita(mie_variabili, permutazione)
print(matrice_output)
