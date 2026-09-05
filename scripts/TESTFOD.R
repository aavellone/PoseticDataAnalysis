rm(list = ls())
library(poseticDataAnalysis)

print("=> Preparazione dei parametri (SEXP) in R...")
# TOLLERANZA
tolerance_r = 1e-9;
sep_r = '_'



el1 = c('e1', 'e5', 'e3', 'e4', 'e2')
dom1 = matrix(c('e1', 'e2',
                'e1', 'e3',
                'e1', 'e4',
                'e2', 'e5'),
              ncol=2,
              byrow=TRUE)
poset1_r = POSet(elements = el1, dom = dom1)
freq_matrix_r = matrix(c(
  0.2, 0.2, 0.4,
  0.2, 0.2, 0.2,
  0.6, 0.6, 0.4),
  ncol = 3,
  byrow = TRUE,
  dimnames = list(c('e1', 'e5', 'e4'),
                  c('GruppoA', 'GruppoB', 'GruppoC')
  )
)

#subpopulation_count_r = c(100, 200, 300);
subpopulation_count_r = NULL

# Dominance, MannWhitneyDominance, MannWhitneyInferentialDominance
#metrics_r = c("Dominance", "MannWhitneyDominance", 'MannWhitneyInferentialDominance')
metrics_r = c("Dominance", "MannWhitneyDominance")

# TOTAL BINS: Obbligatorio per alcune metriche, passiamo 1 di default
total_bins_r = as.integer(0)

# PARAMETRI OPZIONALI: In R corrispondono a NULL
count_r = 10
seed_r = NULL
output_interval_r = NULL


risultato_fod <- FirstOrderDominanceAnalysis(posets = poset1_r,
                                        freq_matrix = freq_matrix_r,
                                        metrics = metrics_r,
                                        subpopulation_count = subpopulation_count_r,
                                        total_bins = total_bins_r,
                                        count = count_r,
                                        seed = seed_r,
                                        output_interval_in_sec = output_interval_r,
                                        sep = sep_r,
                                        tolerance = tolerance_r)


#############################

el2 = c(
  '0000', '0001', '0010', '0011',
  '0100', '0101', '0110', '0111',
  '1000', '1001', '1010', '1011',
  '1100', '1101', '1110', '1111')



dom2 = matrix(c('0000', '0001',
                '0000', '0010',
                '0000', '0100',
                '0000', '1000',

                '0001', '0011',
                '0001', '0101',
                '0001', '1001',

                '0010', '0011',
                '0010', '0110',
                '0010', '1010',

                '0100', '0101',
                '0100', '0110',
                '0100', '1100',

                '0011', '0111',
                '0011', '1011',

                '0101', '0111',
                '0101', '1101',

                '0110', '0111',
                '0110', '1110',

                '1001', '1011',
                '1001', '1101',

                '1010', '1011',
                '1010', '1110',

                '1100', '1101',
                '1100', '1110',

                '1110', '1111',
                '1101', '1111',
                '1011', '1111',
                '0111', '1111'

                ),
              ncol=2,
              byrow=TRUE)
poset1_r_2 = POSet(elements = el2, dom = dom2)
freq_matrix_r_2 = matrix(c(
  0.0912,0.0254,0.0189,0.022,0.0096,
  0.0415,0.0369,0.0512,0.0238,0.0065,
  0.0092,0.0016,0.0038,0.0062,0.0011,
  0.093,0.081,0.0814,0.0669,0.0528,
  0.0384,0.0248,0.0224,0.0159,0.0032,
  0.079,0.0553,0.0452,0.0407,0.0091,
  0.0217,0.0052,0.0049,0.0045,0,
  0.2337,0.2294,0.1995,0.1709,0.1232,
  0.0144,0.0164,0.01,0.006,0.0103,
  0.0405,0.0337,0.0603,0.0517,0.0476,
  0.0022,0.0014,0.0039,0.0028,0.0023,
  0.0679,0.1162,0.1358,0.1826,0.2067,
  0.0183,0.0094,0,0.0041,0.0095,
  0.0555,0.0554,0.0408,0.0485,0.0638,
  0.0038,0.0045,0.004,0.0032,0.0074,
  0.1896,0.3034,0.318,0.3501,0.4468
  ),
  ncol = 5,
  byrow = TRUE,
  dimnames = list(c('0000', '1000', '0100', '1100', '0010', '1010', '0110', '1110', '0001', '1001', '0101', '1101', '0011', '1011', '0111', '1111'),
                  c('Basic', 'Vocat.', 'Short', 'Medium', 'Long')
  )
)

subpopulation_count_r_2 = NULL



metrics_r_2 = c("Dominance", "MannWhitneyDominance")

# TOTAL BINS: Obbligatorio per alcune metriche, passiamo 1 di default
total_bins_r_2 = as.integer(0)

# PARAMETRI OPZIONALI: In R corrispondono a NULL
count_r_2 = NULL
seed_r_2 = NULL
output_interval_r_2 = NULL
tolerance_r_2 = 1e-9;
sep_r_2 = '_'


risultato_fod_2 <- FirstOrderDominanceAnalysis(posets = poset1_r_2,
                                             freq_matrix = freq_matrix_r_2,
                                             metrics = metrics_r_2,
                                             subpopulation_count = subpopulation_count_r,
                                             total_bins = total_bins_r_2,
                                             count = count_r_2,
                                             seed = seed_r_2,
                                             output_interval_in_sec = output_interval_r_2,
                                             sep = sep_r_2,
                                             tolerance = tolerance_r_2)




if(FALSE) {
# Inizializziamo le liste finali per i POSet e per memorizzare eventuali errori rilevati
final_posets_list <- list()
problemi_rilevati <- list()

for (metric_name in names(risultato_fod_2)) {

  # Estraiamo l'oggetto posets (matrici Z) e i quasi-ordini dall'output dell'analisi
  z_matrices_list   <- risultato_fod_2[[metric_name]][["Zs"]]
  quasi_orders_list <- risultato_fod_2[[metric_name]][["quasi.orders"]]
  food_closed <- risultato_fod_2[[metric_name]][["mintr.delta"]]

  # Controlliamo per sicurezza che la metrica contenga dati validi
  if (!is.null(z_matrices_list) && is.list(z_matrices_list)) {

    alpha_names <- names(z_matrices_list)

    # Utilizziamo lapply per scorrere ogni livello di alpha rilevato per la metrica corrente
    posets_for_alpha <- lapply(seq_along(z_matrices_list), function(idx) {

      z_matrix  <- z_matrices_list[[idx]]
      alpha_val <- alpha_names[idx]

      # --- CONTROLLO DI SICUREZZA DI ANTISIMMETRIA SULLA MATRICE Z ---
      M_check <- z_matrix
      diag(M_check) <- 0 # Azzera la diagonale (la riflessività non viola l'antisimmetria)

      # L'operatore '&' tra la matrice e la sua trasposta individua le celle simmetriche (cicli)
      if (any(M_check & t(M_check))) {

        # Estraiamo il quasi-ordine corrispondente a questo specifico alpha per il confronto
        quasi_matrix <- quasi_orders_list[[alpha_val]]

        # Tracciamo il problema per notificarlo all'utente
        problemi_rilevati[[paste(metric_name, alpha_val, sep = " | alpha: ")]] <<- TRUE

        cat("\n=======================================================================\n")
        cat("⚠️ [PROBLEMA INDIVIDUATO]: VIOLAZIONE DELL'ANTISIMMETRIA (RILEVATI CICLI)\n")
        cat("METRICA: ", metric_name, "\n")
        cat("VALORE SOGLIA ALPHA: ", alpha_val, "\n")
        cat("=======================================================================\n")

        cat("\n-----------------------------------------------------------------------\n")
        cat("[2] MATRICE FOOD.CLOSED:\n")
        print(food_closed)
        cat("=======================================================================\n\n")




        if (!is.null(quasi_matrix)) {
          cat("\n[1] MATRICE DEL QUASI-ORDINE GREZZO (Prima del collasso, 0/1):\n")
          # Convertiamo momentaneamente in formato numerico binario per facilitare la lettura visiva in console
          print(apply(quasi_matrix, c(1, 2), function(x) as.integer(x)))

          # Individuiamo esattamente quali elementi si stanno chiudendo in ciclo nel quasi-ordine
          Q_no_diag <- quasi_matrix
          diag(Q_no_diag) <- FALSE
          cicli_quasi <- which(Q_no_diag & t(Q_no_diag), arr.ind = TRUE)

          if (nrow(cicli_quasi) > 0) {
            cat("\n--> Elementi che si dominano reciprocamente nel Quasi-Ordine:\n")
            for (i in 1:nrow(cicli_quasi)) {
              if (cicli_quasi[i, 1] < cicli_quasi[i, 2]) { # Evita di stampare due volte la stessa coppia invertita
                cat(sprintf("    - '%s' <=> '%s'\n",
                            rownames(quasi_matrix)[cicli_quasi[i, 1]],
                            colnames(quasi_matrix)[cicli_quasi[i, 2]]))
              }
            }
          }
        } else {
          cat("\n[!] Matrice del quasi-ordine non trovata nell'oggetto. Assicurati di averla abilitata nell'output.\n")
        }

        cat("\n-----------------------------------------------------------------------\n")
        cat("[2] MATRICE Z QUOZIENTE RISULTANTE (Incompleta / Non Antisimmetrica):\n")
        print(z_matrix)
        cat("=======================================================================\n\n")
      }

      # Se il controllo fallisce, restituiamo un valore nullo o proviamo comunque a generare
      # il POSet (nota: POSetFromIncidence fallirà a sua volta se internamente valida l'antisimmetria)
      tryCatch({
        return(poseticDataAnalysis::POSetFromIncidence(z_matrix))
      }, error = function(e) {
        # Restituisce NULL per questo specifico livello di alpha se l'oggetto S4 non può essere istanziato
        return(NULL)
      })
    })

    # Preserviamo i nomi degli alpha associati
    names(posets_for_alpha) <- alpha_names

    # Rimuoviamo eventuali elementi NULL (quelli che hanno fallito la creazione) se desideri una lista pulita
    # posets_for_alpha <- posets_for_alpha[!sapply(posets_for_alpha, is.null)]

    # Salviamo il risultato per la metrica corrente
    final_posets_list[[metric_name]] <- posets_for_alpha

  } else {
    final_posets_list[[metric_name]] <- list()
  }
}

# --- NOTIFICA FINALE ---
if (length(problemi_rilevati) > 0) {
  stop(sprintf("L'elaborazione è terminata. Sono stati riscontrati e visualizzati problemi di antisimmetria in %d configurazione/i (controlla i log sopra).", length(problemi_rilevati)), call. = FALSE)
} else {
  cat("\n✅ Tutti i controlli completati con successo: nessuna violazione di antisimmetria riscontrata.\n")
}





# Supopniamo che 'fod_output' sia l'oggetto restituito da FirstOrderDominanceAnalysis
# fod_output <- FirstOrderDominanceAnalysis(...)

# Inizializziamo una lista per contenere la nuova struttura con gli oggetti S4 POSet
final_posets_list <- list()

# 1. Iteriamo su tutte le metriche calcolate (es. "Dominance", "MannWhitneyDominance")
for (metric_name in names(risultato_fod_2)) {

  # Estraiamo l'oggetto G (che contiene la lista di matrici Z per ogni alpha)
  z_matrices_list <- risultato_fod_2[[metric_name]][["posets"]]

  # Controlliamo per sicurezza che la metrica contenga dati validi
  if (!is.null(z_matrices_list) && is.list(z_matrices_list)) {

    # 2. Utilizziamo lapply per scorrere ogni matrice di incidenza Z (associata a un alpha)
    # e convertirla in un vero oggetto S4 di classe "POSet"
    posets_for_alpha <- lapply(z_matrices_list, function(z_matrix) {

      # Chiamata esplicita alla funzione di conversione del pacchetto

      M_check <- z_matrix
      diag(M_check) <- 0 # o FALSE se logica

      # L'operatore '&' tra la matrice e la sua trasposta trova le celle simmetriche.
      # Se 'any' è TRUE, significa che per qualche i != j hai sia i <= j che j <= i.
      if (any(M_check & t(M_check))) {
        print(metric_name)
        print(z_matrix)
        stop("The incidence matrix violates the antisymmetry property (cycles detected).", call. = FALSE)

      }


      #poseticDataAnalysis::POSetFromIncidence(z_matrix)

    })

    # Preserviamo i nomi della lista originale (le stringhe che rappresentano i valori di alpha)
    names(posets_for_alpha) <- names(z_matrices_list)

    # Salviamo il risultato per la metrica corrente
    final_posets_list[[metric_name]] <- posets_for_alpha

  } else {
    final_posets_list[[metric_name]] <- list()
  }
}

# Al termine del ciclo, final_posets_list sarà una struttura speculare a quella originale,
# ma final_posets_list[[metric_name]][[as.character(alpha)]] conterrà l'oggetto S4 POSet.

}
