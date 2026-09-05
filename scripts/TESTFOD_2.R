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
metrics_r = c("Dominance")

# TOTAL BINS: Obbligatorio per alcune metriche, passiamo 1 di default
total_bins_r = as.integer(0)

# PARAMETRI OPZIONALI: In R corrispondono a NULL
count_r = 1000000
seed_r = NULL
output_interval_r = 2


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
