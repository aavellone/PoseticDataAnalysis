## Golden reference generator: computes results in R and writes them to
## reference/ for the Python side (check.py) to compare against.
## Run from the _golden/ directory (see run_compare.sh).
suppressMessages(library(poseticDataAnalysis))

dir.create("reference", showWarnings = FALSE)
wf <- function(name) file.path("reference", name)
wm <- function(m, name) write.table(m, wf(name), row.names = FALSE, col.names = FALSE, quote = FALSE)

## --- Shared poset: a,b <= c <= e ; d <= e ---------------------------------
elems <- c("a", "b", "c", "d", "e")
dom <- rbind(c("a", "c"), c("b", "c"), c("c", "e"), c("d", "e"))
pos <- POSet(elements = elems, dom = dom)

writeLines(POSetElements(pos), wf("elements.txt"))
wm(DominanceMatrix(pos), "incidence.txt")
ord <- OrderRelation(pos)
writeLines(sort(apply(ord, 1, function(r) paste(r[1], r[2], sep = "<="))), wf("order.txt"))
writeLines(sort(POSetMaximals(pos)), wf("maximal.txt"))
writeLines(sort(POSetMinimals(pos)), wf("minimal.txt"))

## Exact MRP (deterministic)
emrp <- ExactMRP(pos)
wm(emrp$mrp, "exact_mrp.txt")
writeLines(as.character(emrp$n), wf("exact_mrp_n.txt"))

## MCMC (fixed seed -> bit-identical with Python)
gle <- LEBubleyDyer(pos, seed = "123456789")
wm(LEGet(gle, from_start = TRUE, n = 25), "le_samples.txt")
gmrp <- BubleyDyerMRPGenerator(pos, seed = "987654321")
wm(BubleyDyerMRP(gmrp, n = 200000)$mrp, "bd_mrp.txt")

## Lex / Fuzzy / FOD / dim-reduction
wm(LexMRP(2, 3), "lexmrp.txt")
wm(LexSeparation(2, 3, types = "symmetric")$symmetric, "lexsep_sym.txt")

fdom <- matrix(c(1.0, 0.2, 0.7, 0.8, 1.0, 0.4, 0.3, 0.6, 1.0), nrow = 3, byrow = TRUE)
dimnames(fdom) <- list(c("a", "b", "c"), c("a", "b", "c"))
fs <- FuzzySeparationMinMax(fdom, types = c("symmetric", "asymmetricLower"))
wm(fs$symmetric, "fuzzysep_sym.txt")
wm(fs$asymmetricLower, "fuzzysep_lower.txt")

fpos <- LinearPOSet(c("1", "2", "3"))
freq <- matrix(c(10, 0, 5, 5, 0, 10), nrow = 3, byrow = TRUE)
dimnames(freq) <- list(c("1", "2", "3"), c("G1", "G2"))
fod <- FirstOrderDominanceAnalysis(fpos, freq, metrics = c("Dominance"))
wm(fod$Dominance$mintr.delta, "fod_closed.txt")
wm(fod$Dominance$bin, "fod_bin.txt")

profile <- matrix(c(0,0,0, 1,0,0, 1,1,0, 1,1,1), nrow = 4, byrow = TRUE)
dr <- OptimalBidimensionalEmbedding(profile, c(1, 2, 3, 4), "absolute")
writeLines(format(dr$bestLossValue, digits = 17), wf("dimred_bestloss.txt"))

cat("Golden reference written to reference/\n")
