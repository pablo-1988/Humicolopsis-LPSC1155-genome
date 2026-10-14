# A7.3 — Árbol ultramétrico para CAFE5 a partir del árbol particionado de A1 (Fig. 1c).
# Enraiza en N. crassa, quita CBS 645.76 (sin genoma) y calibra con chronos (edad de la raíz = 1, tiempo relativo).
# Modelo "relaxed", lambda = 0.1: entre 12 combinaciones probadas (correlated/relaxed/discrete/clock x lambda 0.1/1/10)
# es la de menor PHIIC que converge y no deja ramas de longitud ~0.
# Uso: Rscript 04_arbol_ultrametrico.R <arbol_A1> <salida.nwk>
suppressMessages(library(ape))
args <- commandArgs(trailingOnly = TRUE)
tr <- read.tree(args[1])
tr <- drop.tip(tr, "HumcCBS645")
tr <- root(tr, "Neurcras", resolve.root = TRUE)
tr$node.label <- NULL
set.seed(1)
cal <- makeChronosCalib(tr, node = "root", age.min = 1, age.max = 1)
ut <- withCallingHandlers(chronos(tr, lambda = 0.1, model = "relaxed", calibration = cal, quiet = TRUE),
                 warning = function(w) stop("chronos no convergió: ", conditionMessage(w)))
class(ut) <- "phylo"
stopifnot(is.ultrametric(ut, tol = 1e-6), is.binary(ut), min(ut$edge.length) > 1e-3)
write.tree(ut, args[2], digits = 6)
cat("Taxa:", Ntip(ut), "\n")
