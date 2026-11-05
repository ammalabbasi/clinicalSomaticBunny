# Use a pathology purity estimate to pick among PureCN's ranked solutions, instead of bounding PureCN's search.
#   Rscript purecn_select.R <prefix> <sheet_purity> <window>
# Reads <prefix>.rds (a full-range PureCN.R run). Picks the most likely solution (PureCN's ranking) whose purity is
# within +-window of the sheet purity. Writes <prefix>_purity_check.tsv. If the pick is not rank 1, it writes PureCN's
# curation file <prefix>.csv (Purity/Ploidy of the pick, Curated=TRUE) and exits 3, so the caller regenerates the
# outputs with `PureCN.R --rds <prefix>.rds` (no --force), which reorders the solutions through readCurationFile().
# Exit 0: rank 1 kept (inside the window, or no solution inside it: then the check file says so).
#
# Why not --min-purity/--max-purity: on the synthetic test, bounding the search to 0.40-0.70 around a sheet purity of
# 0.55 (truth 0.60, inside the bounds) moved PureCN's refinement of the correct basin to the 0.70 boundary and it
# reported purity 0.70 / ploidy 5.2 (every copy number doubled), while the full-range search found 0.60 / 2.56.
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
sheet <- as.numeric(args[2])
window <- as.numeric(args[3])
suppressPackageStartupMessages(library(PureCN))
ret <- readRDS(paste0(prefix, ".rds"))
res <- ret$results
pur <- vapply(res, function(x) x$purity, numeric(1))
plo <- vapply(res, function(x) x$ploidy, numeric(1))
ll <- vapply(res, function(x) x$total.log.likelihood, numeric(1))
ok <- which(abs(pur - sheet) <= window + 1e-9)
pick <- if (length(ok)) ok[1] else 1L
action <- if (!length(ok)) "no_solution_within_window_kept_rank1" else if (pick == 1L) "rank1_within_window" else
    paste0("selected_rank", pick)
write.table(data.frame(sheet_purity = sheet, window = window,
        rank1_purity = pur[1], rank1_ploidy = round(plo[1], 3),
        selected_rank = pick, selected_purity = pur[pick], selected_ploidy = round(plo[pick], 3),
        n_solutions = length(res), n_within_window = length(ok), action = action,
        candidates = paste(sprintf("%.2f/%.2f(%.0f)", pur, plo, ll), collapse = ";")),
    file = paste0(prefix, "_purity_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("purity check:", action, sprintf("(sheet %.2f, rank1 %.2f/%.2f, selected %.2f/%.2f)\n",
    sheet, pur[1], plo[1], pur[pick], plo[pick]))
if (pick != 1L) {
    cur <- read.csv(paste0(prefix, ".csv"), as.is = TRUE)
    cur$Purity <- pur[pick]
    cur$Ploidy <- plo[pick]
    cur$Curated <- TRUE
    cur$Comment <- paste(c(if (!is.na(cur$Comment) && nzchar(cur$Comment)) cur$Comment,
        sprintf("SOLUTION %d OF %d SELECTED BY SHEET PURITY %.2f", pick, length(res), sheet)), collapse = ";")
    write.csv(cur, paste0(prefix, ".csv"), row.names = FALSE)
    quit(status = 3)
}
