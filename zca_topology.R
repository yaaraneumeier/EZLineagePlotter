# ============================================================================
# zca_topology.R - does each ZCA cluster correspond to the tree topology?
#
# Tree: trees2025_ordering/*_bootstrapped_rerooted.newick (the topology the
# figures draw; rotations do not change it). Clusters: masterlist v84 ZCA,
# identical to the slide map zca_leaf_cluster_map_20260909.csv.
#
# Per cluster (and per parent group of lettered clusters, e.g. 1A+1B -> "1"):
#   clade_drawn     leaves == descendants of one node, rooted as drawn
#                   (the newick's basal node, as ggtree draws it)
#   split_unrooted  leaves or their complement form one side of an edge
#   n_pieces        number of maximal pure clades the cluster splits into
#   mrca_*          size / support of the smallest drawn clade containing it,
#                   and the other-cluster leaves inside it ("intruders")
#
# Output: out/ZCA_trees/zca_topology_clusters.tsv, zca_topology_summary.tsv
# ============================================================================
suppressPackageStartupMessages(library(ape))
dir.create("out/ZCA_trees", showWarnings = FALSE, recursive = TRUE)
KLEIN <- normalizePath(Sys.getenv("EZ_DATA", file.path(Sys.getenv("EZ_ROOT", "."), "data")), mustWork = TRUE)   # run from the repo root
TREES <- c("775-13" = "775_13", "795-09" = "795_09", "841-12" = "841_12", "MM15-127" = "mm15_127",
           "MM16-412" = "mm16_412", "MM16-423" = "mm16_423", "BC15-0267" = "bc15_0267",
           "BC16-0401" = "bc16_0401", "BC16-0545" = "bc16_0545", "BC14-2680" = "bc14_2680")

ml <- read.csv(file.path(KLEIN, "masterlist_lineage_tree_v84_20260923.csv"), colClasses = "character",
               check.names = FALSE, fileEncoding = "latin1")[, c("new_chosen_sr", "Individual Name", "ZCA")]
slides <- read.csv(file.path(KLEIN, "zca_leaf_cluster_map_20260909.csv"), colClasses = "character")

rows <- list()
for (ind in names(TREES)) {
  tr <- read.tree(file.path(KLEIN, "trees2025_ordering", paste0(TREES[[ind]], "_bootstrapped_rerooted.newick")))
  n <- Ntip(tr); tips <- tr$tip.label
  zca <- setNames(ml$ZCA, ml$new_chosen_sr)[tips]
  stopifnot(!anyNA(zca))
  sl <- slides[slides$individual == ind, ]
  stopifnot(setequal(sl$leaf_id, tips), all(setNames(sl$cluster, sl$leaf_id)[tips] == zca))

  # descendant tip sets for every node (tips: themselves), rooted as drawn
  desc <- c(as.list(seq_len(n)), prop.part(tr))            # index = node number
  parent <- integer(n + tr$Nnode); parent[tr$edge[, 2]] <- tr$edge[, 1]
  support <- c(rep(NA, n), suppressWarnings(as.numeric(tr$node.label)))
  sizes <- lengths(desc)

  groups <- split(tips, zca)
  lettered <- grepl("^Cluster [0-9]+[A-Z]$", names(groups))
  parents <- unique(sub("[A-Z]$", "", names(groups)[lettered]))
  for (pg in parents) {
    kids <- names(groups)[startsWith(names(groups), pg) & lettered]
    if (length(kids) > 1 && !pg %in% names(groups)) groups[[paste0(pg, " (", paste(sub(".* ", "", kids), collapse = "+"), ")")]] <- unlist(groups[kids])
  }

  for (g in names(groups)) {
    idx <- match(groups[[g]], tips); k <- length(idx)
    comp <- setdiff(seq_len(n), idx)
    in_g <- vapply(desc, function(d) all(d %in% idx), logical(1))
    exact <- which(in_g & sizes == k)
    clade_drawn <- length(exact) == 1
    split_unrooted <- clade_drawn || k == n ||
      any(vapply(desc, function(d) length(d) == length(comp) && all(d %in% comp), logical(1)))
    # maximal pure clades
    parent_in_g <- ifelse(parent == 0, FALSE, in_g[pmax(parent, 1)])
    pure_max <- which(in_g & !parent_in_g)
    # smallest drawn clade containing the group
    contains <- which(vapply(desc, function(d) all(idx %in% d), logical(1)))
    mrca <- contains[which.min(sizes[contains])]
    intr <- setdiff(desc[[mrca]], idx)
    intr_tab <- if (length(intr)) paste(sprintf("%s:%d", names(table(zca[intr])), as.integer(table(zca[intr]))), collapse = ", ") else ""
    rows[[length(rows) + 1]] <- data.frame(
      individual = ind, cluster = g, n_leaves = k, clade_drawn = clade_drawn, split_unrooted = split_unrooted,
      n_pieces = length(pure_max), piece_sizes = paste(sort(sizes[pure_max], decreasing = TRUE), collapse = "+"),
      mrca_size = sizes[mrca], mrca_support = support[mrca], n_intruders = length(intr),
      intruders_by_cluster = intr_tab,
      intruder_leaves = paste(tips[intr][seq_len(min(12, length(intr)))], collapse = " "),
      stringsAsFactors = FALSE)
  }
}
R <- do.call(rbind, rows)
write.table(R, "out/ZCA_trees/zca_topology_clusters.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
S <- do.call(rbind, lapply(split(R, R$individual), function(d) {
  base <- d[!grepl("\\(", d$cluster), ]
  data.frame(individual = d$individual[1], clusters = nrow(base),
             clades_drawn = sum(base$clade_drawn), splits_unrooted = sum(base$split_unrooted),
             not_clades = paste(base$cluster[!base$split_unrooted], collapse = "; "),
             parent_groups = paste(sprintf("%s:%s", d$cluster[grepl("\\(", d$cluster)],
                                           ifelse(d$split_unrooted[grepl("\\(", d$cluster)], "clade", "not clade")), collapse = "; "))
}))
write.table(S, "out/ZCA_trees/zca_topology_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(S, row.names = FALSE)
cat("\n")
print(R[, c("individual", "cluster", "n_leaves", "clade_drawn", "split_unrooted", "n_pieces", "piece_sizes",
            "mrca_size", "mrca_support", "intruders_by_cluster")], row.names = FALSE)

# ---- convexity: can the partition be made by cutting (k - 1) edges? ----------
# Exact Sankoff parsimony with unit costs (score is root-independent, handles the
# basal trifurcation). Score == k - 1  <=>  every cluster is one connected piece
# of the unrooted tree and the clusters are cut apart by k - 1 edges. The
# backtracked optimal assignment gives those edges; each cut edge's support is
# the node label of its lower node.
sankoff <- function(tr, states) {
  n <- Ntip(tr); S <- sort(unique(states)); m <- length(S)
  post <- rev(postorder(tr)); kids <- split(tr$edge[, 2], tr$edge[, 1])
  cost <- matrix(Inf, n + tr$Nnode, m); cost[cbind(seq_len(n), match(states, S))] <- 0
  for (v in rev(unique(tr$edge[post, 1]))) {         # parents after children
    cost[v, ] <- 0
    for (c in kids[[as.character(v)]]) cost[v, ] <- cost[v, ] + sapply(seq_len(m), function(s) min(cost[c, ] + (seq_len(m) != s)))
  }
  root <- n + 1; st <- integer(n + tr$Nnode); st[root] <- which.min(cost[root, ])
  for (i in rev(post)) { p <- tr$edge[i, 1]; c <- tr$edge[i, 2]
    st[c] <- which.min(cost[c, ] + (seq_len(m) != st[p])) }   # ties keep parent state when optimal
  list(score = min(cost[root, ]), k = m, state = S[st])
}
conv <- list()
for (ind in names(TREES)) {
  tr <- read.tree(file.path(KLEIN, "trees2025_ordering", paste0(TREES[[ind]], "_bootstrapped_rerooted.newick")))
  z <- setNames(ml$ZCA, ml$new_chosen_sr)[tr$tip.label]
  sk <- sankoff(tr, z)
  sup <- c(rep(NA, Ntip(tr)), suppressWarnings(as.numeric(tr$node.label)))
  cut <- which(sk$state[tr$edge[, 1]] != sk$state[tr$edge[, 2]])
  desc_n <- sapply(tr$edge[cut, 2], function(v) if (v <= Ntip(tr)) 1 else Ntip(extract.clade(tr, v)))
  conv[[ind]] <- data.frame(
    individual = ind, clusters = sk$k, parsimony = sk$score, convex = sk$score == sk$k - 1,
    extra_changes = sk$score - (sk$k - 1),
    cut_edges = paste(sprintf("%s|%s (%d leaves, support %s)", sk$state[tr$edge[cut, 1]], sk$state[tr$edge[cut, 2]],
                              desc_n, ifelse(is.na(sup[tr$edge[cut, 2]]), "tip", sprintf("%.2f", sup[tr$edge[cut, 2]]))),
                      collapse = "; "),
    stringsAsFactors = FALSE)
}
C <- do.call(rbind, conv)
write.table(C, "out/ZCA_trees/zca_topology_convexity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n"); print(C[, c("individual", "clusters", "parsimony", "convex", "extra_changes")], row.names = FALSE)
