transitive_closure <- function(edges, nodes = NULL) {
stopifnot(ncol(edges) == 2)
if (is.data.frame(edges))
	edges <- as.matrix(edges)
if (is.null(nodes)) {
	nodes <- unique(as.vector(edges))	
} else {
	stopifnot(all(as.vector(edges) %in% nodes))
}
nnodes <- length(nodes)
if (is.numeric(nodes) && !setequal(nodes, 1:nnodes)) {
	edges <- matrix(as.character(edges), nrow(edges), 2)
	nodes <- as.character(nodes)
}
comp_vec <- integer(nnodes)
if (is.character(nodes)) {
	names(comp_vec) <- nodes	
} 
nedges <- nrow(edges)
ncomps <- 0
for (i in 1:nedges) {
	node1 <- edges[i,1]
	node2 <- edges[i,2]
	val1 <- comp_vec[node1]
	val2 <- comp_vec[node2]	
	if (val1 == 0 && val2 == 0) {
		ncomps <- ncomps + 1
		comp_vec[c(node1, node2)] <- ncomps
	} else if (val1 == 0) {
		comp_vec[node1] <- val2		
	} else if (val2 == 0) {
		comp_vec[node2] <- val1
	}
}
comp_vec[comp_vec == 0] <- seq(ncomps + 1, ncomps + sum(comp_vec == 0))
comp_list <- split(nodes, comp_vec)
comp_size <- sapply(comp_list, length)
comp_list <- comp_list[comp_size > 1]
closure <- lapply(comp_list, combn, m = 2)
closure <- matrix(unlist(closure), ncol = 2, byrow = TRUE)
list(edges = closure, comp_vec = comp_vec, ncomps = max(comp_vec))	
}
