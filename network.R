library(igraph)
source("http://michael.hahsler.net/SMU/ScientificCompR/code/map.R")

g <- barabasi.game(1000, power=1)
layout <- layout.fruchterman.reingold(g)
layout <- layout.sphere(g)
layout <- layout.circle(g)
plot(g, layout=layout, vertex.size=2,
     vertex.label=NA, edge.arrow.size=.2)


degree(g)

B=betweenness(g)
EB=edge.betweenness(g)
plot(density(EB))

plot(g, layout=layout,
     vertex.size=map(betweenness(g),c(1,15)),
     edge.width=map(edge.betweenness(g), c(1,10)),vertex.label=NA,)

eb <- edge.betweenness.community(g)

#member <- community.to.membership(g, eb$merges, step=nrow(eb$merges)-10L+1L)

member=membership(eb)

plot(g,
     vertex.color= rainbow(10, .8, .8, alpha=.8)[member+2],
     vertex.size=5, layout=layout,  vertex.label=NA,
     edge.arrow.size=.2)


ec <- evcent(g)$vector
plot(g, layout=layout, vertex.size=map(ec, c(1,20)), vertex.label=NA, edge.arrow.size=.2)



