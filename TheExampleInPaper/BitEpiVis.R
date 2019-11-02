library(dplyr)
library(RCy3)
library(igraph)
setwd('~/temp/cc/BitEpi')

Color=list(SNP='red',PAIR='blue',TRIPLET='orange',QUADLET='green', OTHER='gray')

# Nodes of the graph are SNPs and Interactions
# Each SNP node could be connected to multiple Interaction Node
# Each Interaction Node is conneced to the SNPs that are involved in that interaction. 
# This function name the interaction nodes by concatinating SNPS with # seprator.
# the 2-SNP, 3-SNP, and 4-SNP names are added as 3 new column to the best dataframe
# For Example if rs123, rs456 and rs789 Interact with each other then
# the Interaction node is called rs123#rs456#rs789
AddInteractionNode = function(data)
{
  x = as.data.frame(t(apply(select(data, SNP, PAIR), 1, sort)))
  data$nP = paste(x$V1, x$V2, sep = "#")
  x = as.data.frame(t(apply(select(data, SNP, TRIPLET_1, TRIPLET_2), 1, sort)))
  data$nT = paste(x$V1, x$V2, x$V3, sep = "#")
  x = as.data.frame(t(apply(select(data, SNP, QUADLET_1, QUADLET_2, QUADLET_3), 1, sort)))
  data$nQ = paste(x$V1, x$V2, x$V3, x$V4, sep = "#")
  return(data)
}

# list all the nodes (1-SNP, 2-SNP, 3-SNP, 4SNP) assing beta to the size and rank them by order
NodeGen = function(dataX)
{
  #1-SNP
  data       = dataX
  data$Node  = data$SNP
  data$order = 1
  data$beta = data$SNP_B
  data$color = Color$SNP
  data       = data[order(-data$SNP_A),]
  data$rank  = seq.int(nrow(data))
  nodes      = select(data, Node, rank, beta, color, order)

  #2-SNP
  data       = dataX
  data$Node  = data$nP
  data$order = 2
  data$beta  = data$PAIR_B
  data$color = Color$PAIR
  data       = data[order(data[,'Node'],-data[,'beta']),]
  data       = data[!duplicated(data$Node),]
  data       = data[order(-data$PAIR_A),]
  data$rank  = seq.int(nrow(data))
  nodes      = rbind(nodes, select(data, Node, rank, beta, color, order))
  
  #3-SNP
  data       = dataX
  data$Node  = data$nT
  data$order = 3
  data$beta  = data$TRIPLET_B
  data$color = Color$TRIPLET
  data       = data[order(data[,'Node'],-data[,'beta']),]
  data       = data[!duplicated(data$Node),]
  data       = data[order(-data$TRIPLET_A),]
  data$rank  = seq.int(nrow(data))
  nodes      = rbind(nodes, select(data, Node, rank, beta, color, order))
  
  #4-SNP
  data       = dataX
  data$Node  = data$nQ
  data$order = 4
  data$beta  = data$QUADLET_B
  data$color = Color$QUADLET
  data       = data[order(data[,'Node'],-data[,'beta']),]
  data       = data[!duplicated(data$Node),]
  data       = data[order(-data$QUADLET_A),]
  data$rank  = seq.int(nrow(data))
  nodes      = rbind(nodes, select(data, Node, rank, beta, color, order))

  return(nodes)
}

# list all edges between interactive nodes (2-SNP, 3-SNP and 4-SNP) and SNP nodes (1-SNP)
EdgeGen = function(data)
{
  edf = data.frame(source=character(), target=character())
  
  for(i in 1:nrow(data)) {
    edf = rbind(edf, data.frame(source=data[i,"nP"], target=data[i,"SNP"]))
    edf = rbind(edf, data.frame(source=data[i,"nP"], target=data[i,"PAIR"]))
    edf = rbind(edf, data.frame(source=data[i,"nT"], target=data[i,"SNP"]))
    edf = rbind(edf, data.frame(source=data[i,"nT"], target=data[i,"TRIPLET_1"]))
    edf = rbind(edf, data.frame(source=data[i,"nT"], target=data[i,"TRIPLET_2"]))
    edf = rbind(edf, data.frame(source=data[i,"nQ"], target=data[i,"SNP"]))
    edf = rbind(edf, data.frame(source=data[i,"nQ"], target=data[i,"QUADLET_1"]))
    edf = rbind(edf, data.frame(source=data[i,"nQ"], target=data[i,"QUADLET_2"]))
    edf = rbind(edf, data.frame(source=data[i,"nQ"], target=data[i,"QUADLET_3"]))
  }
  return(edf)
}

# convert BitEpi Best file to nodes and edges
BestToNodesAndEdges = function(bestFn)
{
  # Read BitEpi "best" file into a data frame
  bestDf = read.csv(bestFn)
  
  bestDf = AddInteractionNode(bestDf)
  
  Nodes = NodeGen(bestDf)
  
  Edges = EdgeGen(bestDf)
  
  return(list(Nodes=Nodes, Edges=Edges))
}

# query nodes and related edges 
QueryGraph = function(Graph, thr, minNodeSize, maxNodeSize)
{
  if(minNodeSize >= maxNodeSize)
  {
    print("minNodeSize is greater or equal maxNodeSize")
    return(NULL,NULL)
  }
  
  allNodes = Graph$Nodes
  allEdges = Graph$Edges
  
  # select nodes to be in the graph
  s1 = allNodes %>% filter(order==1 & allNodes$rank<=thr$SNP)
  s2 = allNodes %>% filter(order==2 & allNodes$rank<=thr$PAIR)
  s3 = allNodes %>% filter(order==3 & allNodes$rank<=thr$TRIPLET)
  s4 = allNodes %>% filter(order==4 & allNodes$rank<=thr$QUADLET)
  selNodes = unique(rbind(s1,s2,s3,s4))
  
  # select interaction nodes
  intNodes = selNodes %>% filter(order>1)
  
  # select all edges for intraction nodes
  intEdges = select(merge(x=allEdges, y=intNodes, by.x='source', by.y='Node'), source, target)
  intEdges = unique(intEdges)
  
  # select all target names for interaction nodes
  tarNames = unique(select(intEdges, target))
  names(tarNames) = 'Node'
  #grab target nodes from all nodes
  tarNodes = merge(x=allNodes, y=tarNames, by='Node')
  
  
  #and merge them to dataset
  selNodes = unique(rbind(selNodes, tarNodes))

  minBeta = min(selNodes$beta)
  maxBeta = max(selNodes$beta)
  rangeBeta = maxBeta - minBeta
  rangeSize = maxNodeSize - minNodeSize
  ratio = rangeSize/rangeBeta
  selNodes$size = ((selNodes$beta - minBeta) * ratio) + minNodeSize;
  
  selNodes[(selNodes$order==1) & (selNodes$rank>thr$SNP),]$color = Color$OTHER
  selNodes[(selNodes$order==1) & (selNodes$rank>thr$SNP),]$size = minNodeSize
  
  return(list(Nodes=selNodes, Edges=intEdges))
}

DoItAll = function(bestFn, thr, minNodeSize, maxNodeSize)
{
  # read best file into a graph
  GraphAll = BestToNodesAndEdges(bestFn)
  
  # query graph
  GraphSelected = QueryGraph(GraphAll, thr, minNodeSize, maxNodeSize)
  
  Edges = GraphSelected$Edges
  Nodes = GraphSelected$Nodes
  
  #plot graph 
  Nodes$label = " "
  network = graph_from_data_frame(d=Edges, directed=FALSE, vertices = Nodes)
  plot(network, vertex.size=V(network)$size, vertex.label=V(network)$Node, vertex.color=V(network)$color, vertex.label=V(network)$label)
  
  cytoscapePing()
  createNetworkFromIgraph(network,"BitEpi Network", title = "BitEpi Graph")
}

thr=list(SNP=3,PAIR=3,TRIPLET=3,QUADLET=3)
minNodeSize = 10
maxNodeSize = 35

#Sort by Alpha and but represent beta as node size in the plot
DoItAll('TheExampleInPaper/data.best.csv', thr, minNodeSize, maxNodeSize)
