addpath('matlab_bgl');

sparseUW = sparse(cMatrixUnweighted);
sparseW = sparse(cMatrixWeighted);

betweenessCentralityUW = betweenness_centrality(sparseUW);
betweenessCentralityW = betweenness_centrality(sparseW);

for i = 1:cellsNumber
    
    neighborstemp = neigborList(i,1:neigborCount(i));
    
    cells(i).strength = sum(sparseW(i,neighborstemp));
    tempNeighbor = 0;
    tempClusterNeighbor = 0;
    
    for j = neighborstemp
        neighbOfNeigb = neigborList(j,1:neigborCount(j));
        neighbOfNeigb(neighbOfNeigb==i)=[];
        clusterTemp = ismember(neighbOfNeigb,neighborstemp);
        
        tempNeighbor = tempNeighbor + (sum(clusterTemp)>0);
        tempClusterNeighbor = tempClusterNeighbor + sum(clusterTemp)/length(clusterTemp);
    end
    
    cells(i).clusteringCoefficient = tempNeighbor / neigborCount(i);
    cells(i).betweenessCentralityUW = betweenessCentralityUW(i);
    cells(i).betweenessCentralityW = betweenessCentralityW(i);
    cells(i).clusteringCoefficient = tempNeighbor / length(neighborstemp);
    cells(i).clusteringFraction = tempClusterNeighbor / length(neighborstemp);
    cells(i).rosetteParameter = tempClusterNeighbor;

end