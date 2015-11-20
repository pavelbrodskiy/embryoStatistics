function [ Assortativity ] = subdivideNetwork( connectionMatrix, regionLabels )

uniqueRegionLabels = unique(regionLabels);

for i = 1:length(uniqueRegionLabels)

    cellsToKeep = find(strcmp(regionLabels, uniqueRegionLabels{i}));
    subConnectionMatricies = connectionMatrix(cellsToKeep, cellsToKeep);  
    Assortativity(i) = assortativity(logical(subConnectionMatricies), 0);

end

end