%function cellPCA( cellMatrix )

clear dataMatrix;

cellMatrix = cells;
frac1 = struct2cell(cellMatrix);
noave = frac1(11,:,:);
noave = noave(:);
noave = cell2mat(noave);

frac2 = struct2cell([noave.neighborAve]);
frac3 = struct2cell([noave.regionalAve]);
frac4 = struct2cell([noave.relativeValue]);
frac1((end-2):end,:,:) = [];

frac1 = permute(frac1,[1 3 2]);
frac2 = permute(frac2,[1 3 2]);
frac3 = permute(frac3,[1 3 2]);
frac4 = permute(frac4,[1 3 2]);
  
rawDataMatrix = [frac1; frac2; frac3; frac4];
[x, y] = size(rawDataMatrix);

for i = x:-1:1
    for j = y:-1:1
        dataMatrix(i,j) = double(rawDataMatrix{i,j});
    end
end

dataMatrix(:,dataMatrix(10,:) == 1) = [];

dataMatrix(end,:) = [];
dataMatrix = dataMatrix';
dataMatrix(:,10) = [];

w = 1./var(dataMatrix);
[wcoeff,score,latent,tsquared,explained] = pca(dataMatrix,'VariableWeights',w);

c3 = wcoeff(:,1:3);
coefforth = inv(diag(std(dataMatrix)))*wcoeff;

colormap(cool(255));

clusterScores = zeros(0);

%clusterScores = dataMatrix(:,2);
for k = 1:length(cells)
    if ~cells(k).edgecells
        clusterScores = [clusterScores cells(k).noAverage.relativeValue.dpERKIntensity];
        %clusterScores = [clusterScores cells(k).dpERKIntensity];
    end
end

percentilesScore = prctile(clusterScores,[33,67]);

binScores(clusterScores < percentilesScore(1)) = 1;
binScores(clusterScores > percentilesScore(2)) = 3;
binScores(binScores == 0) = 2;

for i = 1:5
    for j = (i+1):5
        figure()
        scatter(score(:,i),score(:,j),8,binScores,'fill')
        xlabel(['Principal Component ' int2str(i)])
        ylabel(['Principal Component ' int2str(j)])
        %zlabel('3rd Principal Component')
        title('Principle Component Analysis - Relative EGFR Signaling')
        c = colorbar('SouthOutside');
        colormap(redgreencmap(255));

        print (['PCA relative EGFR ' int2str(i) ' ' int2str(j)],'-r600', '-dtiff');
    end
end
%end

