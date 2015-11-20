function cellPCA( cellMatrix )

identifier = {cellMatrix.subSegmentBelongedTo};
cellFields = fieldnames(cells);
cellFields(ismember(cellFields, 'subSegmentBelongedTo')) = [];
cellFields(ismember(cellFields, 'segmentBelongedTo')) = [];
cellFields(ismember(cellFields, 'noAverage')) = [];

cellAverageFields = fieldnames(cells.noAverage.neighborAve);

uniqueIdentifier = unique(identifier);
uniqueIdentifier(ismember(uniqueIdentifier, 'Head')) = [];
uniqueIdentifier(ismember(uniqueIdentifier, 'Other')) = [];

rawDataMatrix = zeros(0,10);
rowIdentifier = {};
rowValue = zeros(0);

for k = numel(cellFields):-1:1
    firstNumber = 1;
    for i = length(uniqueIdentifier):-1:1
        temp = ismember({cellMatrix.subSegmentBelongedTo}, uniqueIdentifier{i});
        number = sum(temp);
        rawDataMatrix(firstNumber:(firstNumber-1+number),k) = [cellMatrix(temp).(cellFields{k})];
        for m = 1:number
            rowIdentifier(firstNumber+m-1) = {[uniqueIdentifier{i} ' ' num2str(m,['%0' int2str(numel(int2str(number))) 'd'])]};
        rowValue(firstNumber+m-1) = i;
       
        end
        firstNumber = firstNumber + number;
        
    end
end

for k = (numel(cellFields)+1):3:((numel(cellFields)+1)+numel(cellAverageFields)*3)
    
end

w = 1./var(rawDataMatrix);
[wcoeff,score,latent,tsquared,explained] = pca(rawDataMatrix,'VariableWeights',w);

c3 = wcoeff(:,1:3);
coefforth = inv(diag(std(rawDataMatrix)))*wcoeff;

colormap(jet(length(rowIdentifier)));

figure()
scatter(score(:,1),score(:,2),32,rowValue,'fill')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')

title('Principle Component Analysis')
c = colorbar('SouthOutside');
%c.TickLabels = uniqueIdentifier;

print ('PCA Test 2','-r600', '-dtiff');

end

