%function cellDiscriminant( cellMatrix )
cellMatrix = cells;
identifier = {cellMatrix.subSegmentBelongedTo};
cellFields = fieldnames(cellMatrix);
cellFields(ismember(cellFields, 'subSegmentBelongedTo')) = [];
cellFields(ismember(cellFields, 'segmentBelongedTo')) = [];
cellFields(ismember(cellFields, 'noAverage')) = [];

%cellAverageFields = fieldnames(cellMatrix.noAverage.neighborAve);

uniqueIdentifier = unique(identifier);
%uniqueIdentifier(ismember(uniqueIdentifier, 'Head')) = [];
uniqueIdentifier(ismember(uniqueIdentifier, 'Other')) = [];

rawDataMatrix = zeros(0,10);
rowIdentifier = {};
rowValue = zeros(0);
rowIdentifierNumber = {};

for k = numel(cellFields):-1:1
    firstNumber = 1;
    for i = length(uniqueIdentifier):-1:1
        temp = ismember({cellMatrix.subSegmentBelongedTo}, uniqueIdentifier{i});
        number = sum(temp);
        rawDataMatrix(firstNumber:(firstNumber-1+number),k) = [cellMatrix(temp).(cellFields{k})];
        for m = 1:number
            rowIdentifierNumber(firstNumber+m-1) = {[uniqueIdentifier{i} ' ' num2str(m,['%0' int2str(numel(int2str(number))) 'd'])]};
        rowIdentifier(firstNumber+m-1,1) = uniqueIdentifier(i);
        class(firstNumber+m-1,1) = i;
        rowValue(firstNumber+m-1) = i;
       
        end
        firstNumber = firstNumber + number;
        
    end
end

colormap(jet(length(rowIdentifier)));

print ('LDA Test 2','-r600', '-dtiff');

%end

