function [ pTable ] = cellTTest(values, classification)

classList = unique(classification);
[~, classNumber] = size(classList);

pTable = zeros(classNumber);

for i = 1:classNumber
    for j = 1:classNumber
        x = values(ismember(classification, classList(i)));
        y = values(ismember(classification, classList(j)));
        [~, pTemp] = ttest2(x, y);
        pTable(i, j) = pTemp;
    end
end

end

