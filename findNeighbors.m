function [ neigborList, neigborCount, edgeCells ] = findNeigbors( STATS, labelMatrix )

cellNumber = max(labelMatrix(:));
neigborCount = zeros(cellNumber,1);
neigborList = zeros(cellNumber,1);
edgeCells = false(cellNumber, 1);
[u, n] = size(labelMatrix);
labelMatrix(2:(u+1),2:(n+1))=labelMatrix;
labelMatrix(u+2,n+2)=0;
blankLogicMap = false(size(labelMatrix));
logicMap = blankLogicMap;

for i = 1:(cellNumber-1)
    x = STATS(i).BoundingBox;
    y = round(x(2));
    x = round(x(1));
    if x < 2
        x = 2;
    end
    if y < 2
        y = 2;
    end
    temp = STATS(i).Image;
    bwImage = zeros(size(temp)+2);
    [a, b] = size(bwImage);
    bwImage (2:(a-1),2:(b-1)) = temp;
    bwImage = logical(imdilate(bwImage,strel('diamond',1)));
    logicMap((y):(y+a-1),(x):(x+b-1)) = bwImage;
    [x,y] = size(labelMatrix);
    cells = labelMatrix(logicMap(1:x,1:y));
    if max(cells(:)) == cellNumber
        edgeCells(i) = true;
    end
    if ~isempty(cells)
    cells (cells == i) = cellNumber;
    cells = unique(cells);
    cells(cells==cellNumber) = [];
    
    [neigborCount(i), ~] = size(cells);
    neigborList(i, 1:neigborCount(i)) = cells;
    logicMap = blankLogicMap;
    end
end

end

