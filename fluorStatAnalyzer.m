%Color Map Analysis
%This script converts seedwater color map output into data

clear cells

%% Inputs
inputfolder = ['Data' filesep];
%pixelRatio = 2.919708029; % 2.919708029 => 40x
pixelRatio = 4.280821918; % 4.280821918 => 60x
%mapName = 'Embryo 1 cycE PH3 Stage 13 7Nov14 modified.png';
%fluorName1 = 'Embryo 1 PH3 en cycE.tif';
%fluorName2 = 'Embryo 1 GFP en cycE.tif';
%compositeName = 'Composite (RGB).tif';
%segmentationName = 'b.tiff';

mapName = 'Embryo 1 Colormap.tif';
% fluorName1 = 'dpERK.tif';
% fluorName2 = 'engrailed.tif';



compositeName = 'Embryo 1.bmp';
% segmentationName = 'cell boundaries.tif';

% segmentationFluor = imresize(imread([inputfolder segmentationName]),size(label));
dateString = datestr(now);
dateString (dateString == ':') = '.';
mapScale = 3;

% segmentListMaster = {'Head';'Amnioserosa';'Ventral Segments T1 through A1'; ...
%     'Ventral Segments A2 through a5';'Ventral Segments A6 through A9';'C5'; ...
%     'dorsal T1';'dorsal T2';'dorsal T3';'dorsal A1';'dorsal A2';'dorsal A3'; ...
%     'dorsal A4';'dorsal A5'};
% 
% %fbMAT = {''};
% fbMAT = {' front'; ' back'};
% 
% APMAT = {'Anterior'; 'Posterior'};
% segmentListMaster = {'Head', 'Amnioserosa'};

APMAT = {'Anterior'; 'Posterior'};
segmentListMaster = {};
for segA = 1:12
    for AP = 1:2
        %for fb = 1:2
        %tempSTR = ['A' int2str(segA) ' - ' APMAT{AP} fbMAT{fb}]
        tempSTR = ['Seg' int2str(segA) ' - ' APMAT{AP}];
        segmentListMaster = [segmentListMaster {tempSTR}];
        %end
    end
end

% for segA = 1:10
%     for AP = 1:2
%         for fb = 1:2
%         %tempSTR = ['A' int2str(segA) ' - ' APMAT{AP} fbMAT{fb}]
%         tempSTR = ['Seg' int2str(segA) ' - ' APMAT{AP} fbMAT{fb}];
%         segmentListMaster = [segmentListMaster {tempSTR}];
%         end
%     end
% end


%segmentListMaster = {'Head';'Body'};

% Initialize output folder
outputfolder = [inputfolder 'Segmentation Analysis Output ' dateString '\'];
mkdir(outputfolder);

% Import matricies
Map = imread([inputfolder mapName]);
%Map = label;
colorByNumbers = labelMaker(Map);
colorByNumbers(colorByNumbers == 0) = max(colorByNumbers(:)) + 1;
cellsNumber = max(colorByNumbers(:)) - 1;
% compositeMap(:,:,1) = fluorImage1;
% compositeMap(:,:,2) = fluorImage2;
% compositeMap(:,:,3) = segmentationFluor;
compositeMap = imread([inputfolder compositeName]);
fluorImage1 = squeeze(compositeMap(:,:,1));
fluorImage2 = squeeze(compositeMap(:,:,2));

compositeMap = uint8(compositeMap);
fluorImage1 = imresize(fluorImage1, size(colorByNumbers));
fluorImage2 = imresize(fluorImage2, size(colorByNumbers));

%% Obtain Region Properties

STATS = regionprops(colorByNumbers, fluorImage1, 'all');
STATS2 = regionprops(colorByNumbers, fluorImage2, 'MeanIntensity');

Centroid = zeros(cellsNumber, 2);

[neigborList, neigborCount, edgecellss] = findNeigbors(STATS, colorByNumbers);

%% Create cells structure and store statistics


for i = cellsNumber:-1:1
    cells(i).area = STATS(i).Area / pixelRatio^2;
    cells(i).dpERKIntensity = STATS(i).MeanIntensity;
    cells(i).aspectRatio = STATS(i).MajorAxisLength / STATS(i).MinorAxisLength;
    cells(i).perimeter = STATS(i).Perimeter / pixelRatio;
    cells(i).convexArea = STATS(i).ConvexArea / pixelRatio^2;
    cells(i).solidity = cells(i).convexArea / cells(i).area;
    cells(i).majorAxis = STATS(i).MajorAxisLength / pixelRatio;
    cells(i).minorAxis = STATS(i).MinorAxisLength / pixelRatio;
    cells(i).engrailedLevel = STATS2(i).MeanIntensity;
    cells(i).polygonClass = neigborCount(i);
    cells(i).edgecells = edgecellss(i);
    cells(i).noAverage.number = i;
    cells(i).noAverage.centroid = STATS(i).Centroid;
    Centroid(i,:) = STATS(i).Centroid;
end

%% Calculate Regional Averages and Neighbor Statistics

cellFields = fieldnames(cells);

for m = 1:length(cells)
    for i = 2:length(cellFields)
        for j = 1:(i-1)
            if ~strcmp(cellFields{i},'noAverage') && ~strcmp(cellFields{j},'noAverage')
                cotemp = cells(m).(cellFields{i}) ./ cells(m).(cellFields{j});
                cells(m).([cellFields{i} 'over' cellFields{j}]) = cotemp;
            end
        end
    end
end


% for k = 1:numel(cellFields)
%     if ~strcmp(cellFields{k},'noAverage')
%         for i = cellsNumber:-1:1
%             cells(i).noAverage.neighborAve.(cell2mat(cellFields(k))) = cells(neigborList(i,1)).(cell2mat(cellFields(k)));
%             cells(i).noAverage.regionalAve.(cell2mat(cellFields(k))) = cells(neigborList(i,1)).(cell2mat(cellFields(k))) + cells(i).(cell2mat(cellFields(k)));
%             for j = neigborList(i,2:neigborCount(i))
%                 cells(i).noAverage.neighborAve.(cell2mat(cellFields(k))) = cells(i).noAverage.neighborAve.(cell2mat(cellFields(k))) + cells(j).(cell2mat(cellFields(k)));
%                 cells(i).noAverage.regionalAve.(cell2mat(cellFields(k))) = cells(i).noAverage.regionalAve.(cell2mat(cellFields(k))) + cells(j).(cell2mat(cellFields(k)));
%             end
%             cells(i).noAverage.neighborAve.(cell2mat(cellFields(k))) = cells(i).noAverage.neighborAve.(cell2mat(cellFields(k))) ./ neigborCount(i);
%             cells(i).noAverage.regionalAve.(cell2mat(cellFields(k))) = cells(i).noAverage.regionalAve.(cell2mat(cellFields(k))) ./ (neigborCount(i)+1);
%             cells(i).noAverage.relativeValue.(cell2mat(cellFields(k))) = cells(i).(cell2mat(cellFields(k))) ./ cells(i).noAverage.neighborAve.(cell2mat(cellFields(k)));
%         end
%     end
% end



%% Create connection matricies for network statistics

D = pdist(cat(1,Centroid(:,1),Centroid(:,2)));
distance = squareform(D) * pixelRatio;

cMatrixUnweighted = zeros(cellsNumber);
cMatrixWeighted = zeros(cellsNumber);

for i = 1:cellsNumber
    neigbors = neigborList(i,1:neigborCount(i));
    for j = neigbors(neigbors>0)
        cMatrixUnweighted(i,j) = 1;
        cMatrixWeighted(i,j) = distance(i,j);
    end
end

save([outputfolder 'cMatrixUnweighted.mat'], 'cMatrixUnweighted');
save([outputfolder 'cMatrixWeighted.mat'], 'cMatrixWeighted');


%% Perform Segmentation




segmentBelongedTo = cell(cellsNumber, 1);
figureHandle = imshow(compositeMap);

x = 1; y = 1; notDone = true; segmentList = segmentListMaster;
identificationArray = cell(1, cellsNumber);


while notDone && ~isempty(segmentList)
    
    [selection, notDone] = listdlg('ListString',segmentList,'SelectionMode','single','PromptString','Name of Current Segment');
    currentSegment = segmentList{selection};
    segmentList(selection) = [];
    
    if notDone
        [x,y] = ginput();
        [ axyp, ~, ~ ] = polygeom( x, y );
        compositeMap = insertText(compositeMap, axyp(2:3), currentSegment, 'FontSize', 18, 'TextColor', ...
            [255 255 255], 'AnchorPoint', 'Center', 'BoxOpacity', 0);
        figureHandle = imshow(compositeMap);
        in = inpoly(Centroid, [x y]);
        identificationArray(in) = {currentSegment};
    end
    
end

unidentified = double(cellfun('isempty',identificationArray));
unidentified = unidentified .* double(1:cellsNumber);
unidentified(unidentified == 0) = [];

while ~isempty(unidentified)
    for i = unidentified
        identities = cell(1, neigborCount(i));
        for count = 1:neigborCount(i)
            identities(count) = identificationArray(neigborList(i,count));
        end
        identities(cellfun('isempty',identities)) = [];
        [~, numNeighborIdents] = size(identities);
        if numNeighborIdents > 0
            uniIdent = unique(identities);
            count = zeros(length(uniIdent), 1);
            for iy = 1:length(uniIdent)
                count(iy) = length(find(strcmp(uniIdent{iy}, identities)));
            end
            [~, itemp] = max(count);
            identificationArray(i) = uniIdent(itemp);
            unidentified(unidentified == i) = [];
        end
    end
end

for i = cellsNumber:-1:1
    cells(i).segmentBelongedTo = identificationArray{i};
    cells(i).subSegmentBelongedTo = cells(i).segmentBelongedTo;
end

%% Create cell-level macrosegmentation map

for i = length(segmentListMaster):-1:1
    macroSeg(ismember(identificationArray, segmentListMaster{i})) = i;
end

createScatterPlots(Centroid(:,1), Centroid(:,2), macroSeg, 'Segmentation Map', '')


cellFields = fieldnames(cells);

identities = unique(identificationArray);

 for i = 2:length(cellFields)
        for j = 1:(i-1)
            if ~strcmp(cellFields{i},'noAverage') && ~strcmp(cellFields{j},'noAverage')
            cotemp = corrcoef(double([cells.(cellFields{i})]),double([cells.(cellFields{j})]));
            population.(['corr' int2str(i) 'vs' int2str(j)]) = cotemp(1,2);
            end
            end
    end

for currentPop = 1:max(macroSeg);
    for i = 1:length(cellFields)
        if ~strcmp(cellFields{i},'noAverage')
            population(currentPop).(['average' cellFields{i}]) = mean([cells(macroSeg==currentPop).(cellFields{i})]);
            population(currentPop).(['stdev' cellFields{i}]) = std([cells(macroSeg==currentPop).(cellFields{i})]);
        end
    end
    
   
end


%% Run PCA Analysis

cellPCA(cells);

%% Identify sub-segments computationally

% segmentsToSplitEngrailed = {'Body','Head'};
%  segmentsToSplitEngrailed = {'C5'; ...
%      'dorsal T1';'dorsal T2';'dorsal T3';'dorsal A1';'dorsal A2';'dorsal A3'; ...
%      'dorsal A4';'dorsal A5'};
%
% identificationSubArray = cell(1, cellsNumber);
%
% for currentSegmentNumber = 1:length(segmentsToSplitEngrailed)
%     currentSegment = segmentsToSplitEngrailed(currentSegmentNumber);
%     currentCells = find(ismember(identificationArray, currentSegment));
%     activeCells = cells(currentCells);
%
%     X = [];
%
%     for i = length(activeCells):-1:1
%         X(i,1:2) = activeCells(i).noAverage.centroid;
%         X(i,3) = activeCells(i).engrailedLevel;
%     end
%
%     cluster = kmeans(X,2);
%     pCluster = find(ismember(cluster, 1));
%     aCluster = find(ismember(cluster, 2));
%
%     PNumber = [];
%     ANumber = [];
%
%     for i = pCluster'
%         PNumber = cat(2,PNumber,activeCells(i).noAverage.number);
%     end
%
%     for i = aCluster'
%         ANumber = cat(2,ANumber,activeCells(i).noAverage.number);
%     end
%
%     for i = PNumber
%       cells(i).subSegmentBelongedTo = {[cells(i).subSegmentBelongedTo ' P Compartment']};
%     end
%
%     for i = ANumber
%         cells(i).subSegmentBelongedTo = {[cells(i).subSegmentBelongedTo ' A Compartment']};
%     end
%
% end

%% Run all possible p tests

%p tests between all properties for every segment
p_dpERK = cellTTest([cells.dpERKIntensity],[{cells.segmentBelongedTo}]);
p_engrailedLevel = cellTTest([cells.engrailedLevel],[{cells.segmentBelongedTo}]);
p_area = cellTTest([cells.area],[{cells.segmentBelongedTo}]);
p_polygonClass = cellTTest([cells.polygonClass],[{cells.segmentBelongedTo}]);

p_sub_dpERK = cellTTest([cells.dpERKIntensity],[cells.subSegmentBelongedTo]);
p_sub_engrailedLevel = cellTTest([cells.engrailedLevel],[cells.subSegmentBelongedTo]);
p_sub_area = cellTTest([cells.area],[cells.subSegmentBelongedTo]);
p_sub_polygonClass = cellTTest([cells.polygonClass],[cells.subSegmentBelongedTo]);



save([outputfolder 'pTests.mat'], 'p_*');

%% Output Tables

OutputTable = cell(cellsNumber+1, numel(cellFields)+1);
OutputTable(2:(cellsNumber+1),1) = num2cell(1:cellsNumber);

OutputTableNeighbor = OutputTable;
OutputTableRegional = OutputTable;
OutputTableRelative = OutputTable;

cellFields = fieldnames(cells);

for k = numel(cellFields):-1:1
    if ~strcmp(cellFields{k},'noAverage')
        for i = cellsNumber:-1:1
            
            OutputTable(i+1,k+1) = num2cell(cells(i).(cell2mat(cellFields(k))));
            OutputTableNeighbor(i+1,k+1) = num2cell(cells.noAverage.neighborAve(i).(cell2mat(cellFields(k))));
            OutputTableRegional(i+1,k+1) = num2cell(cells.noAverage.regionalAve(i).(cell2mat(cellFields(k))));
            OutputTableRelative(i+1,k+1) = num2cell(cells.noAverage.relativeValue(i).(cell2mat(cellFields(k))));
        end
        OutputTable(1,k+1) = cellFields(k);
        OutputTableNeighbor(1,k+1) = cellFields(k);
        OutputTableRegional(1,k+1) = cellFields(k);
        OutputTableRelative(1,k+1) = cellFields(k);
    end
end

xlswrite([outputfolder 'Cell Level Statistics.xls'], OutputTable, 'Cell Statistics');
xlswrite([outputfolder 'Cell Level Statistics.xls'], OutputTableNeighbor, 'Cell Neighbor Average');
xlswrite([outputfolder 'Cell Level Statistics.xls'], OutputTableRegional, 'Cell Regional Average');
xlswrite([outputfolder 'Cell Level Statistics.xls'], OutputTableRelative, 'Cell Relative Statistics');

%% Create annotated segmentation map
Map = imresize(Map, mapScale);
[a,b,c] = size(Map);

mapPosition = round(Centroid * mapScale);

for i = 1:(cellsNumber - 1)
    Map = insertText(Map, mapPosition(i,:), i, 'FontSize', 10, 'TextColor', ...
        [0 0 0], 'AnchorPoint', 'Center', 'BoxOpacity', 0);
end

imwrite(Map, [outputfolder 'Segmentation Map.png']);
imwrite(compositeMap, [outputfolder 'Macro Segmentation Map.png']);