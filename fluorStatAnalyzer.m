%% fluorStatAnalyzer.m
% Color Map Analysis
% Updated 9/16/15
%
% This script converts seedwater color map output into data
% open fluorescent data which has already been stitched, so we have our
% images and stains, we have cell boundaries and things are already
% segmented so we have label matrix that says where cells are
% Now we do analysis - find which cells are adjacent to which cells,
% expression levels of proteins, shapes/sizes of cells

% want to ask the question: over the course of a compartment, does the
% front have bigger/same/smaller cell areas than back?

%% Inputs - open files
inputfolder = '~/Desktop/Zartman Research/Junior/Fluorescence Analysis/09.1.15 analysis code';
pixelRatio = 2.919708029; % 2.919708029 => 40x
mapName = 'Segmentation Map Modified.png';
fluorName1 = 'dpERK.tif';
fluorName2 = 'engrailed.tif';
compositeName = 'composite.tif';
segmentationName = 'cell boundaries.tif';
segmentationFluor = imread([inputfolder, '/', segmentationName]);
dateString = datestr(now);
dateString (dateString == ':') = '.';
mapScale = 3;

% segmentListMaster = {'Head';'Amnioserosa';'Ventral Segments T1 through A1'; ...
%      'Ventral Segments A2 through a5';'Ventral Segments A6 through A9';'C5'; ...
%      'dorsal T1';'dorsal T2';'dorsal T3';'dorsal A1';'dorsal A2';'dorsal A3'; ...
%      'dorsal A4';'dorsal A5'};

%fbMAT = {''};
%fbMAT = {' front'; ' back'};

%APMAT = {'Anterior'; 'Posterior'};

% for segA = 1:10
%     for AP = 1:2
%         %for fb = 1:2
%             %tempSTR = ['A' int2str(segA) ' - ' APMAT{AP} fbMAT{fb}]
%             tempSTR = ['Seg' int2str(segA) ' - ' APMAT{AP}];
%             segmentListMaster = [segmentListMaster {tempSTR}];
%         %end
%     end
% end

APMAT = {'Anterior'; 'Posterior'};
segmentListMaster = {'Head', 'Amnioserosa'};
for segA = 1:12
    for AP = 1:2
        %for fb = 1:2
        %tempSTR = ['A' int2str(segA) ' - ' APMAT{AP} fbMAT{fb}]
        tempSTR = ['Seg' int2str(segA) ' - ' APMAT{AP}];
        segmentListMaster = [segmentListMaster {tempSTR}];
        %end
    end
end
            
% segmentListMaster = {'Head';'Body'};

% Initialize output folder
outputfolder = [inputfolder, '/' 'Segmentation Analysis Output ' dateString '\'];
mkdir(outputfolder);

%% Import matricies (read in images)
Map = imread([inputfolder, '/', mapName]); % segmentation map
fluorImage1 = imread([inputfolder, '/', fluorName1]);
fluorImage2 = imread([inputfolder, '/', fluorName2]);
colorByNumbers = labelMaker(Map);
cellsNumber = max(colorByNumbers(:)) - 1;
compositeMap = imread([inputfolder, '/', compositeName]);

%% Obtain Region Properties (morphometry, shapes/sizes, centroid)
STATS = regionprops(colorByNumbers, fluorImage1, 'all');
STATS2 = regionprops(colorByNumbers, fluorImage2, 'MeanIntensity');

Centroid = zeros(cellsNumber, 2);

[neigborList, neigborCount, edgecellss] = findNeighbors(STATS, colorByNumbers);

%% Create "cells" data structure and store statistics
for i = cellsNumber:-1:1
    cells(i).area = STATS(i).Area / pixelRatio^2; % focus
    cells(i).dpERKIntensity = STATS(i).MeanIntensity;
    cells(i).aspectRatio = STATS(i).MajorAxisLength / STATS(i).MinorAxisLength;
    cells(i).perimeter = STATS(i).Perimeter / pixelRatio;
    cells(i).convexArea = STATS(i).ConvexArea / pixelRatio^2;
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
% everything in "no average" does not get averaged, probably should be
% separate structure
cellFields = fieldnames(cells);

for k = 1:numel(cellFields)
    if ~strcmp(cellFields{k},'noAverage')
        for i = cellsNumber:-1:1
            cells(i).noAverage.neighborAve.(cell2mat(cellFields(k))) = cells(neigborList(i,1)).(cell2mat(cellFields(k)));
            cells(i).noAverage.regionalAve.(cell2mat(cellFields(k))) = cells(neigborList(i,1)).(cell2mat(cellFields(k))) + cells(i).(cell2mat(cellFields(k)));
            for j = neigborList(i,2:neigborCount(i))
                cells(i).noAverage.neighborAve.(cell2mat(cellFields(k))) = cells(i).noAverage.neighborAve.(cell2mat(cellFields(k))) + cells(j).(cell2mat(cellFields(k)));
                cells(i).noAverage.regionalAve.(cell2mat(cellFields(k))) = cells(i).noAverage.regionalAve.(cell2mat(cellFields(k))) + cells(j).(cell2mat(cellFields(k)));
            end
            cells(i).noAverage.neighborAve.(cell2mat(cellFields(k))) = cells(i).noAverage.neighborAve.(cell2mat(cellFields(k))) ./ neigborCount(i);
            cells(i).noAverage.regionalAve.(cell2mat(cellFields(k))) = cells(i).noAverage.regionalAve.(cell2mat(cellFields(k))) ./ (neigborCount(i)+1);
            cells(i).noAverage.relativeValue.(cell2mat(cellFields(k))) = cells(i).(cell2mat(cellFields(k))) ./ cells(i).noAverage.neighborAve.(cell2mat(cellFields(k)));
        end
    end
end

%% Create connection matricies for network statistics
% which cells are connected to which cells - for network statistics
D = pdist(cat(1,Centroid(:,1),Centroid(:,2)));
distance = squareform(D) * pixelRatio;

cMatrixUnweighted = zeros(cellsNumber);
cMatrixWeighted = zeros(cellsNumber);

for i = 1:cellsNumber
    neigbors = neigborList(i,1:neigborCount(i));
    for j = neigbors
        cMatrixUnweighted(i,j) = 1;
        cMatrixWeighted(i,j) = distance(i,j);
    end
end

save([outputfolder 'cMatrixUnweighted.mat'], 'cMatrixUnweighted');
save([outputfolder 'cMatrixWeighted.mat'], 'cMatrixWeighted');

%% Perform Segmentation
% obtain groups of cells manually, should be replaced with Sarah's script
segmentBelongedTo = cell(cellsNumber, 1);
figureHandle = imshow(compositeMap);

x = 1; 
y = 1; 
notDone = true; 
segmentList = segmentListMaster;
identificationArray = cell(1, cellsNumber);

% Replace with Sarah's tool
while notDone && ~isempty(segmentList)
    
    [selection, notDone] = listdlg('ListString',segmentList,'SelectionMode','single','PromptString','Name of Current Segment');
    currentSegment = segmentList{selection};
    segmentList(selection) = [];
    
    if notDone
        [x,y] = ginputc('Color', 'g', 'ShowPoints', true, 'ConnectPoints', true);
        global storedPoints;
        %line([storedPoints(1,1), storedPoints(end,1)], [storedPoints(1,2), storedPoints(end,2)]);
        [ axyp, ~, ~ ] = polygeom( x, y );
        compositeMap = insertText(compositeMap, axyp(2:3), currentSegment, 'FontSize', 18, 'TextColor', ...
            [255 255 255], 'AnchorPoint', 'Center', 'BoxOpacity', 0);
        compositeMap = insertShape(compositeMap, 'Line', [storedPoints(1,1), storedPoints(1,2), storedPoints(end,1), storedPoints(end,2)], 'Color', 'y');
        % insertShape(..., 'FilledPolygon', ...)
        figureHandle = imshow(compositeMap);
        in = inpoly(Centroid, [x y]);
        identificationArray(in) = {currentSegment};
    end
    
end

% any cell that hasn't been identified into a compartment is assigned to
% the closest compartment - maybe not best strategy
unidentified = double(cellfun('isempty',identificationArray));
unidentified = unidentified .* (1:cellsNumber);
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
    cells(i).segmentBelongedTo = identificationArray{i}; % focus
    cells(i).subSegmentBelongedTo = cells(i).segmentBelongedTo;
end

%% Create cell-level macrosegmentation map
for i = length(segmentListMaster):-1:1
    macroSeg(ismember(identificationArray, segmentListMaster{i})) = i;
end

% put dot in each cell and color sections 
createScatterPlots(Centroid(:,1), Centroid(:,2), macroSeg, 'Segmentation Map', '')

%% Run PCA Analysis
cellPCA(cells);

%% Run all possible p tests
% p tests between all properties for every segment
p_dpERK = cellTTest([cells.dpERKIntensity], {cells.segmentBelongedTo});
p_engrailedLevel = cellTTest([cells.engrailedLevel], {cells.segmentBelongedTo});
p_area = cellTTest([cells.area], {cells.segmentBelongedTo});
p_polygonClass = cellTTest([cells.polygonClass], {cells.segmentBelongedTo});

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