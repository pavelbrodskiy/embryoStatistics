% input scatterX scatterY values title outputFolder
function createScatterPlots(scatterX, scatterY, values, t, outputFolder)

clf;

%figure;

%imshow('Embryo 1 - Posterior Compartment Boundaries.tif');

%hold on

a = jet(256);

[~, newMap] = cmpermute([1,2],a);

colormap(newMap);

%colormap(redgreencmap(256));
%colormap(pmkmp(256));
colormap(jet(256));

%colormap(hot(256));
scatter(scatterX,scatterY,6,values,'fill');
axis ij
axis equal
axis off
title(t)
colorbar('SouthOutside');





print ([outputFolder t],'-r600', '-dtiff');