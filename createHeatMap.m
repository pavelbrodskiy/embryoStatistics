% input scatterX scatterY values title outputFolder
function createHeatMap(colorByNumbers, values, t)

clf;
HeatMap = zeros(size(colorByNumbers));

for i = 1:(max(colorByNumbers(:))-1)
    HeatMap(colorByNumbers==i) = values(i);
end

imwrite(uint16(HeatMap),t);