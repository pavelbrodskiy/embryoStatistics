function [ ColorByNumbers ] = labelMaker( Map )

    MapMatrix = double(Map(:,:,1));
    MapMatrix = MapMatrix + double(Map(:,:,2))*1e3;
    MapMatrix = MapMatrix + double(Map(:,:,3))*1e6;
    [a, b] = size(MapMatrix);
    ColorList = reshape(MapMatrix,a*b,1);
    UniqueColorList = unique(ColorList);
    ColorByNumbers = zeros(a, b);
    [CellSize, ~] = size(UniqueColorList);

    for i = 1:CellSize
        ColorByNumbers = ColorByNumbers + i * double(MapMatrix == UniqueColorList(i));
    end

end

