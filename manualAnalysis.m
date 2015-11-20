for i = 1:3183
    switch cells(i).segmentBelongedTo
        case 'Arf'
            AMANUALidentity(i) = 1;
        case 'Arb'
            AMANUALidentity(i) = 2;
        case 'Pf'
            AMANUALidentity(i) = 3;
        case 'Pb'
            AMANUALidentity(i) = 4;
        otherwise
            AMANUALidentity(i) = 0;
    end
end

AfAnis = [cells(AMANUALidentity==1).aspectRatio];
AbAnis = [cells(AMANUALidentity==2).aspectRatio];
PfAnis = [cells(AMANUALidentity==3).aspectRatio];
PbAnis = [cells(AMANUALidentity==4).aspectRatio];

AfAnis = [cells(AMANUALidentity==1).dpERKIntensity];
AbAnis = [cells(AMANUALidentity==2).dpERKIntensity];
PfAnis = [cells(AMANUALidentity==3).dpERKIntensity];
PbAnis = [cells(AMANUALidentity==4).dpERKIntensity];

AfAnisAve = [];
AbAnisAve = [];
PfAnisAve = [];
PbAnisAve = [];


for i = 1:3183
    switch cells(i).segmentBelongedTo
        case 'Arf'
            AfAnisAve = [AfAnisAve cells(i).noAverage.regionalAve.aspectRatio];
        case 'Arb'
            AbAnisAve = [AbAnisAve cells(i).noAverage.regionalAve.aspectRatio];
        case 'Pf'
            PfAnisAve = [PfAnisAve cells(i).noAverage.regionalAve.aspectRatio];
        case 'Pb'
            PbAnisAve = [PbAnisAve cells(i).noAverage.regionalAve.aspectRatio];
    end
end