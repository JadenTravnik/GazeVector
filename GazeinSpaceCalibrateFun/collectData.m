function [resampVICON,interpDLAB] = collectData(str)

% headersDLAB = importfileDlabHeader([pwd '/raw/Dlab' str '.csv']);
% 
% try
%     iVDLAB = [find(strcmp(headersDLAB,'UTC')),...
%         find(strcmp(headersDLAB,'Dikablis Professional_Eye Data_Original_Left Eye_Pupil X')),...
%         find(strcmp(headersDLAB,'Dikablis Professional_Eye Data_Original_Left Eye_Pupil Y')),...
%         find(strcmp(headersDLAB,'Dikablis Professional_Eye Data_Original_Right Eye_Pupil X')),...
%         find(strcmp(headersDLAB,'Dikablis Professional_Eye Data_Original_Right Eye_Pupil Y'))];
% catch ME
%     ME
% end

dataDLAB = importfileDLAB([pwd '/raw/Dlab' str '.csv']);

interpDLAB = zeros(size(dataDLAB));
interpDLAB(:,1) = dataDLAB(:,1);

for i = 2:size(dataDLAB,2)
    X = dataDLAB(:,i);
    X(X==0) = NaN;
    X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)),'pchip');
    interpDLAB(:,i) = X;
end

clear X i;

%%

headersVICON = importfileViconHeader([pwd '/raw/Vicon' str '.csv']);
dataVICON = importfileVICON([pwd '/raw/Vicon' str '.csv']);

% indexVec should be the index of the following in order:
% wandX, rForehead, rTemple, lForehead, lTemple

if (find(strcmp(headersVICON,'Wand:Tip')))
    iV = [find(strcmp(headersVICON,'Wand:Tip')),...
    find(strcmp(headersVICON,'Dikablis:Head1')),...
    find(strcmp(headersVICON,'Dikablis:Head2')),...
    find(strcmp(headersVICON,'Dikablis:Head3')),...
    find(strcmp(headersVICON,'Dikablis:Head4'))];
else 
    iV = [find(strcmp(headersVICON,'CalWand:Tip')),...
    find(strcmp(headersVICON,'Head:RFRNT')),...
    find(strcmp(headersVICON,'Head:RBCK')),...
    find(strcmp(headersVICON,'Head:LFRNT')),...
    find(strcmp(headersVICON,'Head:LBCK'))];
end
    
sliceVICON = [dataVICON(:,iV(1)) dataVICON(:,iV(1)+1) dataVICON(:,iV(1)+2) ...
            dataVICON(:,iV(2)) dataVICON(:,iV(2)+1) dataVICON(:,iV(2)+2)...
            dataVICON(:,iV(3)) dataVICON(:,iV(3)+1) dataVICON(:,iV(3)+2)...
            dataVICON(:,iV(4)) dataVICON(:,iV(4)+1) dataVICON(:,iV(4)+2)...
            dataVICON(:,iV(5)) dataVICON(:,iV(5)+1) dataVICON(:,iV(5)+2)];

resampVICON = zeros(size(dataDLAB,1),12);

for i = 1:size(sliceVICON,2)
    resampVICON(:,i) = resample(sliceVICON(:,i),size(interpDLAB,1),length(sliceVICON(:,1)));
end

for i = 1:size(resampVICON,2)
    X = resampVICON(:,i);
    X(X==0) = NaN;
    X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)),'pchip');
    resampVICON(:,i) = X;
end

clear i;

% saveAndVisualize;

end