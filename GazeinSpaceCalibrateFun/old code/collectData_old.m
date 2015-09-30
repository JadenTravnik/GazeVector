function [resampVICON,interpDLAB] = collectData(str)

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

% eldata_new.cal.data(1,:) = interpDLAB(:,4); % LEFT X
% eldata_new.cal.data(2,:) = interpDLAB(:,5); % LEFT Y
% eldata_new.cal.data(3,:) = interpDLAB(:,6); % RIGHT X
% eldata_new.cal.data(4,:) = interpDLAB(:,7); % RIGHT Y
% eldata_new.cal.time(1,:) = interpDLAB(:,1); % TIME

% first 300 points after the first 10 will give us static information
% static_range = 10:310;
% eldata_new.static.data(1,:) = interpDLAB(static_range,4); % LEFT X
% eldata_new.static.data(2,:) = interpDLAB(static_range,5); % LEFT Y
% eldata_new.static.data(3,:) = interpDLAB(static_range,6); % RIGHT X
% eldata_new.static.data(4,:) = interpDLAB(static_range,7); % RIGHT Y
% eldata_new.static.time(1,:) = interpDLAB(static_range,1); % TIME

%%

dataVICON = importfileVICON([pwd '/raw/Vicon' str '.csv']);

% indexVec should be the index of the following in order:
% wandX, rForehead, rTemple, lForehead

iV = [24,3,6,9];

sliceVICON = [dataVICON(:,iV(1)) dataVICON(:,iV(1)+1) dataVICON(:,iV(1)+2) ...
    dataVICON(:,iV(2)) dataVICON(:,iV(2)+1) dataVICON(:,iV(2)+2)...
    dataVICON(:,iV(3)) dataVICON(:,iV(3)+1) dataVICON(:,iV(3)+2)...
    dataVICON(:,iV(4)) dataVICON(:,iV(4)+1) dataVICON(:,iV(4)+2)];

resampVICON = zeros(size(dataDLAB,1),12);

for i = 1:size(sliceVICON,2)
    resampVICON(:,i) = resample(sliceVICON(:,i),size(interpDLAB,1),length(sliceVICON(:,1)));
end

% vidata_new.cal.g(1:3,:) = resampVICON(:,1:3)'; % 3D coords of fixation tip
% vidata_new.cal.m1(1:3,:) = resampVICON(:,4:6)'; % 3D coords of marker 1
% vidata_new.cal.m2(1:3,:) = resampVICON(:,7:9)'; % 3D coords of marker 2
% vidata_new.cal.m3(1:3,:) = resampVICON(:,10:12)'; % 3D coords of marker 3
% vidata_new.cal.time = eldata_new.cal.time(1,:); 
% 
% vidata_new.static.g(1:3,:) = resampVICON(static_range,1:3)'; % 3D coords of fixation tip
% vidata_new.static.m1(1:3,:) = resampVICON(static_range,4:6)'; % 3D coords of marker 1
% vidata_new.static.m2(1:3,:) = resampVICON(static_range,7:9)'; % 3D coords of marker 2
% vidata_new.static.m3(1:3,:) = resampVICON(static_range,10:12)'; % 3D coords of marker 3
% vidata_new.static.time = eldata_new.cal.time(1,static_range); 

clear i;

% saveAndVisualize;

end