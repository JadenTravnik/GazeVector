collectData

load data_cal.mat

eldata_BU = eldata;
vidata_BU = vidata;

eldata = eldata_new;
vidata = vidata_new;

%% 

% [par] = GazeinSpaceCalibrateFun(eldata, vidata,opt);