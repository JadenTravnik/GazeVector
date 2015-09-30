function [] = buildGroupModel(names)

allVICON = [];
allDLAB = [];

for i = 1:length(names)
    str = names{i};
    [data.resampVICON,data.interpDLAB] = collectData(names{i});
    allVICON = [allVICON;data.resampVICON];
    allDLAB = [allDLAB;data.interpDLAB];
end

[data.headAngles, data.anglesLeft, data.anglesRight, data.predAnglesLeft, data.predAnglesRight, ...
    data.SSELT, data.SSELP, data.SSERT, data.SSERP, data.mdlPLT, data.mdlPLP, data.mdlPRT, data.mdlPRP]...
    = calcGaze(allVICON,allDLAB);

% save([pwd '/models/group' regexprep(strjoin(names),'[^\w'']','') 'proc.mat'],'-v7.3');

end