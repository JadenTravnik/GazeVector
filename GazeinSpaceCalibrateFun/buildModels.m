function [] = buildModels(str)

[data.resampVICON,data.interpDLAB] = collectData(str);

[data.headAngles, data.anglesLeft, data.anglesRight, data.predAnglesLeft, data.predAnglesRight, ...
    data.SSELT, data.SSELP, data.SSERT, data.SSERP, data.mdlPLT, data.mdlPLP, data.mdlPRT, data.mdlPRP]...
    = calcGaze(data.resampVICON,data.interpDLAB);

save([pwd '/models/' str 'proc.mat']);

end