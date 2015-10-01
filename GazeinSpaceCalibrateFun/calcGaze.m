function allData = calcGaze(resampVICON,interpDLAB)

calibrationPoints = resampVICON(:,1:3);

headAngles = zeros(size(resampVICON,1),2);
anglesLeft = zeros(size(resampVICON,1),2);
anglesRight = zeros(size(resampVICON,1),2);

rForehead = zeros(size(resampVICON,1),3);
rTemple = zeros(size(resampVICON,1),3);
lForehead = zeros(size(resampVICON,1),3);
lTemple = zeros(size(resampVICON,1),3);

eyePosLeft = zeros(size(resampVICON,1),3);
eyePosRight = zeros(size(resampVICON,1),3);

for i=1:size(resampVICON,1);
    rForehead(i,:) = resampVICON(i,4:6);
    rTemple(i,:) = resampVICON(i,7:9);
    lForehead(i,:) = resampVICON(i,10:12);
    lTemple(i,:) = resampVICON(i,13:15);
    
    v1 = cross(lForehead(i,:)-rForehead(i,:), lForehead(i,:)-rTemple(i,:));
    v2 = cross(lForehead(i,:)-rForehead(i,:),v1);
    perpVector = cross(v1,v2);
    headDirection = cross(v1,perpVector);
    headDirection = headDirection / norm(headDirection); % just to make it unit length
    [azimuth,elevation,~] = cart2sph(headDirection(1),headDirection(2),headDirection(3));
    headAngles(i,1:2) = [azimuth elevation];
    
    % The measured eye offset in mm
    eyeOffsetLeft = [50, 60, -30];
    eyeOffsetRight = [50, 60, -88];
    
    % in world coordinates
    eyePosLeft(i,:) = calcEyePosition(lForehead(i,:), rForehead(i,:), rTemple(i,:), eyeOffsetLeft);
	eyePosRight(i,:) = calcEyePosition(lForehead(i,:), rForehead(i,:), rTemple(i,:), eyeOffsetRight);
    
%     pass in headAngles rather than recalculating in the function
    anglesLeft(i,:) = calcAngles(eyePosLeft(i,:), headAngles(i,:), calibrationPoints(i,:));
    anglesRight(i,:) = calcAngles(eyePosRight(i,:), headAngles(i,:), calibrationPoints(i,:));
end

leftEyeX = interpDLAB(:,2);
leftEyeY = interpDLAB(:,3);
rightEyeX = interpDLAB(:,4);
rightEyeY = interpDLAB(:,5);

diffSumLeft = abs(diff(leftEyeX)) + abs(diff(leftEyeY));
diffSumRight = abs(diff(rightEyeX)) + abs(diff(rightEyeY));

makeEyeAnglePupilFigure(anglesLeft,leftEyeX,leftEyeY,anglesRight,rightEyeX,rightEyeY);

% less than the mean of the given signal
% MAGIC constant
stabilityConst = mean(diffSumLeft)+mean(diffSumRight);
stableIDX = (diffSumLeft<stabilityConst)+(diffSumRight<stabilityConst)>-1;

[pLT,mdlPLT] = modelAndPredictAngle(anglesLeft(:,1),[leftEyeX leftEyeY rightEyeX rightEyeY],stableIDX);

% [pLTless,mdlPLTless] = modelAndPredictAngle(anglesLeft(:,1),[leftEyeX leftEyeY],stableIDX);
% figure; 
% plot(oobLoss(mdlPLT,'mode','cumulative'),'k--'); 
% hold on; 
% plot(oobLoss(mdlPLTless,'mode','cumulative'),'r.');
% title('Relative Error on Bagged Trees for single or both eye data');
% legend('left and right','only left');

[pLP,mdlPLP] = modelAndPredictAngle(anglesLeft(:,2),[leftEyeX leftEyeY rightEyeX rightEyeY],stableIDX);
[pRT,mdlPRT] = modelAndPredictAngle(anglesRight(:,1),[leftEyeX leftEyeY rightEyeX rightEyeY],stableIDX);
[pRP,mdlPRP] = modelAndPredictAngle(anglesRight(:,2),[leftEyeX leftEyeY rightEyeX rightEyeY],stableIDX);

predAnglesLeft = [pLT pLP];
predAnglesRight = [pRT pRP];

% Filter the data to remove prediction noise
Hd = designfilt('lowpassiir','FilterOrder',14,'HalfPowerFrequency',0.15,'DesignMethod','butter');
     
predAnglesLeft = filtfilt(Hd,predAnglesLeft);
predAnglesRight = filtfilt(Hd,predAnglesRight);

% % check angle prediction
[~,~,tM,tSTD] = makePredictionFigure(predAnglesLeft,anglesLeft,predAnglesRight,anglesRight);

leftVector = zeros(size(resampVICON,1),3);
rightVector = zeros(size(resampVICON,1),3);

for i=1:size(resampVICON,1)
    [leftX, leftY, leftZ] = sph2cart(predAnglesLeft(i,1)+headAngles(i,1),predAnglesLeft(i,2)+headAngles(i,1),1);
    [rightX, rightY, rightZ] = sph2cart(predAnglesRight(i,1)+headAngles(i,1),predAnglesRight(i,2)+headAngles(i,1),1);

%     focusPoints(i,:) = calcMidPoint(eyePosLeft(i,:), eyePosRight(i,:), [leftX, leftY, leftZ], [rightX, rightY, rightZ]);

    leftTemp = [leftX, leftY, leftZ];
    rightTemp = [rightX, rightY, rightZ];
    
    leftVector(i,:) = leftTemp + eyePosLeft(i,:);
    rightVector(i,:) = rightTemp + eyePosRight(i,:);
end

allData.headAngles = headAngles;
allData.anglesLeft = anglesLeft;
allData.anglesRight = anglesRight;
allData.predAnglesLeft = predAnglesLeft;
allData.predAnglesRight = predAnglesRight;
allData.mdlPLT = mdlPLT;
allData.mdlPLP = mdlPLP;
allData.mdlPRT = mdlPRT;
allData.mdlPRP = mdlPRP;
allData.tM = tM;
allData.tSTD = tSTD;

%% Supplementary Functions

function eyePosTemp = calcEyePosition(lForehead, rForehead, rTemple, eyeOffset)
    v1 = cross(lForehead-rForehead, lForehead-rTemple);
    v2 = cross(lForehead-rForehead, v1);
    A = [v1 / norm(v1); v2/norm(v2); (lForehead-rForehead)/norm(lForehead-rForehead)]';
    eyePosTemp = A*eyeOffset' + lForehead';
end

function angles = calcAngles(eyePos, headAngles, target) 
    eyeTargetVec = target-eyePos;
    [az,el,~] = cart2sph(eyeTargetVec(1),eyeTargetVec(2),eyeTargetVec(3));
    angles = [az - headAngles(1), el - headAngles(2)];
end

function vecN = changeRange(oldMin, oldMax, newMin, newMax, vec)
	oldRange = (oldMax - oldMin);  
	newRange = (newMax - newMin); 
	vecN = (((vec - oldMin) * newRange) / oldRange) + newMin;
end

function [predAngles,mdl] = modelAndPredictAngle(y,X,idx)
    stableX = X(idx,:);
    stabley = y(idx);
    RegTreeTemp = templateTree('Surrogate','On');
    mdl = fitensemble(stableX,stabley,'Bag',60,RegTreeTemp,'type','regression'); 
    predAngles = predict(mdl,X);
end

function figH = makeEyeAnglePupilFigure(anglesLeft,leftEyeX,leftEyeY,anglesRight,rightEyeX,rightEyeY)
    figH = figure; 
    subplot(2,1,1);
    plot(anglesLeft,'LineWidth',4); 
    hold on; 
    plot(-changeRange(0, 384, -2*pi, 2*pi, leftEyeX),'LineWidth',4);
    plot(-changeRange(0, 288, -2*pi, 2*pi, leftEyeY),'LineWidth',4);
    h_legend=legend('LeftEyeTheta','LeftEyePhi','LeftEyeX','LeftEyeY');
    set(h_legend,'FontSize',14,'Location','southeast');
    title('Left Eye Theta/Phi with Left Pupil X/Y') 

    subplot(2,1,2);
    plot(anglesRight,'LineWidth',4); 
    hold on; 
    plot(-changeRange(0, 384, -2*pi, 2*pi, rightEyeX),'LineWidth',4);
    plot(-changeRange(0, 288, -2*pi, 2*pi, rightEyeY),'LineWidth',4);
    h_legend=legend('RightEyeTheta','RightEyePhi','RightEyeX','RightEyeY');
    set(h_legend,'FontSize',14,'Location','southeast');
    title('Right Eye Theta/Phi with Right Pupil X/Y');
end




end