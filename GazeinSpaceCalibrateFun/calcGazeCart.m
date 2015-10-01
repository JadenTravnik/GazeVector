function [predPos,predPosFilt,mdl] = calcGazeCart(resampVICON,interpDLAB)

X = [resampVICON(:,4:15) interpDLAB(:,2:5)];

Hd = designfilt('lowpassiir','FilterOrder',14,'HalfPowerFrequency',0.15,'DesignMethod','butter');

predPos = zeros(size(resampVICON,1),3);
predPosFilt = zeros(size(resampVICON,1),3);
mdl = cell(3,1);

for i=1:3
    [predPos(:,i),mdl{i}] = calcPos(X,resampVICON(:,i));
    predPosFilt(:,i) = filtfilt(Hd,predPos(:,i));
end

figure;
ax(1) = subplot(3,1,1);
plot(resampVICON(:,1),'.');
hold on;
plot(predPos(:,1),'LineWidth',2); 
plot(predPosFilt(:,1),'LineWidth',2); 
title('X Position Accuracy');

ax(2) = subplot(3,1,2);
plot(resampVICON(:,2),'.');
hold on;
plot(predPos(:,2),'LineWidth',2); 
plot(predPosFilt(:,2),'LineWidth',2); 
title('Y Position Accuracy');

ax(3) = subplot(3,1,3);
plot(resampVICON(:,3),'.');
hold on;
plot(predPos(:,3),'LineWidth',2); 
plot(predPosFilt(:,3),'LineWidth',2); 
title('Z Position Accuracy');
linkaxes(ax,'x');

% figure; 
% plot3(resampVICON(:,1),resampVICON(:,2),resampVICON(:,3));
% hold on;
% plot3(predPos(:,1),predPos(:,2),predPos(:,3));

end