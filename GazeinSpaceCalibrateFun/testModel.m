function [SSELT,SSELP,SSERT,SSERP] = testModel(trainStr,testStr)

trainSet = load([pwd '/models/' trainStr 'proc.mat']);
testSet = load([pwd '/models/' testStr 'proc.mat']);

leftEyeX = testSet.data.interpDLAB(:,4);
leftEyeY = testSet.data.interpDLAB(:,5);
rightEyeX = testSet.data.interpDLAB(:,6);
rightEyeY = testSet.data.interpDLAB(:,7);

predAnglesLeft = [predict(trainSet.data.mdlPLT,[leftEyeX leftEyeY]) predict(trainSet.data.mdlPLP,[leftEyeX leftEyeY])];
predAnglesRight = [predict(trainSet.data.mdlPRT,[rightEyeX rightEyeY]) predict(trainSet.data.mdlPRP,[rightEyeX rightEyeY])];

anglesLeft = testSet.data.anglesLeft;
anglesRight = testSet.data.anglesRight;

% check angle prediction
figure
subplot(2,2,1);
plot(predAnglesLeft(:,1),'LineWidth',4); hold on; plot(anglesLeft(:,1),'LineWidth',4);
title('Left Theta Prediction','FontSize',18)
legend('Predicted','Actual');
ylim([min(anglesLeft(:,1)),max(anglesLeft(:,1))]);

subplot(2,2,3);
plot(predAnglesLeft(:,2),'LineWidth',4); hold on; plot(anglesLeft(:,2),'LineWidth',4);
title('Left Phi Prediction','FontSize',18)
legend('Predicted','Actual');
ylim([min(anglesLeft(:,2)),max(anglesLeft(:,2))]);

subplot(2,2,2);
plot(predAnglesRight(:,1),'LineWidth',4); hold on; plot(anglesRight(:,1),'LineWidth',4);
title('Right Theta Prediction','FontSize',18)
legend('Predicted','Actual');
ylim([min(anglesRight(:,1)),max(anglesRight(:,1))]);

subplot(2,2,4);
plot(predAnglesRight(:,2),'LineWidth',4); hold on; plot(anglesRight(:,2),'LineWidth',4);
title('Right Phi Prediction','FontSize',18)
legend('Predicted','Actual');
ylim([min(anglesRight(:,2)),max(anglesRight(:,2))]);

SSELT = sum((predAnglesLeft(:,1) - anglesLeft(:,1)).^2)
SSELP = sum((predAnglesLeft(:,2) - anglesLeft(:,2)).^2)
SSERT = sum((predAnglesRight(:,1) - anglesRight(:,1)).^2)
SSERP = sum((predAnglesRight(:,2) - anglesRight(:,2)).^2)

end