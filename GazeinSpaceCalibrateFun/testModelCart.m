% function errors = testModelCart(trainStr,testStr)
%
trainStr = 'cal2proc';
testStr = 'cal1proc';

trainSet = load([trainStr '.mat']);
testSet = load([testStr '.mat']);

X = [testSet.resampVICON(:,4:15) testSet.interpDLAB(:,2:5)];

Hd = designfilt('lowpassiir','FilterOrder',14,'HalfPowerFrequency',0.15,'DesignMethod','butter');
predPos = zeros(size(testSet.resampVICON,1),3);
predPosFilt = zeros(size(testSet.resampVICON,1),3);
for i=1:3
    predPos(:,i) = predict(trainSet.mdl{i},X);
    predPosFilt(:,i) = filtfilt(Hd,predPos(:,i));
end

makePredFigCart(testSet.resampVICON(:,1:3),predPos,predPosFilt);

MPredError = calcCartError(testSet.resampVICON(:,1:3),predPos);
MPredFiltError = calcCartError(testSet.resampVICON(:,1:3),predPosFilt);
MSelfError = calcCartError(testSet.resampVICON(:,1:3),testSet.predPos);

% 
% end
