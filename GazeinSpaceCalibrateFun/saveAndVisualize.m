%%%

csvwrite('DlabST1_clean.csv',interpDLAB);
% TIME,ProcPupX,ProcPupY,LeftPupX,LeftPupY,RightPupX,RightPupY

csvwrite('ViconST1_clean.csv',resampVICON);
% MarkX,MarkY,MarkZ,Head1X,Head1Y,Head2Z,Head2X,Head2Y,Head2Z,Head3X,Head3Y,Head3Z


%% 

%visualize the resampling
figure
plot(dataVICON(:,12));
hold on;
plot(resampVICON(:,1));
title('Wand Tip X Resampled to Length');

%visualize the interpolation
figure
subplot(2,2,1);
plot(dataDLAB(:,4)','o');
hold on;
plot(interpDLAB(:,4));
title('Left Eye X');

subplot(2,2,2);
plot(dataDLAB(:,6)','o');
hold on;
plot(interpDLAB(:,6));
title('Right Eye X');

subplot(2,2,3);
plot(dataDLAB(:,5)','o');
hold on;
plot(interpDLAB(:,5));
title('Left Eye Y');

subplot(2,2,4);
plot(dataDLAB(:,7)','o');
hold on;
plot(interpDLAB(:,7));
title('Right Eye Y');

%%

figure;
plot(resampVICON(:,1));
hold on;
plot(resampVICON(:,2));
plot(resampVICON(:,3));
plot(interpDLAB(:,6));
plot(interpDLAB(:,7));