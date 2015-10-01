function [figH1,figH2,totalMean,totalSD] = makePredictionFigure(predAnglesLeft,anglesLeft,predAnglesRight,anglesRight)
    figH1 = figure;
    ax1=subplot(2,2,1);
    plot(predAnglesLeft(:,1),'LineWidth',4); hold on; plot(anglesLeft(:,1),'.');
    title('Left Theta Prediction','FontSize',18)
    legend('Predicted','Actual');
    ylim([min(anglesLeft(:,1)),max(anglesLeft(:,1))]);

    ax2=subplot(2,2,3);
    plot(predAnglesLeft(:,2),'LineWidth',4); hold on; plot(anglesLeft(:,2),'.');
    title('Left Phi Prediction','FontSize',18)
    legend('Predicted','Actual');
    ylim([min(anglesLeft(:,2)),max(anglesLeft(:,2))]);

    ax3=subplot(2,2,2);
    plot(predAnglesRight(:,1),'LineWidth',4); hold on; plot(anglesRight(:,1),'.');
    title('Right Theta Prediction','FontSize',18)
    legend('Predicted','Actual');
    ylim([min(anglesRight(:,1)),max(anglesRight(:,1))]);

    ax4=subplot(2,2,4);
    plot(predAnglesRight(:,2),'LineWidth',4); hold on; plot(anglesRight(:,2),'.');
    title('Right Phi Prediction','FontSize',18)
    legend('Predicted','Actual');
    ylim([min(anglesRight(:,2)),max(anglesRight(:,2))]);
    
    figH2 = figure;
    
    mseVec(1) = mean((predAnglesLeft(:,1) - anglesLeft(:,1)).^2);
    sdVec(1) = std((predAnglesLeft(:,1) - anglesLeft(:,1)).^2);

    mseVec(2) = mean((predAnglesLeft(:,2) - anglesLeft(:,2)).^2);
    sdVec(2) = std((predAnglesLeft(:,2) - anglesLeft(:,2)).^2);

    mseVec(3) = mean((predAnglesRight(:,1) - anglesRight(:,1)).^2);
    sdVec(3) = std((predAnglesRight(:,1) - anglesRight(:,1)).^2);

    mseVec(4) = mean((predAnglesRight(:,2) - anglesRight(:,2)).^2);
    sdVec(4) = std((predAnglesRight(:,1) - anglesRight(:,1)).^2);
    
    plot(mseVec);
    
    totalMean = sum(mseVec);
    totalSD = sqrt(sum(sdVec.^2));
end