function fig = makePredFigCart(V,pP,pPF)

fig = figure;
ax(1) = subplot(3,1,1);
plot(V(:,1),'.');
hold on;
plot(pP(:,1),'LineWidth',2); 
plot(pPF(:,1),'LineWidth',2); 
title('X Position Accuracy');
legend('Actual','Predicted','FiltPred');

ax(2) = subplot(3,1,2);
plot(V(:,2),'.');
hold on;
plot(pP(:,2),'LineWidth',2); 
plot(pPF(:,2),'LineWidth',2); 
title('Y Position Accuracy');
legend('Actual','Predicted','FiltPred');

ax(3) = subplot(3,1,3);
plot(V(:,3),'.');
hold on;
plot(pP(:,3),'LineWidth',2); 
plot(pPF(:,3),'LineWidth',2); 
title('Z Position Accuracy');
linkaxes(ax,'x');
legend('Actual','Predicted','FiltPred');

% figure; 
% plot3(V(:,1),V(:,2),V(:,3));
% hold on;
% plot3(pP(:,1),pP(:,2),pP(:,3));

end