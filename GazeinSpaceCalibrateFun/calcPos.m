function [predPos,mdl] = calcPos(X,y)
    RegTreeTemp = templateTree('Surrogate','On');
    mdl = fitensemble(X,y,'Bag',120,RegTreeTemp,'type','regression'); 
    predPos = predict(mdl,X);
end