function error = calcCartError(x1,x2)
    error = sqrt(mean(x1-x2) .^ 2);
end
