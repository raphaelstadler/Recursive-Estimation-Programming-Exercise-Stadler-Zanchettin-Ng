sumError = 0;
for ii = 1:100
    trackErrorNorm = run(2);
    sumError = sumError + trackErrorNorm;
end
sumError = sumError/100;
disp('ERROR:')
disp(sumError)
k,