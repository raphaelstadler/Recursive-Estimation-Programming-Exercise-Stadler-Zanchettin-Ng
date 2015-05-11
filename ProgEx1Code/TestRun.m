k = 20;
sumError = 0;
for ii = 1:k
    trackErrorNorm = run(2);
    sumError = sumError + trackErrorNorm;
end
sumError = sumError/k;
disp('ERROR:')
disp(sumError)
disp('====== two =======')
sumError = 0;
for ii = 1:k
    trackErrorNorm = run(2);
    sumError = sumError + trackErrorNorm;
end
sumError = sumError/k;
disp('ERROR:')
disp(sumError)