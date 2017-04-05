x = linspace(0,50,500)';
y = x+sin(x)+sin(3*x)+0.3*randn(size(x));

gprMdl1 = fitrgp(x,y,'KernelFunction','squaredexponential','Sigma',sigma0);

sigma0 = 0.1;
kparams0 = [3, 6];
gprMdl2 = fitrgp(x,y,'KernelFunction','squaredexponential',...
     'KernelParameters',kparams0,'Sigma',sigma0);
 
ypred1 = resubPredict(gprMdl1);
ypred2 = resubPredict(gprMdl2);

figure();
plot(x,y,'r.');
hold on
plot(x,ypred1,'b');
% plot(x,ypred2,'g');
xlabel('x');
ylabel('y');
% legend({'data','default kernel parameters',...
% 'kparams0 = [3.5,6.2], sigma0 = 0.2'},...
% 'Location','Best');
title('Impact of initial kernel parameter values');
hold off
