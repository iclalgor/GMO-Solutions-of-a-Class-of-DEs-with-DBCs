%
% Copyright (c) 2016, Korhan Gunel
% All rights reserved.
%
% Project Code: OYP-14xxx
% Project Title: Numerical Solutions of Differential Equations with Neural
% Nets
%
% Developer: Korhan Gunel (Adnan Menderes University, Department of Mathematics)
%
% Contact Info: kgunel@adu.edu.tr
%
tic;
clear; close; clc

exampleNo = 2;
xa=0;     ya = 1; 
xb=1;  yb = 0;
h=0.01;
x=xa:h:xb;
y=exact_solution(x);

% PSO parameters
N = 50;                        % # of particle
Max_Iteration  = 100;         % Maximum number of "iterations"
low = -1; % Lower bound
up  = 1;  % Upper bound
m = 6;     % # of neurons in NN
dim = 3*m;   % Size of the particle
PSO_parameters = [N m low up dim Max_Iteration];

for k=1:10
    T = tic;
    [BestCost, BestSol] = pso(N, Max_Iteration,x,xa,ya,xb,yb,low,up,dim);
    Ttrain(k) = toc(T);
    PSO_Solution(k).BestCost=BestCost;
    PSO_Solution(k).BestSol=BestSol;
    if k==1
        bestof_BestCost = PSO_Solution(k).BestCost;
        bestof_BestSol = PSO_Solution(k).BestSol;
    end
    if bestof_BestSol.Cost>PSO_Solution(k).BestSol.Cost
        bestof_BestCost = PSO_Solution(k).BestCost;
        bestof_BestSol = PSO_Solution(k).BestSol;
    end
    
    cf = PSO_Solution(k).BestSol.Position;
    alpha = cf(1:m);
    w = cf(m+1:2*m);
    bias = cf(2*m+1:3*m);
    
    yt = trialy(x,xa,ya,xb,yb,cf); % Trial solutions using GSA optimization for training set
    xt=xa:h/2:xb; % Test Set
    % Remove elements from xe occurred in training set
    xt = xt(~ismember(xt,x));
    % Insert boundary points to the test set
    xt = [xa xt xb];
    T = tic;
    ytt = trialy(xt,xa,ya,xb,yb,cf); % Trial solutions using GSA optimization for test set
    Ttest(k) = toc(T); 
    
    trainMSE(k) = sum((exact_solution(x)-yt).^2)/length(yt);
    testMSE(k) = sum((exact_solution(xt)-ytt).^2)/length(ytt);
    PSO_Solution(k).train.MSE =   trainMSE(k);
    PSO_Solution(k).test.MSE =  testMSE(k);
end

%% Results
alpha = bestof_BestSol.Position(1:m);
weight = bestof_BestSol.Position(m+1:2*m);
bias = bestof_BestSol.Position(2*m+1:3*m);

% Plot cost function
semilogy(bestof_BestCost,'-r','LineWidth',1.25);
title('\fontsize{12}\bf Maliyet Fonksiyonu');
xlabel('\fontsize{12}\bf İterasyon');
ylabel('\fontsize{12}\bf Uygunluk değeri (En iyi değer)');
legend('\fontsize{10}\bf PSO');
fig=gcf;
fig.InvertHardcopy = 'on';
saveas(gcf,['Figure_Exmp_' num2str(exampleNo) '_b.fig']);
print(gcf,['Figure_Exmp_' num2str(exampleNo) '_b.jpg'],'-djpeg','-r300');
print(gcf,['Figure_Exmp_' num2str(exampleNo) '_b.eps'],'-depsc','-r300');

xe = xa:h/10:xb;
ye = exact_solution(xe); % Exact Solution
yt = trialy(x,xa,ya,xb,yb,bestof_BestSol.Position); % Trial solutions using GSA optimization for training set
xt=xa:h/2:xb; % Test Set
% Remove elements from xe occurred in training set
xt = xt(~ismember(xt,x));
% Insert boundary points to the test set
xt = [xa xt xb];
ytt = trialy(xt,xa,ya,xb,yb,bestof_BestSol.Position); % Trial solutions using GSA optimization for test set
figure;
hold on;
plot(x,yt,'ro',xt,ytt,'-.b',xe,ye,'-k','LineWidth',1.25); % Exact Solution for Training Set
title('\fontsize{12}\bf PSO ile optimize edilen Yapay Sinir Ağı Çözümleri');

xlabel('\fontsize{12}\bf x');
ylabel('\fontsize{12}\bf y');
%legend('\fontsize{10}\bf Numerical solution for training set','\fontsize{10}\bf Numerical solution for test set','\fontsize{10}\bf Exact Solution','Location','southeast');
legend('\fontsize{10}\bf Eğitim kümesinin sayısal çözümü','\fontsize{10}\bf Test kümesinin sayısal çözümü','\fontsize{10}\bf Gerçek çözüm','Location','southwest');
hold off;

fig=gcf;
fig.InvertHardcopy = 'on';
saveas(fig,['Figure_Exmp_' num2str(exampleNo) '_a.fig']);
print(gcf,['Figure_Exmp_' num2str(exampleNo) '_a.jpg'],'-djpeg','-r300');
print(gcf,['Figure_Exmp_' num2str(exampleNo) '_a.eps'],'-depsc','-r300');

title = ['The numerical solution of Feed-forward Neural Network trained by PSO for training set in Example ' num2str(exampleNo)];
label = ['lbl:tabloExmp' num2str(exampleNo) '_train'];
displayResults(x,exact_solution(x),yt,title,label,exampleNo,PSO_parameters);

title = ['The numerical solution of Feed-forward Neural Network trained by PSO for test set in Example ' num2str(exampleNo)];
label = ['lbl:tabloExmp' num2str(exampleNo) '_test'];
displayResults(xt,exact_solution(xt),ytt,title,label,exampleNo,PSO_parameters);
clear fig;
save all
fprintf('Mean of MSEs for training set %1.3e ± %1.3e\n',mean(trainMSE),std(trainMSE));
fprintf('Mean of MSEs for test set %1.3e ± %1.3e\n',mean(testMSE),std(testMSE));

wtime = toc;
fprintf ( 1, '  Elapsed time %f seconds to run.\n', wtime );

fprintf('Mean of elapsed time in seconds for training set %1.3e ± %1.3e\n',mean(Ttrain),std(Ttrain));
fprintf('Mean of elapsed time in seconds for test set %1.3e ± %1.3e\n',mean(Ttest),std(Ttest));
