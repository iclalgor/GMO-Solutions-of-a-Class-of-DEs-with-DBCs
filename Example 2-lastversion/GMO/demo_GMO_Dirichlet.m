tic;
clear; close;clc;
exampleNo = 2;
%% Problem
% y'' = f(x,y,y'), a<=x<=b
% y(a) = A;
% y(b) = B;
a = 0; %a-xa
b = 1; %b-xb
A = 1;  %A-ya
B = 0;  %B-yb
h = 0.01;
x = a:h:b;
indisler = randperm(length(x));
%batchSize = 20;
%nBatch = floor(length(x)/batchSize);
%X = zeros(nBatch,batchSize);
%for i=1:nBatch
%    X(i,:) = x( indisler(20*(i-1)+1:20*i));
%end

%% Model
noNeurons = 10;
params.noNeurons = noNeurons;
model = fitnet(noNeurons);
%model = fitnet(noNeurons,'traingd');
model = init(model);
% Pretraining for boundaries
model = train(model,[a b], [A B]);
W = getwb(model);

problem.model =  model;
problem.x = x;
problem.a = a;
problem.b = b;
problem.A = A;
problem.B = B;

% Train model with GMO
run=10;     % Maximum number of the algorithm runnings conducted
np=50;     % Number of search agents (solutions)
maxit=3;  % Maximum number of iterations
epsilon=0; % This parameter must be either set to "eps" (typically for the uni-modal or simple functions) or zero (typically for the multi-modal or complex functions)
lb = -1;
ub = 1;
nx = numel(W); 
varmax=ub*ones(1,nx); % Upper bound defined for the positions which can generally be a desired vector
varmin=lb*ones(1,nx); % Lower bound defined for the positions which can generally be a desired vector
limvel=0.1; % A ratio of the maximum distance in the search space to form the maximum velocity 
velmax=limvel*(varmax(1,1:nx)-varmin(1,1:nx)); % Upper bound defined for the velocities
velmin=-velmax; % Lower bound defined for the velocities
% Definitions of GMO parameters
params.np = np;
params.nx = nx;
params.maxit = maxit;
params.varmax = varmax ;
params.varmin = varmin ;
params.velmax = velmax;
params.velmin = velmin ;
params.epsilon = 1;

z_iter_main=zeros(run,maxit); % her iterasyonda elde edilen
z_final_main=zeros(run,1);
pos_final_main=zeros(run,nx);
%x1=zeros(maxit);
%y1=zeros(maxit);


% for nrun=1:run
%     for i=1:nBatch
%         T = tic;
%         problem.x = X(i,:);
%         [z_iter,z_final,pos_final] = GMO(problem,params);
%         problem.model = setwb(problem.model,pos_final);
%         Ttrain(nrun) = toc(T);
%         GMO_Solution(nrun).z_iter=z_iter;
%         GMO_Solution(nrun).z_final=z_final;
%         GMO_Solution(nrun).pos_final=pos_final;
%     end
%     z_iter_main(nrun,1:maxit)=z_iter(1:maxit);
%     z_final_main(nrun)=z_final;
%     pos_final_main(nrun,1:nx)=pos_final(1:nx);
%     disp(['The best objective function value obtained by GMO = ',num2str(z_final_main(nrun))]);
%     disp(['The best solution obtained by GMO = ','[',num2str(pos_final_main(nrun,1:nx)),']']);
% end

for k=1:10
        %for i=1:nBatch
        T = tic;
        %problem.x = X(i,:); % Training set input
        [z_iter,z_final,pos_final] = GMO(problem,params);
        problem.model = setwb(problem.model,pos_final);
        Ttrain(k) = toc(T);
        GMO_Solution(k).z_iter = z_iter;
        GMO_Solution(k).z_final = z_final;
        GMO_Solution(k).pos_final = pos_final;
        if k==1
            bestof_z_final = GMO_Solution(k).z_final;
            bestof_pos_final = GMO_Solution(k).pos_final;
        end
        if bestof_z_final > GMO_Solution(k).z_final
            bestof_z_final = GMO_Solution(k).z_final;
            bestof_pos_final = GMO_Solution(k).pos_final;
        end
        %end

    % Set trained model with best weights
    model = setwb(model, bestof_pos_final);

   % Generate trial solutions for training set
    yt = trialSolution(model, x,a,b,A,B);% Use trained model for training set
    % Generate trial solutions for test set
    xt = a:h/10:b;%a + (b-a)*rand(1,20); % Test set
    xt = xt(~ismember(xt, x)); % Remove training set points
    xt = [a xt b]; % Insert boundary points
    T=tic;
    ytt = trialSolution(model, xt, a, b, A, B); % Use trained model for test set
    Ttest(k)=toc(T);

    % Calculate mean squared errors
    trainMSE(k) = sum((exactSolution(x) - yt).^2) / length(yt);
    testMSE(k) = sum((exactSolution(xt) - ytt).^2) / length(ytt);
    GMO_Solution(k).train.MSE = trainMSE(k);
    GMO_Solution(k).test.MSE = testMSE(k);
end


%% Results
% Set model with the best weights
[minLoss, minId] = min(z_final_main);
bestWB = pos_final_main(minId, :);
model = setwb(model, bestWB);

% Plot cost function
semilogy(z_iter, '-r', 'LineWidth', 1.25);
title('\fontsize{12}\bf Maliyet Fonksiyonu');
xlabel('\fontsize{12}\bf İterasyon');
ylabel('\fontsize{12}\bf Uygunluk değeri (En iyi değer)');
legend('\fontsize{10}\bf GMO');
% fig = gcf;
% fig.InvertHardcopy = 'on';
% saveas(gcf, ['Figure_Exmp_' num2str(exampleNo) '_b.fig']);
print(gcf, ['Figure_Exmp_' num2str(exampleNo) '_b.jpg'], '-djpeg', '-r300');
print(gcf, ['Figure_Exmp_' num2str(exampleNo) '_b.eps'], '-depsc', '-r300');

% Plot exact solution
xe = a:h:b;
ye = exactSolution(xe); % Exact Solution
figure;
hold on;
plot(xe, ye, '-k', 'LineWidth', 1.25); % Exact Solution for Training Set

% Plot training set solutions
yt = trialSolution(model, x, a, b, A, B);
plot(x, yt, 'ro', 'LineWidth', 1.25); % Numerical solution for training set

% Plot test set solutions
xtest = a:h/2:b;
% Remove elements from xe occurred in training set
xtest = xtest(~ismember(xtest,x));
% Insert boundary points to the test set
xtest = [a xtest b];
ytest = trialSolution(model, xtest, a, b, A, B);
plot(xtest, ytest, 'gx', 'LineWidth', 1.25); % Numerical solution for test set

title('\fontsize{12}\bf GMO ile optimize edilen Yapay Sinir Ağı Çözümleri');
xlabel('\fontsize{12}\bf x');
ylabel('\fontsize{12}\bf y');
legend('Gerçek çözüm', '\fontsize{10}\bf Eğitim kümesinin sayısal çözümü', '\fontsize{10}\bf Test kümesinin sayısal çözümü', 'Location', 'southeast');
hold off;

% fig = gcf;
% fig.InvertHardcopy = 'on';
% saveas(fig, ['Figure_Exmp_' num2str(exampleNo) '_a.fig']);
% print(gcf, ['Figure_Exmp_' num2str(exampleNo) '_a.jpg'], '-djpeg', '-r300');
% print(gcf, ['Figure_Exmp_' num2str(exampleNo) '_a.eps'], '-depsc', '-r300');

%% Reports
title = ['The numerical solution of Feed-forward Neural Network trained by PSO for training set in Example ' num2str(exampleNo)];
label = ['lbl:tabloExmp' num2str(exampleNo) '_train'];
displayResults(x,exactSolution(x),yt,title,label,exampleNo,params);

title = ['The numerical solution of Feed-forward Neural Network trained by PSO for test set in Example ' num2str(exampleNo)];
label = ['lbl:tabloExmp' num2str(exampleNo) '_test'];
displayResults(xtest,exactSolution(xtest),ytest,title,label,exampleNo,params);
clear fig;
save all
fprintf('Mean of MSEs for training set %1.3e ± %1.3e\n',mean(trainMSE),std(trainMSE));
fprintf('Mean of MSEs for test set %1.3e ± %1.3e\n',mean(testMSE),std(testMSE));

wtime = toc;
fprintf(1, '  Elapsed time %f seconds to run.\n', wtime);

fprintf('Mean of elapsed time in seconds for training set %1.3e ± %1.3e\n', mean(Ttrain), std(Ttrain));
fprintf('Mean of elapsed time in seconds for test set %1.3e ± %1.3e\n', mean(Ttest), std(Ttest));
