function [BestCost, BestSol] = pso(N, Max_Iteration,x,xa,ya,xb,yb,low,up,dim)
%CostFunction=@(x,y,positions) MSError(x,y,positions);       

%% PSO Parameters
nVar=dim;            % Number of Decision Variables
VarSize=[1 nVar];   % Size of Decision Variables Matrix
VarMin = low;         % Lower Bound of Variables
VarMax = up;         % Upper Bound of Variables
MaxIt =  Max_Iteration;  % Maximum Number of Iterations
nPop = N;    			 % Population Size (Swarm Size)

% PSO Parameters
%w=1;            % Inertia Weight
%wdamp=0.99;     % Inertia Weight Damping Ratio
%c1=1.5;         % Personal Learning Coefficient
%c2=2.0;         % Global Learning Coefficient


% If you would like to use Constriction Coefficients for PSO,
% uncomment the following block and comment the above set of parameters.

%% Constriction Coefficients
 phi1=2.05;
 phi2=2.05;
 phi=phi1+phi2;
 chi=2/(phi-2+sqrt(phi^2-4*phi));
 w=chi;          % Inertia Weight
 wdamp=1;        % Inertia Weight Damping Ratio
 c1=chi*phi1;    % Personal Learning Coefficient
 c2=chi*phi2;    % Global Learning Coefficient

% Velocity Limits
VelMax=0.2*(VarMax-VarMin);
VelMin=-VelMax;

%% Initialization
empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];

particle=repmat(empty_particle,nPop,1);

GlobalBest.Cost=inf;

for i=1:nPop
    % Initialize Position
    particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize);
    
    % Evaluation
    particle(i).Cost=Cost(particle(i).Position,x,xa,ya,xb,yb);
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        GlobalBest=particle(i).Best;
    end
end

BestCost=zeros(MaxIt,1);

%% PSO Main Loop
for iter=1:MaxIt
    for i=1:nPop
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Evaluation
        particle(i).Cost = Cost(particle(i).Position,x,xa,ya,xb,yb);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                GlobalBest=particle(i).Best;
            end
        end
    end
    
    BestCost(iter)=GlobalBest.Cost;

    bs = basamakSayisi(iter);
    sifir='';
    for ii = 1:7-bs
        sifir=[sifir '0'];
    end
    fprintf('Iteration %s%d : Best Cost = %15.12e\n', sifir,iter,BestCost(iter));
    
    w=w*wdamp;
end

BestSol = GlobalBest;