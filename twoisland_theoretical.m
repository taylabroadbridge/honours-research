%Computing the theoretical results for two island branch length
%clear all
%% Parameters
results3 = [];

lambda = 1;

%Migration rates
%M = 0.5;
R = 0.5;

%Coalescent rates
alpha = 1;
beta = 5;

%Initialising vectors (tree heights, branch lengths, correlations)
tree_height=[];
branch_length=[];
corrL1L2 = [];
corrL1Tau = [];
corrL2Tau = [];
for M = 0.2:0.2:5

N = 50; % number of genes


%for alpha=0:0.2:5


%for N=10:10:50
    
%N=10:10:50;
Rewards = [];
Gamma = sparse(sum(N+1:-1:3),sum(N+1:-1:3)); 

p = 0;

for n = N:-1:2

    %% Constructing Lambda (Gamma) matrix 

    Lambda = sparse(n+1,n+1);

    for j = 1:n
       Lambda(j+1,j) = j*R;
       Lambda(j,j+1) = (n-j+1)*M;
    end



    %% Constructing D matrix (coalescent rate matrix)
    if n == 2
        Gamma(1+p:n+1+p,1+p:n+1+p) = Lambda;
    else
        D = sparse(n+1,n);

        for k = 1:n
            if k<n
            D(k,k) = alpha*nchoosek(n-k+1,2); 
            end

            if k > 1
            D(k+1,k) = beta*nchoosek(k,2); %Problem here
            end

        end
        Gamma(1+p:n+1+p,1+p:n+1+p) = Lambda;
        Gamma(1+p:n+1+p,n+2+p:2*n+1+p) = D;
        p = p + n+1;
        
    end
    
end
full(Gamma);
    len = size(Gamma);
    
    %Populating our full sub-intensity matrix B (joining the Gamma and D
    %matrices)
    for i = 1:len(1)
        
        if i == len-2
            
            Gamma(i,i) = -sum(Gamma(i,:))-alpha;
            
        elseif i == len
            
            Gamma(i,i) = -sum(Gamma(i,:))-beta;
            
        else Gamma(i,i) = -sum(Gamma(i,:));
            
        end
    end


full(Gamma);

sum(Gamma');

inverseGamma = inv(full(Gamma));

%Initial distribution vector
pi = zeros(1,length(inverseGamma));
%pi(1)=1;
for y=2:length(pi)
        pi(y)=0;
end

%Half of the population start in I1, half start in I2
pi(N/2+1)=1; 

I1_rewards=[];
I2_rewards=[];
count=0;

%Calculating rewards for each island for each time step
for i=N:-1:2 
I1_rewards = [I1_rewards i:-1:0];
end

for i=N:-1:2 
I2_rewards = [I2_rewards 0:i];
end

%Storing rewards in a matrix
Rewards = [I1_rewards' I2_rewards'];

%Tree height and branch length
tree_height = [tree_height pi*(-inverseGamma)*ones(1,length(pi))'];
branch_length = [branch_length ;pi*(-inverseGamma)*Rewards];

%Cross moments
EL1L2 = pi*(-inverseGamma)*diag(Rewards(:,1))*(-inverseGamma)*Rewards(:,2)+pi*(-inverseGamma)*diag(Rewards(:,2))*(-inverseGamma)*Rewards(:,1);
EL1Tau = pi*(-inverseGamma)*diag(Rewards(:,1))*(-inverseGamma)*ones(1,length(pi))'+pi*(-inverseGamma)*diag(ones(1,length(pi))')*(-inverseGamma)*Rewards(:,1);
EL2Tau = pi*(-inverseGamma)*diag(Rewards(:,2))*(-inverseGamma)*ones(1,length(pi))'+pi*(-inverseGamma)*diag(ones(1,length(pi))')*(-inverseGamma)*Rewards(:,2);

%Expectations (total branch lengths, tree height)
EL1 = pi*(-inverseGamma)*Rewards(:,1);
EL2 = pi*(-inverseGamma)*Rewards(:,2);
ETau = pi*(-inverseGamma)*ones(1,length(pi))';

%Covariances
CovL1L2 = EL1L2-EL1*EL2;
CovL1Tau = EL1Tau-EL1*ETau;
CovL2Tau = EL2Tau-EL2*ETau;

%EL1Sq = 2*pi*(-inverseGamma)*diag(Rewards(:,1))*(-inverseGamma)*Rewards(:,1);%2*pi*(-inverseGamma)*diag(Rewards(:,1))*(-inverseGamma)*Rewards(:,1);
%EL2Sq = 2*pi*(-inverseGamma)*diag(Rewards(:,2))*(-inverseGamma)*Rewards(:,2);%2*pi*(-inverseGamma)*diag(Rewards(:,2))*(-inverseGamma)*Rewards(:,2);

%Second moments for total branch lengths
EL1Sq = 2*pi*(-inverseGamma)*diag(Rewards(:,1))*(-inverseGamma)*Rewards(:,1);
EL2Sq = 2*pi*(-inverseGamma)*diag(Rewards(:,2))*(-inverseGamma)*Rewards(:,2);

%ETauSq = pi*((-inverseGamma)^(2))*ones(1,length(pi))';
%ETauSq = 2*pi*(-inverseGamma)*diag(ones(1,length(pi))')*(-inverseGamma)*ones(1,length(pi))' %+ pi*(-inverseGamma)*diag(ones(1,length(pi))')*(-inverseGamma)*ones(1,length(pi))'
%ETauSq = 2*pi*(-inverseGamma)*ones(1,length(pi))';

%Variances of total branch lengths
varL1 = EL1Sq-EL1^2;
varL2 = EL2Sq-EL2^2;
%varL1 = 2*pi*((full(Gamma))^(-2))*Rewards(:,1)-(pi*(inverseGamma)*Rewards(:,1))^2;
%varL2 = 2*pi*((full(Gamma))^(-2))*Rewards(:,2)-(pi*(inverseGamma)*Rewards(:,2))^2;
%varTau = ETauSq-ETau^2;
varTau = 2*pi*((full(Gamma))^(-2))*ones(1,length(pi))'-(pi*(inverseGamma)*ones(1,length(pi))')^2;
%varTau = 2*pi*(full(Gamma))^(-2)*ones(1,length(pi))'-(pi*inverseGamma*ones(1,length(pi))')^2

%Storing correlations for each iteration
corrL1L2 = [corrL1L2 (CovL1L2)/(sqrt(varL1*varL2))];
corrL1Tau = [corrL1Tau (CovL1Tau)/sqrt(varL1*varTau)];
corrL2Tau = [corrL2Tau (CovL2Tau)/sqrt(varL2*varTau)];

end
%end

