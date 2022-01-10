%Computing the theoretical results for the seed-bank coalescent
%% Parameters
results3 = [];
lambda = 1;

%Activation and dormancy rates
c = 1;
K = 1;

%Only active population can coalesce
alpha = 1;
beta = 0;

%Initialising vectors (tree heights, branch lengths, correlations)
tree_height=[];
branch_length=[];


%for N=10:10:50
N=10;   
Rewards = [];
Gamma = sparse(sum(N+1:-1:3),sum(N+1:-1:3)); 
p = 0;

for n = N:-1:2

    %% Constructing Lambda (Gamma) matrix 

    Lambda = sparse(n+1,n+1);

    for j = 1:n
       Lambda(j+1,j) = j*K;
       Lambda(j,j+1) = (n-j+1)*c;
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
for y=2:length(pi)
        pi(y)=0;
end

% If all genes are active at time 0
pi(1)=1;


%Storing rewards in a matrix
%Rewards = [I1_rewards' I2_rewards'];
Rewards1=[];

for i = N:-1:2
Rewards1 = [Rewards i:-1:0];
end

%Tree height and branch lengths
tree_height = [tree_height pi*(-inverseGamma)*ones(1,length(pi))'];
branch_length = [branch_length, pi*(-inverseGamma)*Rewards'];

%end
%end
