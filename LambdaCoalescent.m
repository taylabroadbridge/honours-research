%Code to calculate (and plot) theoretical results for Lambda coalescent
%Can be adapted for both psi and beta coalescent

%Initialising
hts = []; %tree heights
lengths = []; %total branch lengths
hts_second_moment = []; %second moment for tree height
lengths_second_moment = []; %second moment for branch length
Pi = []; %Initial distribution vector

%Parameters: psi and beta (depending which coalescent we're looking at)
%alpha=1;
psi = 0.5; 

%Initial population size (number of genes) 
%nn = 20; 


for trials = 100:50:200 %Repeat for initial population sizes 10 to 100
for nn=1:trials 

%Initialising sub-intensity matrix
A = zeros(nn-1,nn-1);

%Reward vector
R = [nn:-1:2]';
%R = diag(Rewards);

%Constructing the sub-intensity matrix given our parameters
for i = 1:nn-1
    
    for j = 1:nn-1-i
        
       A(i,i+j) = g(nn-i+1,j+1,psi);
         
    end
    A(i,i) = -sum(A(i,:)); %diagonal elements
    Pi=[1 zeros(1, length(A)-1)]; %initial distribution vector
    
end
A(i,i) = -g(2,2,psi); %bottom right element (always =-1 in Lambda case)

end

%Storing tree heights
hts(trials)=Pi*(-inv(A))*ones(1,length(Pi))';

%Storing total branch lengths
lengths(trials) = Pi*(-inv(A))*R;

%Calculating second moment for the tree height
hts_second_moment(trials) = 2*Pi*(-inv(A))^2*ones(1,length(Pi))';

end

%Final tree heights, total branch lengths, second moments for pop. sizes 10:10:100
final_heights = nonzeros(hts);
final_branchlengths = nonzeros(lengths);
final_second_moments = nonzeros(hts_second_moment);

%This function determines the value of each element (for our matrix A)
function g_ki = g(k,i,psi)

 %Rates for psi-coalescent
 g_ki = nchoosek(k,i)*psi^(i-2)*(1-psi)^(k-i); %Psi coalescent
 
 %Rates for beta-coalescent
 %g_ki = (nchoosek(k,i))*(beta(i-alpha,k-i+alpha))/(beta(alpha,2-alpha));
end

