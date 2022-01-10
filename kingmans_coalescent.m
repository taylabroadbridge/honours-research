%Algorithm for Kingman's coalescent

%Initialise tree heights and branch lengths
meantreeheight = [];
meanbranchlength = [];

%Inititalising generations
generations = [];
new_generation = [];
treeheight = [];
branchlength = [];

%Iterate through population sizes
for n = 10:2:100

 %Average 1000 samples for each population size
for r = 1:1000

%Populate the original generation (generation k)
for i = 1:n
    for j = 1:n
    generations(j,i) = i;
    end
end


time_counter = 0;
T = []; %time counter
S = []; %store the times of the jumps

k = n;
index = 1;

while k > 1

    
    %generate waiting time until coalescence
    T_k = exprnd(1/nchoosek(k,2));  

    %T vector holds times BETWEEN coalescence
    T = [T T_k];
    
    time_counter = time_counter + T_k;
    
    %S vector holds times OF coalescent events 
    S = [S time_counter];
    
    l = datasample(generations(index,:),1);
    temp = setdiff(generations(index,:),l);   
    m = datasample(temp,1); %since 1 <= l < m <= k. This line chooses random number st. its Neq to l
    
    
    locate = find(generations(index,:) == min(l,m));
    

    index = index + 1;
    
    for i = index:n
    generations(i,locate) = max(l,m);
    end

    k = k-1;
    
   

end

%Calculate tree height
height = sum(T(1:n-1));

L=0;
%for i = n:-1:2
%    L = L + i*T(i);
%end

%Calculate total branch length
L = [n:-1:2]*T';


treeheight(r) = height; %for n=10, this is approx 2*(1-1/n)
branchlength(r) = L; %for n=10, 

end

%Averaging the tree heights and branch lengths
meantreeheight = [meantreeheight mean(treeheight)];
meanbranchlength = [meanbranchlength mean(branchlength)];

end

exptreeheight = [];
expbranchlength = [];

%Comparing to the theoretical
for n=10:2:100
    counter = 0;
    exptreeheight = [exptreeheight 2*(1-1/n)]

    for j=1:n-1
        counter = counter + 1/j ;
    end
    expbranchlength = [expbranchlength 2*counter]
end

