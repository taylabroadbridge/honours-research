A = [];

function[] = Kingman_theoretical(nn) 

A(nn,nn) = -1;

for i = 1:nn
    
    for j = 1:nn
        
        if i==j && i~=nn
            A(i,j) = -nchoosek(n,2);
        end
        
        if j == i+1 
            A(i,j) = nchoosek(n,2);
        end
        
        
    end
end
            
end           
            
            