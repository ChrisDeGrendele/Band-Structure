%ABA Transient effects for Matt

clear all
close all
%Values

alpha = 2;  
beta_1 = 2;
beta_2 = 3;
E = 1 + .000001*i;
%Matrices
Hd = alpha;
Hs_a = beta_1;  
Hs_b = beta_2;

%Containers for storage of solutions
Data_ABA = [];
Datak = [];


[Va, Da] = buildVDV(Hd, Hs_a, E);
[Vb, Db] = buildVDV(Hd, Hs_b, E);


for k = 0:100
    T_ABA = (Va*Da*Va^-1) * ((Vb*Db*Vb^-1) * (Va*Da*Va^-1))^k;
    [V_orig,D_orig] = eig(T_ABA);
    [V,D] = screenout(V_orig,D_orig,Hs_a);
    T_ABA = V*D*V^-1;
    Datak = [Datak,k];
    Data_ABA = [Data_ABA, norm(T_ABA)]; 
end



plot(Datak, Data_ABA)


function [V,D] = screenout(V,D,Hs)

    rank_HS = rank(Hs);
    D_array = diag(D)'; 
    %Sort D_array and I which is a permuation array
    [D_array, I] = sort(D_array, 'ComparisonMethod', 'abs');
    %sorts V columns corresponding to D
    temp = V;
    for i = 1:length(I)
        V(:,i) = temp(:,I(i));
    end
    %Index of half sorted array
    index2 = length(D_array)/2;
    %Sets all the data points above half and below half-Rank_HS = 0
    for i = 1:length(D_array)
        if  i > index2
            D_array(i) = 0;
        elseif i < index2-rank_HS
            D_array(i) = 0;
        end
    end
    

    %Diagonlize array
    D = diag(D_array);
end

function [V1,D1] = buildVDV(Hd, Hs, e)
    num = size(Hs,1);
    A = [e*eye(num) - Hd, -Hs; eye(num), zeros(num)];
    B = [Hs', zeros(num); zeros(num), eye(num)];
    [V,D] = eig(A,B); %AV = BVD
   
    V1=V; %Set outputs values
    D1=D;
end
    
    
    






