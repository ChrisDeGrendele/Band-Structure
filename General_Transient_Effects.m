%A SIMPLE TWO BAND STRUCTURE
alpha = 0;  %Try 0, try -1/2.
beta_1 = -.3;
beta_2 = -1;
Hd = [alpha, beta_1;beta_1,alpha];
Hs = [0,beta_2;beta_2,0];
%E = 0.7 + .000001*i;
E = 0.7;

[V,D] = buildVDV(Hd,Hs,E);
%Te = Build_Te(Hd,Hs,E);
%[V,D] = eig(Te);



Data_k = [];
Data_norm = [];

for k = 1:20
    Data_norm = [Data_norm,  norm( V* (D^k) *V^-1 )];
    Data_k = [Data_k, k];
end
plot(Data_k, Data_norm)
    

%Below is the Ta Tb stuff
%We're going to have 1 Hd and 2 Hs For an A(BA)^k molecule

%Hd = 
%Hs_a = 
%Hs_b = 

%[Va, Da] = buildVDV(Hd, Hs_a, E);
%[Vb, Db] = buildVDV(Hd, Hs_b, E);

%Data_ABA = [];
%Datak = [];
%T_Ba = (Vb*Db*Vb^-1) * (Va*Da*Va^-1);
%T_A = (Va*Da*Va^-1);
%T_Batok = 1; %Starts at 1

%for k = 0:8
%    Teff = T_A * T_batok;  
%    Datak = [Datak,k];
%    Data_ABA = [Data_ABA, norm(Teff)];
%    T_batok = T_batok*T_ba; %Each k multiply in a new ab. This saves computing power 
%end

%plot(Datak, Data_ABA)






function [V1,D1] = buildVDV(Hd, Hs, e)
    num = size(Hs,1);
    
    A = [e*eye(num) - Hd, -Hs; eye(num), zeros(num)];
    B = [Hs', zeros(num); zeros(num), eye(num)];
    [V,D] = eig(A,B); %AV = BVD
    
    rank_HS = rank(Hs);
    D_array = diag(D)';
    [D_array, I] = sort(D_array);
    
    %sorts V correspoding to D
    temp = V;
    for i = 1:length(I)
        V(:,i) = temp(:,I(i));
    end
    
    
    index1 = size(D_array,2)/2 - rank_HS +1;
    index2 = size(D_array,2)/2 + rank_HS;
    
    %Sets all the data points not in 2r to 0
    for i = 1:length(D_array)
        if i < index1 || i > index2
            D_array(i) = 0;
        end
    end
    
    D = diag(D_array);
    %Screen out all eigen values >1
    for i = 1:length(D)
        if (real(D(i,i))^2 + imag(D(i,i))^2)^(1/2) > 1.0001
            D(i,i) = 0;
        end
    end
    
    V1 = V;
    D1 = D;
end
    
    
    






