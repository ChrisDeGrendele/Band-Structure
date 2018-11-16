%alpha = 0.8;
%beta_1 = -.4;
%beta_2 = -1.3;

%Test system 1
%Hd = [alpha, beta_1;beta_1,-alpha];
%Hs = [0,beta_2;beta_2,0];
%min_epsilon = 1e-3;
%ideal_spacing = .02;
%E = [-4,4];

%[Data_x, Data_y] = amr_KvsE2(Hd,Hs,E,min_epsilon, ideal_spacing);



function [Data_x, Data_y] = amr_KvsE(Hd,Hs,E,min_epsilon, ideal_spacing)

a = E(1,1);
b = E(1,2);
rank_HS = rank(Hs);
%Storage for solutions
Data_x = [];
Data_y = [];


epsilon = .1; %Starting Epsilon
e = a-epsilon; %Starting energy

while e < b
    e = e+epsilon;   %Start at a. Because we start at e = a-epsilon.
    Hsdagger = (Hs)' ;
    length = size(Hsdagger,1);
    for j = 1:length  %Calculate Hermitian Conjugate
        for k = 1:length
            Hsdagger(j,k) = conj(Hsdagger(j,k));
        end
    end
    
    %BUILD AV = BVD (GENERALIZED EIGEN VALUE PROBLEM)
    A = [e*eye(length) - Hd, -Hs; eye(length), zeros(length)];
    B = [Hsdagger, zeros(length); zeros(length), eye(length)];
    [~,D] = eig(A,B); %AV = BVD
    
    D_array = [];
    for j = 1:size(D,2)
        D_array = [D_array, D(j,j) ];
    end
    
    D_array = sort(D_array);
    
    for j = 2:size(D_array,2)
        diff = abs( abs(log( D_array(j) ) /(1i)) - abs( log( D_array(j) )/(1i) ) );
        if epsilon <= min_epsilon
            bad_data = false;
            break
        elseif diff <= ideal_spacing
            bad_data = true;  %Bad => the gap may be too big. 
            epsilon = epsilon/2; %Halfs step
            break
        elseif 2*diff < ideal_spacing
            bad_data = true; %Bad => TOO refined
            epsilon = epsilon*2; %Halfs step
    
            break
        else
            bad_data = false;
        end
    end
    
    if bad_data == false 
        
        temp = [];
        for j = 1:size(D,2)
            temp = [temp, D(j,j)];
        end
        temp = sort(temp);
        for j = 1:rank_HS
            index1 = size(D,2)/2 - j +1;
            index2 = size(D,2)/2 +j;

            Data_x = [Data_x, log( temp(index1) )/(1i), log( temp(index2) )/(1i)];
            Data_y = [Data_y, e, e];    
        end
        
        e = e + epsilon; %Only step forward if data is good
    end

            
end
end