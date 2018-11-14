function [Data_x, Data_y] = amr_KvsE(Hd,Hs,E,min_epsilon, min_spacing)

a = E(1,1);
b = E(1,2);
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
    
    
    for j = 2:size(D,2)
        diff = abs( abs( log( D(j,j) )/(1i)) - abs( log( D(j-1,j-1) )/(1i)));
        if min_epsilon <= min_spacing
            bad_data = true;
            break
        elseif diff > min_spacing
            bad_data = true;  %Bad => the gap may be too big. 
            epsilon = epsilon/2; %Halfs step
            e = e-epsilon; %reverts back.
            break
        elseif 2*diff < min_spacing
            bad_data = true; %Bad => TOO refined
            epsilon = epsilon*2; %Halfs step
            e = e-epsilon; %reverts back.
            break
        else
            bad_data = false;
        end
    end
    
    if bad_data == false
        for j = 1:size(D,2)
            Data_x = [Data_x, (log( D(j,j) )/(1i)) ]; %STORE DATA
            Data_y = [Data_y, e];
        end
        e = e + epsilon;
    end

            
end