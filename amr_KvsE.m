%alpha = 0;
%beta = -.5;
%Hd = [alpha, beta, 0, 0, 0, beta ;beta,alpha,beta, 0, 0, 0; 0, ...
%    beta, alpha, beta, 0, 0; 0, 0, beta, alpha, beta, 0; 0, 0, 0, ...
%    beta, alpha, beta; beta, 0, 0, 0, beta, alpha];
%Hs = zeros(6);
%Hs(1,3) = beta;


%min_epsilon = 1e-3;
%ideal_spacing = .02;

%Data_x, Data_y] = amr_KvsE2(Hd,Hs,E,min_epsilon, ideal_spacing);

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
    num = size(Hsdagger,1);
    
    %BUILD AV = BVD (GENERALIZED EIGEN VALUE PROBLEM)
    A = [e*eye(num) - Hd, -Hs; eye(num), zeros(num)];
    B = [Hsdagger, zeros(num); zeros(num), eye(num)];
    [V,D] = eig(A,B); %AV = BVD
    

    %%%%%%%%%
    D_array = [];
    Simple_D = [];
    
        for j = 1:size(D,2)
            D_array = [D_array, D(j,j)];
        end
        
        D_array = sort(D_array, 'ComparisonMethod', 'abs');
        
        for j = 1:rank_HS
            index1 = size(D_array,2)/2 - j +1;
            index2 = size(D_array,2)/2 +j;
            
            Simple_D = [Simple_D, log( D_array(index1) )/(1i), log( D_array(index2) )/(1i)];    
        end
    
    %%%%%%%%
    
   
    
    for j = 2:size(Simple_D,2)
        diff = abs( log( Simple_D(j) ) /(1i) -  log( Simple_D(j-1) )/(1i)  );
        if epsilon <= min_epsilon
            bad_data = false;
        elseif diff >= ideal_spacing
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
        
        for j = 1:size(Simple_D,2)
            Data_x = [Data_x, Simple_D(j) ];
            Data_y = [Data_y, e];
        end
        %temp = sort(temp);
        %for j = 1:rank_HS
        %    index1 = size(D,2)/2 - j +1;
        %    index2 = size(D,2)/2 +j;

        %    Data_x = [Data_x, log( temp(index1) )/(1i), log( temp(index2) )/(1i)];
        %    Data_y = [Data_y, e, e];    
        %end
        
        e = e + epsilon; %Only step forward if data is good
    end

            
end
end