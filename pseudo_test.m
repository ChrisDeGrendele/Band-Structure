

alpha = .8;
beta_1 = -.4;
beta_2 = -1.3;
E = [-4,-3];
epsilon = 5e-1;
min_epsilon = 5e-1;
ideal_spacing = 4e-1;

Data_x = [];
Data_y = [];

Hd = [alpha, beta_1;beta_1,-alpha];
Hs = [0,beta_2;beta_2,0];

%Hd = [alpha, 0, beta_2, 0;
%      0, alpha, beta_1, beta_2;
%      beta_2, beta_1, -alpha, 0;
%      0, beta_2, 0, -alpha];
%Hs = [0, 0, beta_2, beta_1;
%      0, 0, 0, beta_2;
%      0, 0, 0, 0;
%      0, 0, 0, 0;];
  


rank_HS = rank(Hs);

T_E = Build_Te(Hd, Hs, E);
[V,D] = eig(T_E);
fprintf("V:\n")
abs(V)

fprintf("D:\n")
abs(D)

D_array = [];
Simple_D = []; %Not sure what this is for...

%append diagonol values to an array
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


Imag_k = [];
Imag_E = [];
Real_k = [];
Real_E = [];

for e = 1:size(Data_x,2)
    if abs(imag( Data_x(e) )) > 1.e-6
        Imag_k = [Imag_k, imag(Data_x(e)) ];
        Imag_E = [Imag_E, Data_y(e)];
    else
        Real_k = [Real_k, Data_x(e)];
        Real_E = [Real_E, Data_y(e)];
    end
end



%Sorts data
[Real_k, sortIndex] = sort(Real_k);
Real_E = Real_E(sortIndex);


[Imag_k, sortIndex] = sort(Imag_k);
Imag_E = Imag_E(sortIndex);

%Screens out zeros
i = 1;
eps_forzero = .01;
while i < size(Real_k,2)
    
    if i > size(Real_k,2) %End of array reached
        break
    elseif abs( Real_k(i) ) - eps_forzero < 0
        Real_E(i) = [];
        Real_k(i) = [];
    elseif abs( Real_k(i) ) - eps_forzero > 0
        i = i+1; %Only step forward if condition holds
    end
    
end

% A copy of above, just for Imaginary values
%maybe make this into a function for cleaner code, but it's only used 
%twice so idk
i = 1;
eps_forzero = .01;
while i < size(Imag_k,2)
    
    if i > size(Imag_k,2) %End of array reached
        break
    elseif abs( Imag_k(i) ) - eps_forzero < 0
        Imag_E(i) = [];
        Imag_k(i) = [];
    elseif abs( Imag_k(i) ) - eps_forzero > 0
        i = i+1; %Only step forward if condition holds
    end
    
end
    



figure(1);
ax1=subplot('Position',[.4 .3 .3 .3]); %Position = [Left bottom width height];
scatter(abs(Real_k), Real_E, '.')
title('Real k v. E')
ylim([a,b])

hold on
ax2=subplot('Position',[.1 .3 .3 .3]);
set(ax2,'YTick',[],'XTick',[]);
scatter( -abs(Imag_k), Imag_E, '.')
title('Imaginary k v. E')
ylim([a,b])

%hold off
%plot(abs(Real_k), Real_E)


save_val = true;
if save_val == true
    save('/Users/chris/Desktop/assym.mat','Imag_E','Imag_k', 'Real_E', 'Real_k');
    %saveas(1,'/Users/chris/Documents/GitHub/Band-Structure/pictures/meta.jpg')
end

return
%Below is stuff for calculating branch point 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

Imag_k = abs(Imag_k);
Real_k = abs(Real_k);

Pure_Imag_k = Imag_k;
Pure_Imag_E = Imag_E;


%MAX LOOP VALUE:
input1 = input('Input rough bounds for E of max loop value [lower,higher] = ');
low = input1(1,1);
high = input1(1,2);

%Screen out all the zeros (again)
stop = 0;
while stop ==0
stop =1;
for j = 1:size(Imag_E,2)
    if size(Imag_k,2) < j
        stop=0;
        break
    elseif Imag_k(j) == 0 
        Imag_k(j) = [];
        Imag_E(j) = [];
        stop = 0;
    elseif (Imag_E(j) > high) | (Imag_E(j) < low)
        Imag_k(j) = [];
        Imag_E(j) = [];
        stop = 0;
    end
end
end

[min_k,loc] = min(Imag_k);
E_upper = Imag_E(loc);


%LOWER LOOP VALUE
input2 = input('Input rough bounds for E of min loop value [lower,higher] = ');
low = input2(1,1);
high = input2(1,2);
Imag_k = Pure_Imag_k;
Imag_E = Pure_Imag_E;

%Screen out all the zeros (again)
stop = 0;
while stop ==0
stop =1;
for j = 1:size(Imag_E,2)
    if size(Imag_k,2) < j
        stop=0;
        break
    elseif Imag_k(j) == 0 
        Imag_k(j) = [];
        Imag_E(j) = [];
        stop = 0;
    elseif (Imag_E(j) > high) | (Imag_E(j) < low)
        Imag_k(j) = [];
        Imag_E(j) = [];
        stop = 0;
    end
end
end
[min_k,loc] = min(Imag_k);
E_lower = Imag_E(loc);


%BRANCH POINT
input3 = input('Input rough bounds for the branch point [lower,higher] = ');
low = input3(1,1);
high = input3(1,2);
Imag_k = Pure_Imag_k;
Imag_E = Pure_Imag_E;

%Screen out all the zeros (again)
i = 1;
eps_forzero = .01;
while i < size(Imag_k,2)
    
    if i > size(Imag_k,2) %End of array reached
        break
    elseif abs( Imag_k(i) ) - eps_forzero < 0
        Imag_E(i) = [];
        Imag_k(i) = [];
    elseif abs( Imag_k(i) ) - eps_forzero > 0
        i = i+1; %Only step forward if condition holds
    end
    
end
    

%Calculate branch point
[K_branch,loc] = max(Imag_k);
E_branch = Imag_E(loc);

fprintf("The branch point is: (%f, %f)\n", K_branch, E_branch )
fprintf("E_upper = %f, E_lower = %f\n", E_upper, E_lower)

%Returns a value [0,1] of position of branch where 0.5 is center (E value)
Branch_pos = (E_branch - E_lower)/(E_upper - E_lower);

%Returns a ratio of width/height. Bigger = longer, smaller = "stubbier"
wid_hei_ratio = K_branch/(E_upper - E_lower);

fprintf("Normalied Branch position = %f\n", Branch_pos)
fprintf("Ratio of width/height = %f\n", wid_hei_ratio)






  
