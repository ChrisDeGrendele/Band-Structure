clear all
%alpha = -3;
%beta_1 = -2.2;  %THIS GIVES WEIRD STUFF?
%beta_2 = 1.5;

alpha = 0;
beta_1 = -1;
beta_2 = -3;

%A SIMPLE TWO BAND STRUCTURE
Hd = [alpha, beta_1;beta_1,-alpha];
Hs = [0,beta_2;0,0];

%Test System 2
%Hd = [alpha, beta_1;beta_1,alpha];
%Hs = [beta_2,0;0,0];

%Test System 3 
%alpha = 0.69;
%beta_1 = -1.2;
%beta_2 = -.7;
%Hd = [alpha, beta_1;beta_1,alpha];
%Hs = [beta_2,0;0,0];
%Hd = [alpha, beta_1;beta_1,-alpha];
%Hs = [0,beta_2;beta_2,0];

%Test System 4
%alpha = 0.69;
%beta_1 = -1.35;
%beta_2 = .92;
%Hd = [alpha, beta_1;beta_1,alpha];
%Hs = [beta_2,0;0,0];
%Hd = [alpha, beta_1;beta_1,-alpha];
%Hs = [0,beta_2;beta_2,0];

%Test System 5
%alpha = 0.21;
%beta_1 = -2.43;
%beta_2 = .13;
%Hd = [alpha, beta_1;beta_1,alpha];
%Hs = [beta_2,0;0,0];
%Hd = [alpha, beta_1;beta_1,-alpha];
%Hs = [0,beta_2;beta_2,0];

%Model System 3 PARA
%alpha = 0;
%beta = -.5;
%Hd = [alpha, beta, 0, 0, 0, beta ;beta,alpha,beta, 0, 0, 0; 0, beta, ...
%    alpha, beta, 0, 0; 0, 0, beta, alpha, beta, 0; 0, 0, 0, beta, alpha,...
%    beta; beta, 0, 0, 0, beta, alpha];
%Hs = zeros(6);
%Hs(1,4) = beta;

%Model System 4 META
alpha = 0;
beta = -.5;
Hd = [alpha, beta, 0, 0, 0, beta ;beta,alpha,beta, 0, 0, 0; 0, ...
    beta, alpha, beta, 0, 0; 0, 0, beta, alpha, beta, 0; 0, 0, 0, ...
    beta, alpha, beta; beta, 0, 0, 0, beta, alpha];
Hs = zeros(6);
Hs(1,3) = beta;

%Model System 5 ORTHO
%alpha = 0;
%beta = -.5;
%Hd = [alpha, beta, 0, 0, 0, beta ;beta,alpha,beta, 0, 0, 0; 0, ...
%    beta, alpha, beta, 0, 0; 0, 0, beta, alpha, beta, 0; 0, 0, 0, ...
%    beta, alpha, beta; beta, 0, 0, 0, beta, alpha];
%Hs = zeros(6);
%Hs(1,2) = beta;

%Model System 6
%Hd = -.5* [0,1,0,0,0,1;1,0,1,0,0,0;0,1,0,1,0,0;0,0,1,0,1,0;0,0,0,1,0,1;1,0,0,0,1,0];
%Hs = zeros(6);
%Hs(2,1) = -.5;
%Hs(3,4) = -.5;
%Hs(5,4)= -.5;
%Hs(6,1) = -.5;

%Data
%Hd = [.8,-.4;-.4,-.8];
%Hs = [0,-1.3;-1.3,0];
%E = [-10,10];

%Hd = [0];
%Hs = [-1/2];
E=[-5,5];
a = E(1,1);
b = E(1,2);

  
    
%DELETED E LOOP HERE. 
%CALL FUNCTION THAT WILL BUILD DATA_X AND DATA_Y
%WHERE DATA X = COMPLEX_K AND DATA Y = COMPLEX_Y

%Build imaginary and real data to plot. 
%This will plot a point with imag to imag
min_epsilon = 1e-3;
ideal_spacing = .02;


[Data_x, Data_y] = amr_KvsE(Hd,Hs,E,min_epsilon, ideal_spacing);
%And only purely real to real

Imag_k=  [];
Imag_E = [];
Real_k = [];
Real_E = [];
for e = 1:size(Data_x,2)
    if imag( Data_x(e) ) ~= 0
        Imag_k = [Imag_k, imag(Data_x(e)) ];
        Imag_E = [Imag_E, Data_y(e)];
    elseif imag( Data_x(e) ) == 0
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


save_val = true;
if save_val == true
    save('/Users/chris/Documents/GitHub/Band-Structure/pictures/meta.mat','Imag_k', ...
        'Imag_E','Real_k','Real_E');
    saveas(1,'/Users/chris/Documents/GitHub/Band-Structure/pictures/meta.jpg')
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






