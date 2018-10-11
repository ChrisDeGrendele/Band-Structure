clear all
alpha =-.3
beta_1 = -2.2
beta_2 = 1.5
%Test system 1
%Hd = [alpha, beta_1;beta_1,-alpha];
%Hs = [0,beta_2;beta_2,0];

%Test System 2
Hd = [alpha, beta_1;beta_1,alpha];
Hs = [beta_2,0;0,1];

%Data
%Hd = [.8,-.4;-.4,-.8];
%Hs = [0,-1.3;-1.3,0];
%E = [-10,10];

%Hd = [0];
%Hs = [-1/2];
E=[-5,5];

    a = E(1,1);
    b = E(1,2);
  
    
for e = a:.1:b 
    Te = Build_Te(Hd,Hs,e);
    [V,D] = eig(Te);

        
        %#Number 1 2d plot of E vs k
    [count, count2] = size(D);
    
    if (a==e)
        for j = 1:count
           Data_x(1,j) = log( D(j,j) )/(1i);
           Data_y(1,j) = e;
           j = j+1;
        end 
    else 
        for j = 1:count
            Data_x = [Data_x, (log( D(j,j) )/(1i)) ];
            Data_y = [Data_y, e];
        end
    end
            
end

figure;
subplot(2,2,2);
scatter(real(Data_x), Data_y)
title('Real k v. E')

subplot(2,2,1);
scatter(imag(Data_x), Data_y)
title('Imaginary k v. E')


%scatter(real(Data_x), Data_y)
%figure
%scatter(imag(Data_x), Data_y)

