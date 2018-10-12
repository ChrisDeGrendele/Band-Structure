clear all
alpha =-3
beta_1 = -2.2
beta_2 = 1.5
%Test system 1
Hd = [alpha, beta_1;beta_1,-alpha];
Hs = [0,beta_2;beta_2,0];

%Test System 2
%Hd = [alpha, beta_1;beta_1,alpha];
%Hs = [beta_2,0;0,1];

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

%SCREEN OUT ZEROS
Real_k = real(Data_x);
Real_E = Data_y; 
sizee = size(Real_k,2);
stop=0
while stop == 0 
stop = 1;
for j = 1:sizee
    if size(Real_k,2) < j
        break
    elseif abs( Real_k(j) ) - .01 < 0 
        Real_k(j) = [];
        Real_E(j) = [];
        stop=0;
    end
end

end


Imag_k = imag(Data_x);
Imag_E = Data_y;
sizee = size(Imag_k,2);
stop=0;
n=1


while stop ==0 
stop = 1;
for j = 1:sizee
    if size(Imag_k,2) < j
        break
    elseif abs( Imag_k(j) ) - .01 < 0 
        Imag_k(j) = [];
        Imag_E(j) = [];
        stop=0;
    end
end

end



figure;
subplot(2,2,2);
scatter(abs(Real_k), Real_E)
title('Real k v. E')

subplot(2,2,1);
scatter( -abs(Imag_k), Imag_E)
title('Imaginary k v. E')

Imag_k = abs(Imag_k);
Real_k = abs(Real_k)'

input = input('Input rough bounds for the branch point [lower,higher] = ');
low = input(1,1);
high = input(1,2);

%Sort into order
[sortedX, sortIndex] = sort(Imag_k);
sortedY = Imag_E(sortIndex);


for j = 1:size(Imag_E,2)
    if Imag_E(j) > high
        Imag_k(j) = 0;
        stop=0;
    elseif Imag_E(j) < low
        Imag_k(j) = 0;
        stop=0;
    end
end

while stop ==0
stop =1;
for j = 1:size(Imag_E,2)
    if size(Imag_k,2) < j
        stop=0;
        break
    elseif Imag_k(j) ==0 
        Imag_k(j) = [];
        Imag_E(j) = [];
        stop = 0;
    end
end
end


Number = size(Imag_k,1); 
Slope=0;
for j = 1:(Number-1)
    if j ==1
        Slope = (Imag_E(j+1) -Imag_E(j)) / (Imag_k(j+1) - Imag_k(j));
    else
        Slope = [Slope, (Imag_E(j+1) -Imag_E(j)) / (Imag_k(j+1) - Imag_k(j))];
    end
end

[Max_Slope,Ind] = max(Slope); 

fprintf("The branch point is: (%f, %f)", Imag_k(Ind+1), Imag_E(Ind+1) )

[XX,YY] = max(Imag_k);
YYY = Imag_E(YY);

XX
YYY



