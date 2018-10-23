clear all
alpha = -3;
beta_1 = -2.2;  %THIS GIVES WEIRD STUFF?
beta_2 = 1.5;

%alpha = 0.8;
%beta_1 = -.4;
%beta_2 = -1.3;

%Test system 1
Hd = [alpha, beta_1;beta_1,-alpha];
Hs = [0,beta_2;beta_2,0];

%Test System 2
%Hd = [alpha, beta_1;beta_1,alpha];
%Hs = [beta_2,0;0,0];

%Data
%Hd = [.8,-.4;-.4,-.8];
%Hs = [0,-1.3;-1.3,0];
%E = [-10,10];

%Hd = [0];
%Hs = [-1/2];
E=[-10,10];

    a = E(1,1);
    b = E(1,2);
  
    
for e = a:.01:b 
    
    Hsdagger = (Hs)' ;
    length = size(Hsdagger,1);

    %Calculate Inverse Hermitian Conjugate
    for j = 1:length
        for k = 1:length
            Hsdagger(j,k) = conj(Hsdagger(j,k));
        end
    end
    
    %BUILD AV = BVD (GENERALIZED EIGEN VALUE PROBLEM)
    A = [e*eye(length) - Hd, -Hs; eye(length), zeros(length)];
    B = [Hsdagger, zeros(length); zeros(length), eye(length)];
    [V,D] = eig(A,B); %AV = BVD

        
    %Take the 0 and infinity eigen values out.
  %  for i = 1:size(D)
  %      diag =sort( diag(D)' );
        
   %     if abs( diag(1,i)  ) < .01
            
        
    %end
  
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

%Build imaginary and real data to plot. 
%This will plot a point with imag to imag
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


%SCREEN OUT 0 , pi, pi/2 for real data

stop=0;
[Real_k, sortIndex] = sort(Real_k);
Real_E = Real_E(sortIndex);
while stop == 0 
for j = 1:size(Real_k,2)
    if size(Real_k,2) < j
        break
    elseif abs( Real_k(j) ) - .01 < 0   %OUT WITH 0 
        Real_k(j) = [];
        Real_E(j) = [];
        Real_k(end) = [];   %Screens out zeros and infty eigen values.
        Real_E(end) = [];
        break
    elseif abs (Real_k(j) - pi*.5)  < .01   %OUT WITH PI/2
        Real_k(j) = [];
        Real_E(j) = [];
        stop=0;
    elseif abs( Real_k(j)  - pi)  < 0.01     %OUT WITH PI
        Real_k(j) = [];
        Real_E(j) = [];
        stop=0;
    else
        stop = 1;
    end
end

end


%Sorts data
[Imag_k, sortIndex] = sort(Imag_k);
Imag_E = Imag_E(sortIndex);
%Same as above, copy-pasted, screens out zeros for imag data now
stop=0;
while stop ==0 
for j = 1:size(Imag_k,2)
    flag = 0;
    if size(Imag_k,2) < j
        flag =1; 
        break
    elseif abs( Imag_k(j) ) < .1
        Imag_k(j) = [];
        Imag_E(j) = [];
        Imag_k(end) = []; %Screens out zeros and corresponding "inifity"
        Imag_E(end) = []; 
        flag = 1;
        break
    else
        stop=1;
    end
    if flag == 1; %Work around to break out of outer loop and restart process
        break
    end
end
end



figure;
a=subplot('Position',[.4 .3 .3 .3]); %Position = [Left bottom width height];
scatter(abs(Real_k), Real_E, '.')
title('Real k v. E')
hold on
b=subplot('Position',[.1 .3 .3 .3]);
set(b,'YTick',[],'XTick',[]);
scatter( -abs(Imag_k), Imag_E, '.')
title('Imaginary k v. E')


Imag_k = abs(Imag_k);
Real_k = abs(Real_k);

input = input('Input rough bounds for the branch point [lower,higher] = ');
low = input(1,1);
high = input(1,2);



%Screen out big / small values for calculating branch point, set them = 0
for j = 1:size(Imag_E,2)
    if Imag_E(j) > high
        Imag_k(j) = 0;
        stop=0;
    elseif Imag_E(j) < low
        Imag_k(j) = 0;
        stop=0;
    end
end

%Screen out all the zeros (again)
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
    end
end
end

%Calculate branch point
[Branch_k,loc] = max(Imag_k);
Branch_E = Imag_E(loc);

fprintf("The branch point is: (%f, %f)\n", Branch_k, Branch_E )




