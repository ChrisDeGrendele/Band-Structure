close all;
clear all;
%This is a modification of Transfunc.m in order to test results
%This is #4 on "Project Ideas"
%This is now modified to try pseudo inverse.


%No difference between this example inv() and pinv()
%alpha = .8;
%beta_1 = -.4;
%beta_2 = -1.3;
%Hd = [alpha, beta_1;beta_1,-alpha];
%Hs = [0,beta_2;beta_2,0];


%Josh's Example: Vertical lines? Doesn't work...
alpha = .8;
beta_1 = -.4;
beta_2 = -1.3;
Hd = [alpha, 0, beta_2, 0;
      0, alpha, beta_1, beta_2;
      beta_2, beta_1, -alpha, 0;
      0, beta_2, 0, -alpha];
Hs = [0, 0, beta_2, beta_1;
      0, 0, 0, beta_2;
      0, 0, 0, 0;
      0, 0, 0, 0;];



E_real = [-5,5];
dx = .001;



a = E_real(1,1);
b = E_real(1,2);
z=1;

%PSEUDO VS REAL INV???
%Hsdagger = (inv(Hs))' ;
Hsdagger = (pinv(Hs))';

[m,n] = size(Hsdagger);
I = eye(m,n);

%Calculate Inverse Hermitian Conjugate
for j = 1:m
    for k = 1:n 
        Hsdagger(j,k) = conj(Hsdagger(j,k));
    end
end


%resolution = 400;%Resolution
%output_size = [3000 1800];%Size in pixels

%Build T(E) real
for j = a:dx:b
    First = Hsdagger*(j*I-Hd);
    Second = -Hsdagger*Hs;
    Third = I;
    Fourth = zeros(m,n);
    T = [First, Second; Third, Fourth];
    
    e = eig(T); %returns column vector of eigen values
    e2 = ones(size(e))*j;
    e = [e(:),e2(:)];
    
    %Creates lambda as a matrix
       %First column is eigen values
       %Second Column is the corresponding E
    if j == a
        lambda = e;
    else    
        lambda = [lambda;e];
    end
   
end 
%Calculates k
lambda(:,1) = log(lambda(:,1))/i;

real_k = [];
real_e = [];
imag_k = [];
imag_e = [];

for j = 1:size(lambda,1)
    if isreal( lambda(j,1) )
        real_k = [real_k, abs(real( lambda(j,1) ))];
        real_e = [real_e, lambda(j,2)];
    else
        imag_k = [imag_k, abs(imag( lambda(j,1) ))];
        imag_e = [imag_e, lambda(j,2)];
    end
end

scatter(real_k, real_e)
figure(2)
scatter(imag_k, imag_e)
    
    %scatter3( real(lambda(:,1)), imag(lambda(:,1)), lambda(:,2),'filled')
%scatter(real(lambda(:,1)), lambda(:,2))
