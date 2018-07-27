
%This is a modification of Transfunc.m in order to test results
%This is #4 on "Project Ideas"

Hd = input('Enter Matrix Hd as nxn >')
Hs = input('Enter Matrix Hs as nxn >')
E_real = input('Enter range as a 1x2 matrix of E real>')


a = E_real(1,1);
b = E_real(1,2);
z=1;

Hsdagger = (inv(Hs))' ;

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
for j = a:b
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
 lambda(:,1) = log(lambda(:,1))/i;
 scatter3( real(lambda(:,1)), imag(lambda(:,1)), lambda(:,2),'filled') 
