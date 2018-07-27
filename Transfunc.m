%function [out1,out2] = function_name(in1,in2)
%function [Te] = Build_Te()

%Test


Hd = input('Enter Matrix Hd as nxn >')
Hs = input('Enter Matrix Hs as nxn >')
E_real = input('Enter range as a 1x2 matrix of E real>')
%E_imag = input('Enter range as a 1x2 matrix of E imaginary>')


a = E_real(1,1);
b = E_real(1,2);
%c = E_imag(1,1);
%d = E_imag(1,2);
z=a;
%zz=c;

Hsdagger1 = (inv(Hs)) ;
Hsdagger = (Hsdagger1)' ;

[m,n] = size(Hsdagger);
I = eye(m,n);

%Calculate Inverse Hermitian Conjugate
for j = 1:m
    for k = 1:n 
        Hsdagger(j,k) = conj(Hsdagger(j,k));
    end
end


resolution = 400;%Resolution
output_size = [3000 1800];%Size in pixels

%Build T(E) real
for e = a:b
    First = Hsdagger*((e+.0001i)*I-Hd);
    Second = -Hsdagger*Hs;
    Third = I;
    Fourth = zeros(m,n);
    T = [First, Second; Third, Fourth];
    
 
    
    opts.npts=100;
    %opts.levels=-7.5:0.5:-2;
    %opts.ax=[-1 0.2 -1.05 0];
    %opts.proj_lev=Inf;
    %opts.colour=1;
    %opts.thick_lines=1;
    %opts.scale_equal=1;
    %opts.grid=0;
    %opts.dim=1;
    %opts.no_ews=0;
    %opts.no_psa=0;
    opts.fov=1;
    %opts.direct=1;
    %opts.print_plot=1;
    %eigtool(T,opts)
    %set(gcf,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
    %print(compose("Treal_%d.png",z),'-dpng',['-r' num2str(resolution)]);
    %saveas(fig, compose("Treal_%",z), 'png')
    z=z+1;
    
    
    
    
    %Copied from below
    
    
   %Finding P(a) which is max of ||A^k|| vs k
    [V,D] = eig(T);
    %invv = inv(V);
    k = 1;
    
    %Find D' (eigen values > 1 = 0)
    [XX,YY] = size(D);
    i =1;
    while i <= XX
    j=1;
       while j <=YY
          Mag = ( (real( D(i,j) ))^2 + (imag( D(i,j) ))^2 )^(1/2);
          if Mag < 1.0
              Dprime(i,j) = D(i,j);
          else                  
              Dprime(i,j) = 0;
          end
          j = j+1;
        end
        i = i+1;
    end
    
        k=1;
        while k < 8
        T_k = norm( V* (Dprime^k) *V^-1 );        
        stuff(1,k) = k;
        stuff(2,k) = T_k;
        stuff(3,k) = e;
        k = k+1; %intervals
        end
        
        if e == a
            XXX = stuff(1,:);
            ZZZ = stuff(2,:);
            YYY = stuff(3,:);
        else
           XXX = [XXX, stuff(1,:)];
           ZZZ = [ZZZ, stuff(2,:)];
           YYY = [YYY, stuff(3,:)];
        end
    
    
   %     if e == a
    %        scatter3(XXX,ZZZ,YYY,'*')
     %       hold all
      %      line(XXX,ZZZ,YYY)
       % else 
        %    scatter3(XXX,ZZZ,YYY,'*')
         %   line(XXX,ZZZ,YYY)
       % end
        
        
end 

plot3(XXX,YYY,ZZZ,'.','markersize',12)
grid on
[xi,yi] = meshgrid(1:8, -3:3);
zi = griddata(XXX,YYY,ZZZ,xi,yi);
surf(xi,yi,zi)
xlabel('k');
ylabel('E');
zlabel('||A^k||');
[c,h] = contour(xi,yi,zi,16);
clabel(c,h)
xlabel('k');
ylabel('E');









   %Finding P(a) which is max of ||A^k|| vs k
    [V,D] = eig(T);
    %invv = inv(V);
    k = 1;
    
    %Find D' (eigen values > 1 = 0)
    [XX,YY] = size(D);
    i =1;
    while i <= XX
    j=1;
       while j <=YY
          Mag = ( (real( D(i,j) ))^2 + (imag( D(i,j) ))^2 )^(1/2);
          if Mag < 1.0
              Dprime(i,j) = D(i,j);
          else                  
              Dprime(i,j) = 0;
          end
          j = j+1;
        end
        i = i+1;
    end
     
    
    while k < 7
        T_k = norm( V* (Dprime^k) *V^-1 );        
        pl(k,1) = k;
        pl(k,2) = T_k;
        k = k+1;
    end
    
    %scatter( pl(:,1), pl(:,2) )
    P_a = max ( pl(:,2) ); 
        


%Attempt at Kreiss Constant
%N = size(T);
%z = 2;
%e =1;
%while abs(e) > 1.0
%    
%    x(z-1) = (z-1)* norm( (z*eye(N) - T)');
%    if z > 2
%        e = x(z-1) - x(z-2);
%    else
%        z=2;
%    end
    
%    z=z+1;
%end
%K = max(x);
%Top = (N*K*2.71828);
%fprintf("The Kreiss Constant is %e. The range should be [%e,%e]\n", K,K,Top)



    
    
    




%Build T(E) Imaginary 
    %for j = c:d
     % First = Hsdagger*((j*i)*I-Hd);
     % Second = -Hsdagger*Hs;
      %Third = I;
      %Fourth = zeros(m,n);
      %T = [First, Second; Third, Fourth];
      %opts.npts=100;
      %opts.levels=-7.5:0.5:-2;
      %opts.ax=[-1 0.2 -1.05 0];
      %opts.proj_lev=Inf;
      %opts.colour=1;
      %opts.thick_lines=1;
      %opts.scale_equal=1;
      %opts.grid=0;
      %opts.dim=1;
      %opts.no_ews=0;
      %opts.no_psa=0;
      %opts.fov=1;
      %opts.direct=1;
      %opts.print_plot=1;
      %eigtool(T, opts)
      %set(gcf,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
      %print(compose("Timag_%d.png",zz),'-dpng',['-r' num2str(resolution)]);
      %saveas(fig, compose("Timag_%",zz), 'png')
      %zz=zz+1;
    %end 
  
%end

