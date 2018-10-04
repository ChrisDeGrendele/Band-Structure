%Data
Hd = [.8,-.4;-.4,-.8];
Hs = [0,-1.3;-1.3,0];
E = [-3,3];

%Hd = [0];
%Hs = [-1/2];
%E=[-2,2];
    a = E(1,1);
    b = E(1,2);
  
    
   %%%%%  ||A^k|| vs. E vs. k.  %%%%%%
for e = a:b 
    Te = Build_Te(Hd,Hs,e);
    Build_Pseudo(Te,e)
    [V,D] = eig(Te);
    
   
%Build D' (eigen values > 1 = 0)
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
         %This is the plot of ||A^k|| vs E vs k
        while k < 20
        T_k = norm( V* (Dprime^k) *V^-1 );
                if (e == a) & (k==1)
                    XXX = k;
                    YYY = e;
                    ZZZ = T_k;                    
                else
                    XXX = [XXX,k];
                    YYY = [YYY,e];
                    ZZZ = [ZZZ,T_k];                    
                end
        k = k+1; 
        end
        
end

%Visulizations:
figure
plot3(XXX,YYY,ZZZ,'.','markersize',12)
grid on
[xi,yi] = meshgrid(1:20, a:b);
zi = griddata(XXX,YYY,ZZZ,xi,yi);
figure
surf(xi,yi,zi)
xlabel('k');
ylabel('E');
zlabel('||A^k||');
figure
[c,h] = contour(xi,yi,zi,16); %last number is number of contour lines
clabel(c,h)
xlabel('k');
ylabel('E');







   %Finding P(a) which is max of ||A^k|| vs k
   %Note this is not addapted specifically for different E values yet. Not 
   %yet needed. This will run for the final e value provided above.
    [V,D] = eig(Te);
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

