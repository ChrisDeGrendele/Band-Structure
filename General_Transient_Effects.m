
%A SIMPLE TWO BAND STRUCTURE
alpha = 0;
beta_1 = -.3;
beta_2 = -1;
Hd = [alpha, beta_1;beta_1,-alpha];
Hs = [0,beta_2;0,0];
E = 1;

T = buildT(Hd,Hs,E);

[V,D] = eig(T);


for i = 1:length(D)
    Mag = ( (real( D(i,i) ))^2 + (imag( D(i,i) ))^2 )^(1/2);
    if Mag > 1
        D(i,i) =0;
    end
end

Data_k = [];
Data_norm = [];

for k = 1:7
    Data_norm = [Data_norm,  norm( V* (D^k) *V^-1 )];
    Data_k = [Data_k, k];
end
plot(Data_k, Data_norm)
    





function T = buildT(Hd, Hs, e)
    num = size(Hs,1);
    
    A = [e*eye(num) - Hd, -Hs; eye(num), zeros(num)];
    B = [Hs', zeros(num); zeros(num), eye(num)];
    [V,D] = eig(A,B); %AV = BVD
    
    rank_HS = rank(Hs);
    D_array = diag(D)';
    [D_array, I] = sort(D_array);
    
    %sorts V correspoding to D
    temp = V;
    for i = 1:length(I)
        V(:,i) = temp(:,I(i));
    end
    
    index1 = size(D_array,2)/2 - rank_HS +1;
    index2 = size(D_array,2)/2 + rank_HS;
    
    %Sets all the data points not in 2r to 0
    for i = 1:length(D_array)
        if i < index1
            D_array(i) = 0;
        elseif i > index2
            D_array(i) = 0;
        end
    end
    
    T = V * diag(D_array) * inv(V);
    
end
    
    
    






