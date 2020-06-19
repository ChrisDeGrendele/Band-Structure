%Eigen Vector Test 2/24/20 w/ Josh

%parameters
alpha = .8;
beta_1 = -.4;
beta_2 = -1.3;
Hd = [alpha, beta_1;beta_1,-alpha];
Hs = [0,beta_2;beta_2,0];
E_value = -4;

%Traditional Calculation
Te = Build_Te(Hd, Hs, E_value);
[V_orig, D_orig] = eig(Te);

%Generalized Calculation
Hsdagger = (Hs)'
num = size(Hsdagger,1);
A = [e*eye(num) - Hd, -Hs; eye(num), zeros(num)];
B = [Hsdagger, zeros(num); zeros(num), eye(num)];
[V_gen,D_gen] = eig(A,B);

V_orig, V_gen

D_orig, D_gen