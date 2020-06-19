close all
clear all

%======================TRIAL SYSTEMS======================================

%asymmetric trial, matt email june 10th 2020
%energy should be -10,10
%alpha = 0.3;
%beta_1 = -0.3;
%beta_2 = -1.0;
%alpha_core = 0.5;
%beta_core = -3.2;
%beta_3 = 0.5;

%H1:::::
%Hd  = [alpha, beta_1, 0;
%       beta_1, alpha,  0;
%       0,     0,       alpha_core];

%H2::::: 
%Hd = [alpha,  beta_1, 0;
%      beta_1, alpha,  beta_3;
%      0,      beta_3  beta_core];

%Hs = [0, beta_2, 0;
%     0, 0,      0;
%      0, 0,      beta_core];



%%%%%%%%%%%%%%%%%%%%%%%%


% Asymmetric Band test with ABA like SSH model
%alpha = 0;
%beta_1 = -3;
%beta_2 = -1;

%alpha_p = 5;
%beta_1_p = 2;
%beta_2_p = 3;

%H_SSH = [alpha, beta_1;beta_1,-alpha];
%V = [0,beta_2;0,0];

%H_SSH_p = [alpha_p, beta_1_p;beta_1_p,-alpha_p];
%V_p = [0,beta_2_p;0,0];

%Hd = [H_SSH, conj(V)'; V, H_SSH_p];
%Hs = [zeros(2), V_p; zeros(2), zeros(2)];





%%%%%%IDK what this one is below%%%%%%
%alpha = 0.9;
%beta_sig = -0.25;
%beta_pi = 1.2;
%Hd = [alpha, 0, 0; 0, alpha, 0; 0, 0, alpha];
%Hs = [beta_sig, 0, 0; 0, beta_pi, 0; 0, 0, beta_pi];

%beta_pi2 = 0.6;
%Hd = [alpha, 0, 0; 0, alpha, 0; 0, 0, alpha];
%Hd = [Hd, zeros(3);zeros(3), Hd];
%Hs = [beta_sig, 0, 0; 0, beta_pi, 0; 0, 0, beta_pi2];
%Hs2 = [beta_sig, 0, 0; 0, 0.4*beta_pi, 0; 0, 0, 0.3*beta_pi2];
%Hs = [zeros(3),Hs;zeros(3),Hs2];

%======================OTHER SYSTEMS======================================

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


%Model System 6
%Hd = -.5* [0,1,0,0,0,1;1,0,1,0,0,0;0,1,0,1,0,0;0,0,1,0,1,0;0,0,0,1,0,1;1,0,0,0,1,0];
%Hs = zeros(6);
%Hs(2,1) = -.5;
%Hs(3,4) = -.5;
%Hs(5,4)= -.5;
%Hs(6,1) = -.5;

%alpha = 0.8;
%beta_1 = -.4;
%beta_2 = -1.3;

%Test system 1
%Hd = [alpha, beta_1;beta_1,-alpha];
%Hs = [0,beta_2;beta_2,0];


%======================END OF OTHER SYSTEMS===============================


%======================PAPER SYSTEMS======================================

%SSH "Simple Example"
alpha = 0;
beta_1 = -0.3;
beta_2 = -1.0;
systems(1).Hd = [alpha, beta_1;beta_1,-alpha];
systems(1).Hs = [0,beta_2;0,0];
systems(1).E  =[-8,8];
systems(1).name = "simple2bandex.mat";


% PARA
alpha = 0;
beta = -.5;
systems(2).Hd = [alpha, beta, 0, 0, 0, beta ;beta,alpha,beta, 0, 0, 0; 0, ...
                beta, alpha, beta, 0, 0; 0, 0, beta, alpha, beta, 0; 0, 0, 0, ...
                beta, alpha, beta; beta, 0, 0, 0, beta, alpha];
systems(2).Hs = zeros(6);
systems(2).Hs(1,4) = beta;
systems(2).E = [-4,4];
systems(2).name = "Para.mat";


% META
alpha = 0;
beta = -.5;
systems(3).Hd = [alpha, beta, 0, 0, 0, beta ;beta,alpha,beta, 0, 0, 0; 0, ...
                beta, alpha, beta, 0, 0; 0, 0, beta, alpha, beta, 0; 0, 0, 0, ...
                beta, alpha, beta; beta, 0, 0, 0, beta, alpha];
systems(3).Hs = zeros(6);
systems(3).Hs(1,3) = beta;
systems(3).E = [-4,4];
systems(3).name = "Meta.mat";


% ORTHO
alpha = 0;
beta = -.5;
systems(4).Hd = [alpha, beta, 0, 0, 0, beta ;beta,alpha,beta, 0, 0, 0; 0, ...
                beta, alpha, beta, 0, 0; 0, 0, beta, alpha, beta, 0; 0, 0, 0, ...
                beta, alpha, beta; beta, 0, 0, 0, beta, alpha];
systems(4).Hs = zeros(6);
systems(4).Hs(1,2) = beta;
systems(4).E = [-4,4];
systems(4).name = "Ortho.mat";


%Asymetric Complex Bands 
systems(5).Hd = zeros(10);
systems(5).Hd(1,2)= -.5;
systems(5).Hd(2,1) = -.5;
systems(5).Hd(1,9)= -.5;
systems(5).Hd(9,1) = -.5;
systems(5).Hd(2,3)= -.5;
systems(5).Hd(3,2) = -.5;
systems(5).Hd(3,10)= -.5;
systems(5).Hd(10,3) = -.5;
systems(5).Hd(4,5)= -.5;
systems(5).Hd(5,4) = -.5;
systems(5).Hd(4,10)= -.5;
systems(5).Hd(10,4) = -.5;
systems(5).Hd(5,6)= -.5;
systems(5).Hd(6,5) = -.5;
systems(5).Hd(6,7)= -.5;
systems(5).Hd(7,6) = -.5;
systems(5).Hd(7,8)= -.5;
systems(5).Hd(8,7) = -.5;
systems(5).Hd(8,9)= -.5;
systems(5).Hd(9,8) = -.5;
systems(5).Hd(9,10)= -.5;
systems(5).Hd(10,9) = -.5;
systems(5).Hs = zeros(10);
systems(5).Hs(4,7) = -.5;
systems(5).E = [-8,8];
systems(5).name= "Assymetric.mat";


%Comb Molecule
systems(6).Hd = [0,-0.5;-0.5,0];
systems(6).Hs = [-0.5, 0;0,0];
systems(6).E = [-4,4];
systems(6).name = "CombMolecule.mat";


% Siloxane
systems(7).Hd = zeros(8);
systems(7).Hd(1,1) = 4.08;
systems(7).Hd(1,5) = -3;
systems(7).Hd(1,8) = 5.4;
systems(7).Hd(2,2) = 10;
systems(7).Hd(2,6) = -1.4;
systems(7).Hd(3,3) = 10;
systems(7).Hd(3,7) = -1.4;
systems(7).Hd(4,4) = 10;
systems(7).Hd(4,5) = 5;
systems(7).Hd(4,8) = 6;
systems(7).Hd(5,1) = -3;
systems(7).Hd(5,4) = 5;
systems(7).Hd(5,5) = -1;
systems(7).Hd(6,2) = -1.4;
systems(7).Hd(6,6) = -16;
systems(7).Hd(7,3) = -1.4;
systems(7).Hd(7,7) = -16;
systems(7).Hd(8,1) = 5.4;
systems(7).Hd(8,4) = 6;
systems(7).Hd(8,8) = -16;
systems(7).Hs = zeros(8);
systems(7).Hs(1,4) = -3;
systems(7).Hs(1,8) = 5;
systems(7).Hs(2,6) = -1.4;
systems(7).Hs(3,7) = -1.4;
systems(7).Hs(4,5) = 5.4;
systems(7).Hs(4,8) = 6;
systems(7).Hs(5,5) = -0.15;
systems(7).Hs(5,8) = .35;
systems(7).Hs(6,6) = -0.13;
systems(7).Hs(7,7) = -0.13;
systems(7).Hs(8,5)= .35;
systems(7).Hs(8,8) = .56;
systems(7).E = [-50,50];
systems(7).name = "siloxane2.mat";

%Vertical Band 1
alpha = 0;
beta_1 = -0.3;
beta_2 = -1.0;
epsilon = 0.1;
systems(8).Hd = [alpha, beta_1; beta_1, alpha];
systems(8).Hs = [0, beta_2; epsilon, 0];
systems(8).E = [-5,5];
systems(8).name = "VerticalBand.1.mat";

%Vertical Band 2
alpha = 0;
beta_1 = -0.3;
beta_2 = -1.0;
epsilon = 0.01;
systems(9).Hd = [alpha, beta_1; beta_1, alpha];
systems(9).Hs = [0, beta_2; epsilon, 0];
systems(9).E = [-5,5];
systems(9).name = "VerticalBand.01.mat";

%Vertical Band 3
alpha = 0;
beta_1 = -0.3;
beta_2 = -1.0;
epsilon = 0.001;
systems(10).Hd = [alpha, beta_1; beta_1, alpha];
systems(10).Hs = [0, beta_2; epsilon, 0];
systems(10).E = [-5,5];
systems(10).name = "VerticalBand.001.mat";

%Vertical Band 4
alpha = 0;
beta_1 = -0.3;
beta_2 = -1.0;
epsilon = 0.0001;
systems(11).Hd = [alpha, beta_1; beta_1, alpha];
systems(11).Hs = [0, beta_2; epsilon, 0];
systems(11).E = [-5,5];
systems(11).name = "VerticalBand.0001.mat";




%======================END OF PAPER SYSTEMS===============================

for system = systems
    simulation(system)
end


function simulation(system, dosave, doplot, dobranchpoint)

    save_directory = "/Users/chris/Documents/GitHub/Band-Structure/newdata/";
    min_epsilon = 5e-3;
    ideal_spacing = 4e-4;
    
    
    if nargin == 1
        dosave = true;
        doplot = false;
        dobranchpoint = false;
    elseif nargin ~= 4
        error("Define all the args: dosave, doplot, dobranch, or none of them.")
    end

    %a = E(1,1);
    %b = E(1,2);


    
    %Data x -> Pure k values
    %Data y -> Pure E values
    [Data_x, Data_y] = amr_KvsE(system.Hd, system.Hs, system.E, min_epsilon, ideal_spacing);

    Imag_k = [];
    Imag_E = [];
    Real_k = [];
    Real_E = [];

    for e = 1:size(Data_x,2)
        if abs(imag( Data_x(e) )) > 1.e-6
            Imag_k = [Imag_k, imag(Data_x(e)) ];
            Imag_E = [Imag_E, Data_y(e)];
        else
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


    if doplot
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
        hold off
    end


    if dosave
        save(save_directory + system.name,'Imag_E','Imag_k', 'Real_E', 'Real_k');
    end
    
    %Below is stuff for calculating branch point
    %--------------------------------------------------------------------------
    if dobranchpoint

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
        while stop == 0
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
    end
end
