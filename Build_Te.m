function [Te] = Build_Te(Hd, Hs, E)

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

        First = Hsdagger*((E+.0001i)*I-Hd);
        Second = -Hsdagger*Hs;
        Third = I;
        Fourth = zeros(m,n);
        Te = [First, Second; Third, Fourth];
   
end
