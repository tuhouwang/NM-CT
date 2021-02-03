function [kr, eigvectorw, eigvectorb] = EigenValueVector(Nw, Nb, ...
                interface, bottom, kw, kb, rhow, rhob, Lowerboundary)

    D1  = DerivationMatrix(Nw + 1);
    D2  = DerivationMatrix(Nb + 1);

    A = 4.0 / interface ^ 2 * ConvolutionMatrix(ChebTransFFT(Nw, rhow)) * ...
        D1 * ConvolutionMatrix(ChebTransFFT(Nw, 1.0 ./ rhow)) * D1 +...
        ConvolutionMatrix( ChebTransFFT(Nw, kw .^ 2));

    B = 4.0 / (bottom-interface) ^ 2 * ConvolutionMatrix(ChebTransFFT(Nb, rhob)) * ...
        D2 * ConvolutionMatrix(ChebTransFFT(Nb, 1.0 ./ rhob)) * D2 +...
        ConvolutionMatrix( ChebTransFFT(Nb, kb .^ 2) );

    U = zeros(Nw+Nb+2, Nw+Nb+2);
    U(1 :Nw-1,    1      :Nw-1)    = A(1:Nw-1, 1 :Nw-1);
    U(1 :Nw-1,    Nw+Nb-1:Nw+Nb)   = A(1:Nw-1, Nw:Nw+1);
    U(Nw:Nw+Nb-2, Nw     :Nw+Nb-2) = B(1:Nb-1, 1 :Nb-1);
    U(Nw:Nw+Nb-2, Nw+Nb+1:Nw+Nb+2) = B(1:Nb-1, Nb:Nb+1);
    
    %upper boundary
    U(Nw+Nb-1, 1      :Nw-1 ) = 1.0;
    U(Nw+Nb-1, Nw+Nb-1:Nw+Nb) = 1.0;
    
    %lower boundary
    bottom_boundary = (-1.0) .^ (0:Nb);
    if(Lowerboundary == 1)
        bottom_boundary = bottom_boundary * D2;
    end
    U(Nw+Nb+2, Nw     :Nw+Nb-2) = bottom_boundary(1 :Nb-1);
    U(Nw+Nb+2, Nw+Nb+1:Nw+Nb+2) = bottom_boundary(Nb:Nb+1);
    
    %first interface boundary
    U(Nw+Nb, 1      :Nw-1   ) = (-1.0).^(0 : Nw-2);
    U(Nw+Nb, Nw     :Nw+Nb-2) = -1.0;
    U(Nw+Nb, Nw+Nb-1:Nw+Nb  ) = (-1.0).^(Nw-1 : Nw);
    U(Nw+Nb, Nw+Nb+1:Nw+Nb+2) = -1.0;
    
    %second interface boundary
    Pu =  1 / rhow(Nw+1) / interface * ((-1.0) .^ (0 : Nw)) * D1;
    Pd = -1 / rhob(1)    / (bottom - interface) * ones(1, Nb+1) * D2;

    U(Nw+Nb+1, 1      :Nw-1   ) = Pu(1  : Nw-1);
    U(Nw+Nb+1, Nw     :Nw+Nb-2) = Pd(1  : Nb-1);
    U(Nw+Nb+1, Nw+Nb-1:Nw+Nb  ) = Pu(Nw : Nw+1);
    U(Nw+Nb+1, Nw+Nb+1:Nw+Nb+2) = Pd(Nb : Nb+1);
    
    %blocking
    L11 = U(1       : Nw+Nb-2, 1       : Nw+Nb-2);
    L12 = U(1       : Nw+Nb-2, Nw+Nb-1 : Nw+Nb+2);
    L21 = U(Nw+Nb-1 : Nw+Nb+2, 1       : Nw+Nb-2);
    L22 = U(Nw+Nb-1 : Nw+Nb+2, Nw+Nb-1 : Nw+Nb+2);

    L = L11 - L12 * (L22 \ L21);
    [v, k2] = eig(L);

    v2 = - (L22 \ L21) * v;

    eigvectorw = [v(1  : Nw-1,    :); v2(1:2, :)];
    eigvectorb = [v(Nw : Nw+Nb-2, :); v2(3:4, :)];

    kr       = sqrt(diag(k2));
    [~, ind] = sort(real(kr), 'descend');
    kr       = kr(ind);

    eigvectorw = eigvectorw(:, ind);
    eigvectorb = eigvectorb(:, ind);

end
