function a = Normalization(eigvectorw, eigvectorb, nmodes, ...
                           rhow, rhob, interface, bottom)

    Nw = size(eigvectorw, 1) - 1;
    Nb = size(eigvectorb, 1) - 1;
    Rw = ConvolutionMatrix(ChebTransFFT(Nw, 1./rhow));
    Rb = ConvolutionMatrix(ChebTransFFT(Nb, 1./rhob));

    P1      = zeros(1, Nw+1);
    P2      = zeros(1, Nb+1);

    k       = 0 : 2 : Nw;
    P1(k+1) = -2 ./ (k .^ 2 - 1);
    k       = 0 : 2 : Nb;
    P2(k+1) = -2 ./ (k .^ 2 - 1);
    
    f1      = zeros(Nw + 1, nmodes);
    f2      = zeros(Nb + 1, nmodes);
    
    for j = 1 : nmodes
        f1(:, j) = ConvolutionMatrix(eigvectorw(:, j)) * eigvectorw(:, j);
        f2(:, j) = ConvolutionMatrix(eigvectorb(:, j)) * eigvectorb(:, j);
    end
    
    a = sqrt( P1 * Rw * f1 * interface / 2 + ...
        P2 * Rb * f2 * (bottom - interface) / 2 );

end
