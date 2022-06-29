function a = Normalization(eigvector, nmodes, rho, dep, Layers)

    a = zeros(1, nmodes);
    for i = 1 : Layers
        N = size(eigvector{i}, 1) - 1;
        R = ConvolutionMatrix(ChebTransFFT(N, 1 ./ rho{i}));
        P      = zeros(1, N+1);
        k      = 0 : 2 : N;
        P(k+1) = -2 ./ (k .^ 2 - 1);
        f      = zeros(N + 1, nmodes);

        for j = 1 : nmodes
            f(:, j) = ConvolutionMatrix(eigvector{i}(:, j)) * eigvector{i}(:, j);
        end

        a = a + P * R * f * (dep{i}(end) - dep{i}(1)) * 0.5;
    end
end
