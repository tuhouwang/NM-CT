function fk = ChebTransFFT(N, fx)

%  The function computes the Chebyshev tranforms by FFT in the nodes:
%  x_j = cos(j pi/N) from Physical space to Spectral space.
%  
%  Input:
%  N:           Degree of polynomials
%  fx:          vector to be transformed

%  Output:
%  fk:          Chebyshev coefficients of fx
%--------------------------------------------------------------------------

    fk = fft([fx; flipud(fx(2:N))]);              % Extend and compute fft
    fk = fk(1:N+1).*[0.5; ones(N-1,1); 0.5]/N;    % fk contains Chebyshev
                                                  % coefficients of fx
end
