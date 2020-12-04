function T = ChebPoly(N, x)
%   Chebyshev polynomials of degree N.
%   T = ChebPoly(N, x) returns the Chebyshev polynomials for degree k=0:N 
%   at points x. The i-th row and j-th column of T stores T_j(x_i).
   
    T = zeros(length(x), N+1);
    for k = 0 : N
        T(:, k+1) = cos( k * acos(x) );
    end

end