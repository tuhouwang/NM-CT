function D  = DerivationMatrix(n)

    D = zeros(n, n);
    
    for k = 1 : n
        j = k + 1 : 2 : n;
        D(k, j) = 2 * j - 2;
    end
    
    D(1, :) = D(1, :) * 0.5;

end
