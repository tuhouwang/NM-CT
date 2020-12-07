function C  = ConvolutionMatrix(v)

    n = length(v);
    C = zeros(n, n);
    
    for i = 1 : n
        for k = 1 : n
            
            j = k - i + 1;       
            if (j <= n && j >= 1)
                C(k,i) = C(k,i) + v(j);
            end
            
            j = i - k + 1;
            if (j <= n && j >= 1)
                C(k,i) = C(k,i) + v(j);
            end
            
            j = i + k - 1 ;
            if (j <= n && j > i)
                C(k,i) = C(k,i) + v(j);
            end  
            
            C(k,i) = C(k,i) * 0.5; 
        end
    end

end
