function [kr, eigvector] = EigenValueVector(Ns, Layers, dep, k, rho, kh, ...
                                            rhoh, Lowerboundary)

    D = cell(Layers, 1);
    U = zeros(sum(Ns+1), sum(Ns+1));

    if(Lowerboundary == 'A')
        n = 1;
        for i = 1 : Layers
            D(i) = {DerivationMatrix(Ns(i) + 1)};

            A = 4.0 / (dep{i}(end) - dep{i}(1) ) ^ 2 * ConvolutionMatrix(ChebTransFFT(Ns(i), rho{i})) * ...
                D{i} * ConvolutionMatrix(ChebTransFFT(Ns(i), 1.0 ./ rho{i})) * D{i} + ...
                ConvolutionMatrix(ChebTransFFT(Ns(i), k{i} .^ 2));

            U(n:n+Ns(i), n:n+Ns(i)) = A - kh ^ 2 * eye(Ns(i)+1);
            n = n + Ns(i) + 1;
        end

        V = zeros(sum(Ns+1));
        W = eye  (sum(Ns+1));
        
        % boundary condition
        n = 1;
        for i = 1 : Layers - 1
            % sound pressure is continuous
            U(n+Ns(i)-1,         n:n+Ns(i)          ) = (-1.0).^(0:Ns(i));
            U(n+Ns(i)-1, n+Ns(i)+1:n+Ns(i)+Ns(i+1)+1) = -1.0;
            W(n+Ns(i)-1,                 n+Ns(i)-1  ) =  0;
            
            Pu = -1 / rho{i}(end) / (dep{i}(end) - dep{i}(1)) * ((-1.0).^(0 : Ns(i)))  * D{i};
            Pd =  1 / rho{i+1}(1) / (dep{i+1}(end) - dep{i+1}(1)) * ones(1, Ns(i+1)+1) * D{i+1};
            % normal velocity is continuous
            U(n+Ns(i),         n:n+Ns(i)          ) = Pu;
            U(n+Ns(i), n+Ns(i)+1:n+Ns(i)+Ns(i+1)+1) = Pd;
            W(n+Ns(i),           n+Ns(i)          ) = 0;
            
            n = n + Ns(i) + 1;
        end

        % upper boundary
        U(end-1, 1:Ns(1)+1) = 1.0;
        W(end-1, end-1)     = 0.0;
        % the most important lower boundary
        bottom = (-1.0) .^ (0:Ns(end));
        U(end, end-Ns(end):end) = -2i * rhoh / rho{end}(end) / (dep{end}(end) - dep{end}(1)) * bottom * D{end};
        V(end, end-Ns(end):end) = bottom;
        W(end, end)             = 0;

        [v, kz] = polyeig(U, V, W);
        kr  = sqrt(kh ^ 2 - kz .^ 2);
        ind = find(real(kz) >= 0 & abs(kr) < max(abs(k{1})));
        kr  = kr(  ind);
        v   = v (:,ind);
        [~, ind] = sort(real(kr), 'descend');
        kr  = kr(  ind);
        v   = v (:,ind);

        eigvector = cell(Layers, 1);
        n = 1;
        for i = 1 : Layers
            eigvector(i) = {v(n:n+Ns(i), :)};
            n = n + Ns(i) + 1;
        end
        
    else
        % ideal seafloor
        n = 1;
        for i = 1 : Layers
            D(i) = {DerivationMatrix(Ns(i) + 1)};

            A = 4.0 / ( dep{i}(end) - dep{i}(1) ) ^ 2 * ConvolutionMatrix(ChebTransFFT(Ns(i), rho{i})) * ...
                D{i} * ConvolutionMatrix(ChebTransFFT(Ns(i), 1.0 ./ rho{i})) * D{i} + ...
                ConvolutionMatrix( ChebTransFFT(Ns(i), k{i} .^ 2));

            U(n:n+Ns(i)-2,               n:n+Ns(i)-2    ) = A(1:Ns(i)-1,     1:Ns(i)-1);
            U(n:n+Ns(i)-2, sum(Ns-1)+2*i-1:sum(Ns-1)+2*i) = A(1:Ns(i)-1, Ns(i):Ns(i)+1);
            n = n + Ns(i) - 1;
        end
        % boundary condition
        n = 1;
        for i = 1 : Layers - 1
            Pu = -1 / rho{i}(end) / (dep{i}(end) - dep{i}(1)) * ((-1.0).^(0 : Ns(i)))  * D{i};
            Pd =  1 / rho{i+1}(1) / (dep{i+1}(end) - dep{i+1}(1)) * ones(1, Ns(i+1)+1) * D{i+1};
            % sound pressure is continuous
            U(sum(Ns-1)+2*i-1,               n:n+Ns(i)-2    ) = (-1.0).^(0:Ns(i)-2);
            U(sum(Ns-1)+2*i-1, sum(Ns-1)+2*i-1:sum(Ns-1)+2*i) = (-1.0).^(Ns(i)-1:Ns(i));
            % normal velocity is continuous
            U(sum(Ns-1)+2*i,               n:n+Ns(i)-2    ) = Pu(1    :Ns(i)-1);
            U(sum(Ns-1)+2*i, sum(Ns-1)+2*i-1:sum(Ns-1)+2*i) = Pu(Ns(i):Ns(i)+1);

            n = n + Ns(i) - 1;
            % sound pressure is continuous
            U(sum(Ns-1)+2*i-1,               n:n+Ns(i+1)-2    ) = -1;
            U(sum(Ns-1)+2*i-1, sum(Ns-1)+2*i+1:sum(Ns-1)+2*i+2) = -1;
            % normal velocity is continuous
            U(sum(Ns-1)+2*i,               n:n+Ns(i+1)-2    ) = Pd(1      :Ns(i+1)-1);
            U(sum(Ns-1)+2*i, sum(Ns-1)+2*i+1:sum(Ns-1)+2*i+2) = Pd(Ns(i+1):Ns(i+1)+1);
        end

        % upper boundary, pressure-free boundary
        U(end-1, 1          :    Ns(1)-1) = 1.0;
        U(end-1, sum(Ns-1)+1:sum(Ns-1)+2) = 1.0;
        % lower boundary, perfectly free / rigid
        Low = (-1.0) .^ (0 : Ns(end));
        if(Lowerboundary == 'R')
            Low = Low * D{end};
        end

        U(end, sum(Ns-1)-Ns(end)+2:sum(Ns-1)) = Low(1      :Ns(end)-1);
        U(end,               end-1:end)       = Low(Ns(end):Ns(end)+1);

        % blocking
        L11 = U(1          :sum(Ns-1), 1          :sum(Ns-1));
        L12 = U(1          :sum(Ns-1), sum(Ns-1)+1:sum(Ns+1));
        L21 = U(sum(Ns-1)+1:sum(Ns+1), 1          :sum(Ns-1));
        L22 = U(sum(Ns-1)+1:sum(Ns+1), sum(Ns-1)+1:sum(Ns+1));

        L = L11 - L12 * (L22 \ L21);
        [v, k2] = eig(L);

        v2 = - (L22 \ L21) * v;

        eigvector = cell(Layers, 1);

        n = 1;
        for i = 1 : Layers
            eigvector(i) = {[v(n:n+Ns(i)-2, :); v2(2*i-1:2*i,:)]};
            n = n + Ns(i) - 1;
        end

        kr = sqrt(diag(k2));
        [~, ind] = sort(real(kr), 'descend');
        kr       = kr(ind);
        for i = 1 : Layers
            eigvector(i) = {eigvector{i}(:, ind)};
        end        
    end
end
