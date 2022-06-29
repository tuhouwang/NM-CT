function [nr, r, rhozs, k, kh] = Initialization(Layers, Ns, freq, ...
                 rmax, dr, zs, dep, c, rho, alpha, interface, ch, alphah)

    w  = 2 * pi * freq;
    r  = 0 : dr : rmax;
    nr = length(r);

    for i = 1 : Layers
        if( zs <= interface(i) )
            x = cos((0 : Ns(i)) * pi / Ns(i))';
            z = ((dep{i}(1) + dep{i}(end)) / (dep{i}(end) -  dep{i}(1)) ...
                      - x) * (dep{i}(end) - dep{i}(1)) * 0.5;
        
            rhozs = interp1(z, rho{i}, zs, 'linear');
            break
        end
    end
    
    k = cell(Layers, 1);  
    for i = 1 : Layers
        k(i) = {w ./ c{i} .* (1.0 + 1i * alpha{i} / (40.0 * pi * log10(exp(1.0))))};
    end

    kh = w / ch * (1.0 + 1i * alphah / (40.0 * pi * log10(exp(1.0))));

end
