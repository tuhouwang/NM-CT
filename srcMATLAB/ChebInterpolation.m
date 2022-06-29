function [c, rho, alpha] = ChebInterpolation(dep, c, rho, alpha, Layers, Ns)

    for i = 1 : Layers
        x = cos((0 : Ns(i)) * pi / Ns(i))';
        z = ((dep{i}(1) + dep{i}(end)) / (dep{i}(end) - dep{i}(1)) - x) ...
           * (dep{i}(end) - dep{i}(1)) * 0.5;

        c(i)     = {interp1(dep{i}, c{i},     z, 'linear', 'extrap')};
        rho(i)   = {interp1(dep{i}, rho{i},   z, 'linear', 'extrap')};
        alpha(i) = {interp1(dep{i}, alpha{i}, z, 'linear', 'extrap')};
    end
    
end
