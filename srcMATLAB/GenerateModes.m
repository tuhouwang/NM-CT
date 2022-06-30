function [psi, psizs, z] = GenerateModes(eigvector, nmodes, dz, zs, ...
                   rho, rhoh, kr, kh, Lowerboundary, dep, Layers)

    z   = 0 : dz : dep{end}(end);    
    zm  = [];
    psi = [];
    
    for i = 1 : Layers
        zi = dep{i}(1) : dz : dep{i}(end);
        xi = -2 / (dep{i}(end) - dep{i}(1)) * zi + ...
                  (dep{i}(end) + dep{i}(1)) / (dep{i}(end) - dep{i}(1));
        if(isempty(zm) == 1)
           zm = zi;
           psi= InvChebTrans(eigvector{i}, xi);
        elseif(zm(end) == zi(1))   
           zm  = [zm, zi(2:end)];
           p   = InvChebTrans(eigvector{i}, xi);
           psi = [psi; p(2:end,:)];             
        else
           zm  = [zm, zi];
           psi = [psi; InvChebTrans(eigvector{i}, xi)];          
        end  
    end
    
    psi  = interp1(zm, psi, z, 'linear', 'extrap');
    norm = Normalization(eigvector, nmodes, rho, dep, Layers);
    
    if(Lowerboundary == 'A')
        norm = norm + 0.5 ./ rhoh * psi(end, :) .^ 2  ./ sqrt(kr .^ 2 - kh ^ 2).';
    end         

    psi   = psi * diag(sqrt(1 ./ norm));
    psizs = interp1(z, psi, zs, 'linear', 'extrap');
    
end
