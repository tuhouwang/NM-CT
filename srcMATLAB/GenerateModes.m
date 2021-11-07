function [psi, psizs, z] = GenerateModes(eigvectorw, eigvectorb, nmodes, dz, ...
                                        zs, rhow, rhob, interface, bottom)

    zt1  = 0 : dz : interface;
    zt2  = interface : dz : bottom;
    z    = 0 : dz : bottom;

    xt1  = -2 / interface * zt1 + 1;
    xt2  = -2 / (bottom - interface) * zt2 + ...
           (bottom + interface) / (bottom - interface);

    psi1 = InvChebTrans(eigvectorw, xt1);
    psi2 = InvChebTrans(eigvectorb, xt2);
    
    if(zt1(end) == zt2(1))
        psi  = [psi1(1 : length(xt1) - 1, :); psi2];
    else
        zt   = [zt1; zt2];
        psi  = [psi1; psi2];
        psi  = interp1(zt,psi,z,'linear','extrap');
    end

    norm = Normalization(eigvectorw, eigvectorb, ...
                 nmodes, rhow, rhob, interface, bottom);

    psi   = psi * diag(1 ./ norm);
    psizs = interp1(z, psi, zs, 'linear');
    
end
