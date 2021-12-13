function [psi, psizs, z] = GenerateModes(eigvectorw, eigvectorb, nmodes, dz, ...
                  zs, rhow, rhob, rhoh, kr, kh, Lowerboundary, interface, Hb)

    zt1  = 0 : dz : interface;
    zt2  = interface : dz : Hb;
    z    = 0 : dz : Hb;

    xt1  = -2 / interface * zt1 + 1;
    xt2  = -2 / (Hb - interface) * zt2 + ...
           (Hb + interface) / (Hb - interface);

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
                 nmodes, rhow, rhob, interface, Hb);
             
    if(Lowerboundary == 'A')
        norm = norm + 0.5 ./ rhoh * psi(end, :) .^ 2  ./ sqrt(kr .^ 2 - kh ^ 2).';
    end

    psi   = psi * diag( sqrt(1 ./ norm) );
    psizs = interp1(z, psi, zs, 'linear');
    
end
