function [nr, r, rhozs, kw, kb, w] = Initialization(Nw, Nb, freq, rmax, ...
            dr, zs, rhow, rhob, cw, cb, alphaw, alphab, interface, bottom)

    w       = 2 * pi * freq;
    r       = dr : dr : rmax;
    nr      = length(r);

    if(zs <= interface)
        x1     = cos( (0 : Nw) * pi / Nw)';
        z1     = (1.0 - x1) * interface / 2;
        rhozs  = interp1(z1, rhow, zs, 'linear');
    else
        x2     = cos( (0 : Nb) * pi / Nb)';
        z2     = (1.0 - x2) * (bottom - interface) / 2 + interface;
        rhozs  = interp1(z2, rhob, zs, 'linear');
    end

    kw = w ./ cw .* (1.0 + 1i * alphaw / (40.0 * pi * log10(exp(1.0))));
    kb = w ./ cb .* (1.0 + 1i * alphab / (40.0 * pi * log10(exp(1.0))));

end
