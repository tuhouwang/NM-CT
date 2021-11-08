function [casename, Nw, Nb, cpmax, freq, zs, zr, rmax, dr, interface, ...
    Hb, dz, Lowerboundary, tlmin, tlmax, depw, cw, rhow, alphaw, ...
    depb, cb, rhob, alphab, ch, rhoh, alphah] = ReadEnvParameter(env_file)

    fid           = fopen(env_file);
    casename      = fgetl(fid);
    Nw            = fscanf(fid, '%d', 1);
    Nb            = fscanf(fid, '%d', 1);
    cpmax         = fscanf(fid, '%f', 1);
    freq          = fscanf(fid, '%f', 1);
    zs            = fscanf(fid, '%f', 1);
    zr            = fscanf(fid, '%f', 1);
    rmax          = fscanf(fid, '%f', 1);
    dr            = fscanf(fid, '%f', 1);
    interface     = fscanf(fid, '%f', 1);
    Hb            = fscanf(fid, '%f', 1);
    dz            = fscanf(fid, '%f', 1);
    Lowerboundary = fscanf(fid, '%s', 1);
    tlmin         = fscanf(fid, '%f', 1);
    tlmax         = fscanf(fid, '%f', 1);
    nw            = fscanf(fid, '%d', 1);
    nb            = fscanf(fid, '%d', 1);
    
    if (interface > 0.0 && interface < Hb && Nw > 2 && Nb > 2)
        WProfile  = fscanf(fid, '%f %f', [4, nw]);
        depw      = WProfile(1, 1:nw);
        cw        = WProfile(2, 1:nw);
        rhow      = WProfile(3, 1:nw);
        alphaw    = WProfile(4, 1:nw);
        BProfile  = fscanf(fid, '%f %f', [4, nb]);
        depb      = BProfile(1, 1:nb);
        cb        = BProfile(2, 1:nb);
        rhob      = BProfile(3, 1:nb);
        alphab    = BProfile(4, 1:nb);
    else
        error('Error! h must greater than 0 and less than H!');
    end
    
    if (Lowerboundary ~= 'V' && Lowerboundary ~= 'R' && Lowerboundary ~= 'A')
        disp('Error! The lower boundary must be vaccum, rigid or halfspace!');
        error('Please input 0 or 1 in Lowerboundary!');
    end
    
    if (Lowerboundary == 'A')
        halfspace = fscanf(fid, '%f', 4);
        ch     = halfspace(2);
        rhoh   = halfspace(3);
        alphah = halfspace(4);
    else
        ch     = 0.0;
        rhoh   = 0.0;
        alphah = 0.0;
    end
    
    % Check the input underwater sound profile    
    if (Nw < 2 || Nb < 2)
        error('Nw and Nb must greater than 2!');
    end
    
    if (depw(1) ~= 0.0 || depw(nw) ~= interface ||...
            depb(1) ~= interface || depb(nb) ~= Hb)
        error('Error! input sound profile is unsuitable!');
    end

    if ((rmax / dr - floor(rmax / dr)) ~= 0)
        error('Please reinput the dr and rmax !');
    end

    if (zs <= 0  || zs >= Hb || zr <= 0 || zr >= Hb)
        error('zs and zr must be greater than 0 and less than H!');
    end
    
    if (tlmin >= tlmax)
        error('tlmin must less than tlmax!');
    end    

    fclose(fid);

end
