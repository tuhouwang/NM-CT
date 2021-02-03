function [casename, Nw, Nb, cpmax, freq, zs, zr, rmax, dr, interface, ...
    bottom, dz, Lowerboundary, tlmin, tlmax, depw, cw, rhow, alphaw, ...
    depb, cb, rhob, alphab] = ReadEnvParameter(env_file)

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
    bottom        = fscanf(fid, '%f', 1);
    dz            = fscanf(fid, '%f', 1);
    Lowerboundary = fscanf(fid, '%d', 1);
    tlmin         = fscanf(fid, '%f', 1);
    tlmax         = fscanf(fid, '%f', 1);
    nw            = fscanf(fid, '%d', 1);
    nb            = fscanf(fid, '%d', 1);
    
    if (interface > 0.0 && interface < bottom && Nw > 2 && Nb > 2)
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
        error('Error! h must greater than 0 and less than H !');
    end
    
    % Check the input underwater sound profile
    if (depw(1) ~= 0.0 || depw(nw) ~= interface ||...
            depb(1) ~= interface || depb(nb) ~= bottom)
        error('Error! input sound profile is unsuitable !');
    end

    if ((interface / dz - floor(interface / dz)) ~=0 ||...
            (bottom / dz - floor(bottom / dz)) ~=0)
        error('Error! The input dz unsuitable !');
    end

    if ((rmax / dr - floor(rmax / dr)) ~= 0)
        error('Please reinput the dr and rmax !');
    end

    if (Lowerboundary ~= 0 && Lowerboundary ~= 1)
        disp('Error! The lower boundary must be rigid or soft !');
        error('Please input 0 or 1 in Lowerboundary !');
    end

    if (zs <= 0  || zs >= bottom || zr <= 0 || zr >= bottom)
        error('zs and zr must be greater than 0 and less than H !');
    end
    
    if (tlmin >= tlmax)
        error('tlmin must less than tlmax !');
    end    

    fclose(fid);

end
