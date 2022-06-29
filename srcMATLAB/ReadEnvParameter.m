function [casename,Layers, Ns, cpmax, freq, zs, zr, rmax, dr, interface, dz, tlmin, tlmax, ...
          dep, c, rho, alpha, ch, rhoh, alphah, Lowerboundary] = ReadEnvParameter(env_file)

    fid           = fopen(env_file);
    casename      = fgetl(fid);
    Layers        = fscanf(fid, '%d', 1);
    Ns            = fscanf(fid, '%d', Layers);
    cpmax         = fscanf(fid, '%f', 1);
    freq          = fscanf(fid, '%f', 1);
    zs            = fscanf(fid, '%f', 1);
    zr            = fscanf(fid, '%f', 1);
    rmax          = fscanf(fid, '%f', 1);
    dr            = fscanf(fid, '%f', 1);
    interface     = fscanf(fid, '%f', Layers);
    dz            = fscanf(fid, '%f', 1);
    tlmin         = fscanf(fid, '%f', 1);
    tlmax         = fscanf(fid, '%f', 1);
    nprofile      = fscanf(fid, '%d', Layers);
    
    dep   = cell(Layers,1);
    c     = cell(Layers,1);   
    rho   = cell(Layers,1);
    alpha = cell(Layers,1);
    
    for i = 1 : Layers
        if (i < Layers && interface(i) > 0.0 && ...
            interface(i) < interface(i+1) && nprofile(i) >= 2)
            Profile     = fscanf(fid, '%f %f', [4, nprofile(i)]);
            dep(i)      = {Profile(1, 1:nprofile(i))};
            c(i)        = {Profile(2, 1:nprofile(i))};
            rho(i)      = {Profile(3, 1:nprofile(i))};
            alpha(i)    = {Profile(4, 1:nprofile(i))};
        elseif(interface(i) > 0.0 && nprofile(i) >= 2)
            Profile     = fscanf(fid, '%f %f', [4, nprofile(i)]);
            dep(i)      = {Profile(1, 1:nprofile(i))};
            c(i)        = {Profile(2, 1:nprofile(i))};
            rho(i)      = {Profile(3, 1:nprofile(i))};
            alpha(i)    = {Profile(4, 1:nprofile(i))};
        else
            error('Error! h must greater than 0 and less than H!');
        end
    end
    
    Lowerboundary = fscanf(fid, '%s', 1);
    
    if (Lowerboundary ~= 'V' && Lowerboundary ~= 'R' && Lowerboundary ~= 'A')
        disp('Error! The lower boundary must be vaccum, rigid or halfspace!');
        error('Please set V, R or A in Lowerboundary!');
    end
    
    if (Lowerboundary == 'A')
        halfspace = fscanf(fid, '%f', 3);
        ch     = halfspace(1);
        rhoh   = halfspace(2);
        alphah = halfspace(3);
    else
        ch     = 0.0;
        rhoh   = 0.0;
        alphah = 0.0;
    end
    
    % Check the input underwater sound profile    
    for i = 1 : Layers
        if (dep{i}(end) ~= interface(i))
            error('Error! input sound profile is unsuitable!');
        end
        if (Ns(i) < 2)
            error('N must greater than 2!');
        end
    end
    
    if ( interface(end) / dz - floor( interface(end) / dz ) ~= 0 )
        error('Error! The input dz unsuitable!');
    end  

    if ((rmax / dr - floor(rmax / dr)) ~= 0)
        error('Please reinput the dr and rmax!');
    end

    if (zs <= 0  || zs >= interface(end) || zr <= 0 || zr >= interface(end))
        error('zs and zr must be greater than 0 and less than H!');
    end
    
    if (tlmin >= tlmax)
        error('tlmin must less than tlmax!');
    end    

    fclose(fid);

end
