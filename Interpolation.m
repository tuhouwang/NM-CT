function [c,rho,alpha] = Interpolation(dep,c,rho,alpha,N,s,t)

    x = cos( (0 : N) * pi / N )';
    z = ( (t + s) / (t - s) - x ) * (t - s) / 2.0;

    c     = interp1(dep,  c,  z,'linear');
    rho   = interp1(dep, rho, z,'linear');
    alpha = interp1(dep,alpha,z,'linear');

end