function [tl,tl_zr] = SynthesizeSoundField(r,z,kr,rhozs,psizs,psi,zr)

    bessel = besselh(0,1,kr * r);
    p      = psi * diag( psizs ) * bessel * 1i * pi / rhozs;
    tl     = -20 * log10( abs(p) );
    tl_zr  = interp1(z,tl,zr,'linear');

end
