% Ocean acoustic normal modes.

% Copyright (C) 2020 Houwang Tu
% -------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify it |
% under the terms of the GNU General Public License as published by the   |
% Free Software Foundation, either version 3 of the License, or (at your  |
% option) any later version.                                              |
%                                                                         |
% This code is distributed in the hope that it will be useful, but without|
% any warranty; without even the implied warranty of merchantability or   |
% fitness for a particular purpose. See the GNU General Public License for|
% more details.                                                           |
%                                                                         |
% You should have received a copy of the GNU General Public License along |
% with this program. If not, see <http://www.gnu.org/licenses/>.          |
%                                                                         |
% Originally developed as part of the author's article (H.Tu, Y.Wang, Q.  |
% Lan et al., A Chebyshev-Tau spectral method for normal modes of         |
% underwater sound propagation with a layered marine environment, Journal |
% of Sound and Vibration, https://doi.org/10.1016/j.jsv.2020.115784) under|
% the supervision of Prof. Yongxian Wang, National University of Defense  |
% Technology, China.                                                      |
%																		  |
% This Matlab/Scilab style code computes the layered and range-independent|
% modal acoustic field using the Chebyshev-Tau spectral method based on   |
% the normal modes.                                                       |
% -------------------------------------------------------------------------
% edit 'input.txt';
clear;
close all;
clc;
tic;

[casename,Nw,Nb,cpmax,dr,zs,zr,rmax,freq,interface,bottom,dz,...
    Lowerboundary,tlmin,tlmax,depw,cw,rhow,alphaw,...
    depb,cb,rhob,alphab] = ReadEnvParameter('input.txt');

[cw,rhow,alphaw] = Interpolation(depw,cw,rhow,alphaw,Nw,  0,   interface);
[cb,rhob,alphab] = Interpolation(depb,cb,rhob,alphab,Nb,interface,bottom);
%----------------------------------------------
[nr,r,rhozs,kw,kb,w] = Initialization(Nw,Nb,freq,rmax,dr,zs,...
    rhow,rhob,cw,cb,alphaw,alphab,interface,bottom);

[kr,eigvectorw,eigvectorb] = EigenValueVector(Nw,Nb,interface,...
    bottom,kw,kb,rhow,rhob,Lowerboundary);

[nmodes,kr,eigvectorw,eigvectorb] = NumofModes(w,kr,...
    eigvectorw,eigvectorb,cpmax);

[psi,psizs,z] = GenerateModes(eigvectorw,eigvectorb,nmodes,dz,...
    zs,rhow,rhob,interface,bottom);

[tl,tl_zr] = SynthesizeSoundField(r,z,kr,rhozs,psizs,psi,zr);

%---------------Show the results-------------------
ShowWavenumbers(kr,casename);
% ShowMode(psi,z);
ShowTLcurve(r,zr,tl_zr);
ShowSoundField(r,z,tl,tlmin,tlmax,casename,interface);
% SaveSoundField('tl.bin',tlmin,tlmax,r,z,tl);
toc;
%--------------------------------------------------------------------------
function [casename,Nw,Nb,cpmax,dr,zs,zr,rmax,freq,interface,...
    bottom,dz,Lowerboundary,tlmin,tlmax,depw,cw,rhow,alphaw,...
    depb,cb,rhob,alphab] = ReadEnvParameter(env_file)

fid           = fopen(env_file);
casename      = fgetl(fid);
Nw            = fscanf(fid,'%d',1);
Nb            = fscanf(fid,'%d',1);
cpmax         = fscanf(fid,'%f',1);
freq          = fscanf(fid,'%f',1);
zs            = fscanf(fid,'%f',1);
zr            = fscanf(fid,'%f',1);
rmax          = fscanf(fid,'%f',1);
dr            = fscanf(fid,'%f',1);
interface     = fscanf(fid,'%f',1);
bottom        = fscanf(fid,'%f',1);
dz            = fscanf(fid,'%f',1);
Lowerboundary = fscanf(fid,'%d',1);
tlmin         = fscanf(fid,'%f',1);
tlmax         = fscanf(fid,'%f',1);
nw            = fscanf(fid,'%d',1);
nb            = fscanf(fid,'%d',1);
if(interface > 0.0 && interface < bottom && Nw > 2 && Nb > 2)
    W_Profile    = fscanf(fid,'%f %f',[4,nw]);
    depw(1:nw)   = W_Profile(1,1:nw);
    cw(1:nw)     = W_Profile(2,1:nw);
    rhow(1:nw)   = W_Profile(3,1:nw);
    alphaw(1:nw) = W_Profile(4,1:nw);
    B_Profile    = fscanf(fid,'%f %f',[4,nb]);
    depb         = B_Profile(1,1:nb);
    cb(1:nb)     = B_Profile(2,1:nb);
    rhob(1:nb)   = B_Profile(3,1:nb);
    alphab(1:nb) = B_Profile(4,1:nb);
else
    error('Error! h must less than H and greater than 0!');
end
% Check the input underwater sound profile
if(depw(1) ~= 0.0 || depw(nw) ~= interface ||...
        depb(1) ~= interface || depb(nb) ~= bottom)
    error('Error! input sound profile is unsuitable!');
end

if((interface / dz - floor(interface / dz)) ~=0 ||...
        (bottom / dz - floor(bottom / dz)) ~=0)
    error('Error! The input dz unsuitable!');
end

if((rmax / dr - floor(rmax / dr)) ~= 0)
    error('Please reinput the dr and rmax!');
end

if(Lowerboundary ~= 0 && Lowerboundary ~= 1)
    disp('Error! The lower boundary must be rigid or soft!');
    error('Please input 0 or 1 in Lowerboundary!');
end

if(tlmin >= tlmax)
    error('tlmin must less than tlmax!');
end

fclose(fid);

end

function [c,rho,alpha] = Interpolation(dep,c,rho,alpha,N,s,t)

x = cos( (0 : N) * pi / N )';
z = ( (t + s) / (t - s) - x ) * (t - s) / 2.0;

c     = interp1(dep,c,z,'linear');
rho   = interp1(dep,rho,z,'linear');
alpha = interp1(dep,alpha,z,'linear');

end

function [nr,r,rhozs,kw,kb,w] = Initialization(Nw,Nb,freq,rmax,...
    dr,zs,rhow,rhob,cw,cb,alphaw,alphab,interface,bottom)

w       = 2 * pi * freq;
r       = dr : dr : rmax;
nr      = length(r);

if(zs <= interface)
    x1     = cos( (0 : Nw) * pi / Nw)';
    z1     = (1.0 - x1) * interface / 2;
    rhozs  = interp1(z1,rhow,zs,'linear');
else
    x2     = cos( (0 : Nb) * pi / Nb)';
    z2     = (1.0 - x2) * (bottom - interface) / 2 + interface;
    rhozs  = interp1(z2,rhob,zs,'linear');
end

kw = w ./ cw .* (1.0+1i*alphaw/(40.0*pi*log10(exp(1.0))));
kb = w ./ cb .* (1.0+1i*alphab/(40.0*pi*log10(exp(1.0))));

end

function D  = DerivationMatrix(n)

D = zeros(n, n);
for k = 1 : n
    j = k+1 : 2 : n;
    D(k, j) = 2 * j - 2;
end
D(1, :) = D(1, :) / 2;

end

function C  = ConvolutionMatrix(v)

n = length(v);
C = zeros(n, n);
for k = 1 : n
    for i = 1 : n
        for j = 1 : n
            if ((i-1 + j-1) == (k-1))
                C(k,i) = C(k,i) + v(j) / 2;
            end
            if (abs(i - j) == (k-1))
                C(k,i) = C(k,i) + v(j) / 2;
            end
        end
    end
end

end

function T  = ChebPolyValue(n,z)

m = length(z);

T = zeros(m, n);
for k = 0 : n-1
    T(:, k+1) = cos( k * acos(z) );
end

end

function fx = InvChebTrans(fk,z)

n = size(fk, 1);
T = ChebPolyValue(n, z);
fx = T * fk;

end

function fk = ChebTrans(fx,z)

n = length(fx);
T = ChebPolyValue(n, z);

fx(1) = fx(1) / 2;
fx(n) = fx(n) / 2;

fk = T' * fx * (2 / (n-1));

fk(1) = fk(1) / 2;
fk(n) = fk(n) / 2;

end

function [kr,eigvectorw,eigvectorb] = EigenValueVector(Nw,Nb,...
    interface,bottom,kw,kb,rhow,rhob,Lowerboundary)

D1  = DerivationMatrix(Nw+1);
D2  = DerivationMatrix(Nb+1);
x1  = cos( (0 : Nw) * pi / Nw )';
x2  = cos( (0 : Nb) * pi / Nb )';

A =4.0 / interface^2 * ConvolutionMatrix(ChebTrans(rhow, x1)) *...
    D1 * ConvolutionMatrix(ChebTrans(1./rhow, x1)) * D1 +...
    ConvolutionMatrix( ChebTrans(kw.^2, x1 ) );
B =4.0 / (bottom-interface)^2 * ConvolutionMatrix(ChebTrans(rhob, x2))...
    *D2 * ConvolutionMatrix(ChebTrans(1./rhob, x2)) * D2 +...
    ConvolutionMatrix( ChebTrans(kb.^2, x2 ) );

U = zeros(Nw+Nb+2,Nw+Nb+2);
U(1:Nw-1,1:Nw-1)              = A(1:Nw-1,1:Nw-1);
U(1:Nw-1,Nw+Nb-1:Nw+Nb)       = A(1:Nw-1,Nw:Nw+1);
U(Nw:Nw+Nb-2,Nw:Nw+Nb-2)      = B(1:Nb-1,1:Nb-1);
U(Nw:Nw+Nb-2,Nw+Nb+1:Nw+Nb+2) = B(1:Nb-1,Nb:Nb+1);
%upper boundary
U(Nw+Nb-1,1:Nw-1)        = 1.0;
U(Nw+Nb-1,Nw+Nb-1:Nw+Nb) = 1.0;
%lower boundary
bottom_boundary = (-1.0) .^ (0:Nb);
if(Lowerboundary == 1)
    bottom_boundary = bottom_boundary * D2;
end
U(Nw+Nb+2,Nw:Nw+Nb-2)      = bottom_boundary(1 :Nb-1);
U(Nw+Nb+2,Nw+Nb+1:Nw+Nb+2) = bottom_boundary(Nb:Nb+1);
%first interface boundary
U(Nw+Nb,1:Nw-1)          = (-1.0).^(0 : Nw-2);
U(Nw+Nb,Nw:Nw+Nb-2)      = -1.0;
U(Nw+Nb,Nw+Nb-1:Nw+Nb)   = (-1.0).^(Nw-1 : Nw);
U(Nw+Nb,Nw+Nb+1:Nw+Nb+2) = -1.0;
%second interface boundary 2
Pu = 1 / rhow(Nw+1) / interface * ((-1.0) .^ (0 : Nw)) * D1;
Pd = -1/ rhob(1) / (bottom - interface) * ones(1,Nb+1) * D2;

U(Nw+Nb+1,1:Nw-1)          = Pu(1  : Nw-1);
U(Nw+Nb+1,Nw:Nw+Nb-2)      = Pd(1  : Nb-1);
U(Nw+Nb+1,Nw+Nb-1:Nw+Nb)   = Pu(Nw : Nw+1);
U(Nw+Nb+1,Nw+Nb+1:Nw+Nb+2) = Pd(Nb : Nb+1);
%blocking
L11 = U(1 : Nw+Nb-2,            1 : Nw+Nb-2);
L12 = U(1 : Nw+Nb-2,      Nw+Nb-1 : Nw+Nb+2);
L21 = U(Nw+Nb-1 : Nw+Nb+2,      1 : Nw+Nb-2);
L22 = U(Nw+Nb-1 : Nw+Nb+2,Nw+Nb-1 : Nw+Nb+2);

L = L11 - L12 * (L22 \ L21);
[v,k2] = eig(L);

v2 = - (L22 \ L21) * v;

eigvectorw = [v(1  : Nw-1, :)    ; v2(1:2, :)];
eigvectorb = [v(Nw : Nw+Nb-2, :) ; v2(3:4, :)];

k2      = sqrt(diag(k2));
[~,ind] = sort(real(k2),'descend');
kr      = k2(ind);

eigvectorw = eigvectorw(:,ind);
eigvectorb = eigvectorb(:,ind);

end

function [nmodes,kr,eigvectorw,eigvectorb] = NumofModes...
    (w,kr,eigvectorw,eigvectorb,cpmax)

cp = w ./ real(kr);
nmodes  = 0;
for i = 1 : length(kr) - 1
    if(cp(i) <= cpmax)
        nmodes = i;
    end
end

if(nmodes == 0)
    error('Incorrect maximum phase speed input!');
end

kr = kr(1 : nmodes);
eigvectorw = eigvectorw(:, 1:nmodes);
eigvectorb = eigvectorb(:, 1:nmodes);

end

function a = Normalization(eigvectorw,eigvectorb,nmodes,...
    rhow,rhob,interface,bottom)

Nw = size(eigvectorw, 1) - 1;
Nb = size(eigvectorb, 1) - 1;
x1 = cos( (0 : Nw) * pi / Nw)';
x2 = cos( (0 : Nb) * pi / Nb)';
Rw = ConvolutionMatrix(ChebTrans(1./rhow, x1));
Rb = ConvolutionMatrix(ChebTrans(1./rhob, x2));

a = zeros(nmodes, 1);

P1      = zeros(1, Nw+1);
P2      = zeros(1, Nb+1);

k       = 0 : 2 : Nw;
P1(k+1) = -2 ./ (k.^2 - 1);
k       = 0 : 2 : Nb;
P2(k+1) = -2 ./ (k.^2 - 1);

for j = 1 : nmodes
    
    f1 = ConvolutionMatrix(eigvectorw(:,j)) * eigvectorw(:,j);
    f2 = ConvolutionMatrix(eigvectorb(:,j)) * eigvectorb(:,j);
    f1 = Rw * f1;
    f2 = Rb * f2;
    
    a(j) = sqrt( P1 * f1 * interface / 2 + ...
        P2 * f2 * (bottom - interface) / 2 );
end

end

function [psi,psizs,z] = GenerateModes(eigvectorw,eigvectorb,nmodes,dz,...
    zs,rhow,rhob,interface,bottom)

zt1  = 0 : dz : interface;
zt2  = interface : dz : bottom;
z    = 0 : dz : bottom;

xt1  = -2 / interface * zt1 + 1;
xt2  = -2 / (bottom - interface) * zt2 +...
    (bottom + interface) / (bottom - interface);

psi1 = InvChebTrans(eigvectorw, xt1);
psi2 = InvChebTrans(eigvectorb, xt2);
psi  = [psi1(1 : length(xt1) - 1, :); psi2];

a    = Normalization(eigvectorw,eigvectorb,...
    nmodes,rhow,rhob,interface,bottom);

psi  = psi * diag(1 ./ a);

psizs = interp1(z,psi,zs,'linear');
end

function [tl,tl_zr] = SynthesizeSoundField(r,z,kr,rhozs,psizs,psi,zr)

bessel = besselh(0,1,kr * r);
p      = psi * diag( psizs ) * bessel * 1i * pi / rhozs;
tl     = -20 * log10( abs(p) );
tl_zr  = interp1(z,tl,zr,'linear');

end

function ShowWavenumbers(kr,casename)

disp('plot the modal wavenumbers!');
figure;
plot(real(kr),imag(kr),'r*');title(casename);
xlabel('Real Wavenumber (1/m)');
ylabel('Imaginary Wavenumber (1/m)');
set(gca,'FontSize',16,'FontName','Times New Roman');

end

function ShowMode(psi,z)

figure;
mode_num = input('What mode number do you want to plot?:');
plot(imag(psi(:,mode_num)),z,'k--','LineWidth',1) ;hold on;
plot(real(psi(:,mode_num)),z,'r-','LineWidth',0.5);
set(gca,'YDir','reverse');ylabel( 'Depth (m)');
set(gca,'FontSize',16,'FontName','Times New Roman');

end

function ShowTLcurve(r,zr,tl_zr)

figure;
disp('plot the transmission loss curve at zr!');
plot(r,tl_zr,'b-','LineWidth',1.5);
set(gca,'YDir','reverse');
xlabel( 'Range (m)'); ylabel('TL (dB)');
title(['Depth=',num2str(zr),'m']);
set(gca,'FontSize',16,'FontName','Times New Roman');

end

function ShowSoundField(r,z,tl,tlmin,tlmax,casename,interface)

figure;
disp('plot the transmission loss field!');
pcolor( r, z, tl); hold on;
plot(r,interface*ones(length(r)),'k--','Linewidth',1.5);
caxis( [tlmin tlmax] ); colormap( flipud(jet) );
shading flat; colorbar; view( 0, -90 );
xlabel( 'Range (m)')  ; ylabel( 'Depth (m)');
colorbar( 'YDir', 'Reverse' );title(casename);
set(gca,'FontSize',16,'FontName','Times New Roman');

end

function SaveSoundField(filename,tlmin,tlmax,r,z,tl)

tl_fid = fopen(filename,'w');
fwrite(tl_fid,length(z),'int32');
fwrite(tl_fid,length(r),'int32');
fwrite(tl_fid,tlmin,'double');
fwrite(tl_fid,tlmax,'double');
fwrite(tl_fid,z, 'double');
fwrite(tl_fid,r, 'double');
fwrite(tl_fid,tl,'double');
fclose(tl_fid);
end