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

[casename, Nw, Nb, cpmax, freq, zs, zr, rmax, dr, interface, bottom, ...
    dz, Lowerboundary, tlmin, tlmax, depw, cw, rhow, alphaw, ...
    depb, cb, rhob, alphab] = ReadEnvParameter('input.txt');

[cw, rhow, alphaw] = Interpolation(depw, cw, rhow, alphaw, Nw,  0,   interface);
[cb, rhob, alphab] = Interpolation(depb, cb, rhob, alphab, Nb, interface, bottom);
%----------------------------------------------
[nr, r, rhozs, kw, kb, w] = Initialization(Nw, Nb, freq, rmax, dr, zs, ...
    rhow, rhob, cw, cb, alphaw, alphab, interface, bottom);

[kr, eigvectorw, eigvectorb] = EigenValueVector(Nw, Nb, interface, ...
    bottom, kw, kb, rhow, rhob, Lowerboundary);

[nmodes, kr, eigvectorw, eigvectorb] = NumOfModes(w, kr, ...
    eigvectorw, eigvectorb, cpmax);

[psi, psizs, z] = GenerateModes(eigvectorw, eigvectorb, nmodes, dz, ...
    zs, rhow, rhob, interface, bottom);

[tl, tl_zr] = SynthesizeSoundField(r, z, kr, rhozs, psizs, psi, zr);

%---------------Show the results-------------------
ShowWavenumbers(kr, casename);
ShowTLcurve(r, zr, tl_zr);
ShowSoundField(r, z, tl, tlmin, tlmax, casename, interface);
% ShowMode(psi, z);
% SaveSoundField('tl.bin', tlmin, tlmax, r, z, tl);
toc;
