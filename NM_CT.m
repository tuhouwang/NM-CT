%-------------------------------------------------------------------
% Compute the 2 layers of range-independent modal acoustic field    |
% using the Chebyshev-Tau spectral method based on the normal modes.|
% This code was written by Tu Houwang at 06/13/2020 from NUDT.      |
% The water and sediment can be specified in different 'N' by user. |
% ------------------------------------------------------------------
edit 'input.txt';
clear;
close all;
clc;
tic;

[casename,Nw,Nb,cpmax,dr,zs,zr,rmax,freq,interface,bottom,dz,...
    Lowerboundary,tlmin,tlmax,rhow,rhob,alphaw,alphab,depw,cw,...
    depb,cb]= Read_input('input.txt');

cw     = Interpolation(depw,cw,Nw,0,interface);
cb     = Interpolation(depb,cb,Nb,interface,bottom);
rhow   = Interpolation(depw,rhow,Nw,0,interface);
rhob   = Interpolation(depb,rhob,Nb,interface,bottom);
alphaw = Interpolation(depw,alphaw,Nw,0,interface);
alphab = Interpolation(depb,alphab,Nb,interface,bottom);

%----------------------------------------------

[nr,r,rho_zs,kw,kb]=Initialization(Nw,Nb,freq,rmax,dr,zs,...
    rhow,rhob,cw,cb,alphaw,alphab,interface,bottom);

[kr,eigvector_w,eigvector_b]=Eigen_value_vector(Nw,Nb,interface,...
    bottom,kw,kb,rhow,rhob,Lowerboundary);

[n_modes,kr,eigvector_w,eigvector_b]=Modes_number(freq,kr,...
    eigvector_w,eigvector_b,cpmax);

[psi,z]  = Generate_modes(eigvector_w,eigvector_b,n_modes,dz,...
    rhow,rhob,interface,bottom);

[tl,tl_zr] = Calculate_tl_field(n_modes,nr,r,kr,zs,rho_zs,dz,psi,zr);

%---------------print the results-------------------
Print_wavenumbers(kr,casename);
% Print_mode(psi,z);
Print_tl_zr(r,zr,tl_zr);
Print_tl_field(r,z,tl,tlmin,tlmax,casename,interface);
% Write_to_file('tl.bin',tlmin,tlmax,r,z,tl);
toc;
%--------------------------------------------------------------------------
function [casename,Nw,Nb,cpmax,dr,zs,zr,rmax,freq,interface,...
    bottom,dz,Lowerboundary,tlmin,tlmax,rhow,rhob,...
    alphaw,alphab,depw,cw,depb,cb] = Read_input(env_file)

    fid          = fopen(env_file);
    casename     = fgetl(fid);
    Nw           = fscanf(fid,'%d',1);
    Nb           = fscanf(fid,'%d',1);
    cpmax        = fscanf(fid,'%f',1);
    freq         = fscanf(fid,'%f',1);
    zs           = fscanf(fid,'%f',1);
    zr           = fscanf(fid,'%f',1);
    rmax         = fscanf(fid,'%f',1);
    dr           = fscanf(fid,'%f',1);
    interface    = fscanf(fid,'%f',1);
    bottom       = fscanf(fid,'%f',1);
    dz           = fscanf(fid,'%f',1);
    Lowerboundary= fscanf(fid,'%d',1);
    tlmin        = fscanf(fid,'%f',1);
    tlmax        = fscanf(fid,'%f',1);
    nw           = fscanf(fid,'%d',1);
    nb           = fscanf(fid,'%d',1);
    if(interface>0.0 && interface<bottom && Nw>2 && Nb>2)
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
    if(depw(1) ~=0.0 || depw(nw)~=interface ||...
            depb(1)~=interface || depb(nb)~=bottom)
        error('Error! input sound profile is unsuitable!');
    end

    if((interface/dz-floor(interface/dz))~=0 ||...
            (bottom/dz-floor(bottom/dz))~=0)
        error('Error! The input d unsuitable!');
    end

    if((rmax/dr-floor(rmax/dr))~=0)
        error('Please reinput the dr and rmax!');
    end

    if(Lowerboundary~=0 && Lowerboundary~=1)
        disp('Error! The lower boundary must be rigid or soft!');
        error('Please input 0 or 1 in Lowerboundary!');
    end

    if(tlmin >= tlmax)
        error('tlmin must less than tlmax!');
    end

    fclose(fid);

end

function c = Interpolation(depc,c,N,s,t)

    x  = cos((0:N)*pi/N)';
    z  = ((t+s)/(t-s) - x).*(t-s)/2.0;
    c  = interp1(depc,c,z,'linear','extrap');

end

function [nr,r,rho_zs,kw,kb] = Initialization(Nw,Nb,freq,rmax,...
    dr,zs,rhow,rhob,cw,cb,alphaw,alphab,interface,bottom)
    
    w       = 2 * pi * freq;
    nr      = rmax/dr;
    r       = dr:dr:rmax;

    x1   = cos((0 : Nw) * pi / Nw)';
    x2   = cos((0 : Nb) * pi / Nb)';
    if(zs <= interface)
        z1    = (1.0 - x1) * interface/2;
        rho_zs=interp1(z1,rhow,zs,'linear');
    else
        z1    = (1.0 - x2) * (bottom-interface)/2 +interface;
        rho_zs=interp1(z1,rhob,zs,'linear');
    end

    kw = w./cw.*(1.0+1i.*alphaw/(40.0*pi*log10(exp(1.0))));
    kb = w./cb.*(1.0+1i.*alphab/(40.0*pi*log10(exp(1.0))));
    
end

function D  = DeriveMatrix(n)

    D = zeros(n, n);
    for k = 1 : n
        j = k+1 : 2 : n;
        D(k, j) = 2 * j - 2;
    end
    D(1,:) = D(1,:) / 2;
        
end

function T  = ChebPolynomialValue(n,z)

    m = length(z);

    T = zeros(m, n);
    for k = 0 : n-1
        T(:, k+1) = cos( k * acos(z) );
    end

end

function fx = invCheb(fk,z)

    n = size(fk, 1);
    T = ChebPolynomialValue(n, z);
    fx = T * fk;
    
end

function C  = Convolution(v)

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

function fk = Cheb(fx,z)

    n = length(fx);
    T = ChebPolynomialValue(n, z);

    fx(1) = fx(1) / 2;
    fx(n) = fx(n) / 2;

    fk = T.' * fx .* (2 / (n-1));

    fk(1) = fk(1) / 2;
    fk(n) = fk(n) / 2;
    
end

function [kr,eigvector_w,eigvector_b] = Eigen_value_vector(Nw,Nb,...
    interface,bottom,kw,kb,rhow,rhob,Lowerboundary)

    D1  = DeriveMatrix(Nw+1);
    D2  = DeriveMatrix(Nb+1);
    
    x1  = cos( (0 : Nw) * pi / Nw )';
    x2  = cos( (0 : Nb) * pi / Nb )';
   
    A =4.0/interface^2.*Convolution(Cheb(rhow,x1))*...
        D1*Convolution(Cheb(1./rhow,x1))*D1+...
        Convolution( Cheb(kw.^2, x1 ) );
    B =4.0/(bottom-interface)^2.*Convolution(Cheb(rhob,x2))*...
        D2*Convolution(Cheb(1./rhob,x2))*D2+...
        Convolution( Cheb(kb.^2, x2 ) );
    
    ALL=zeros(Nw+Nb+2,Nw+Nb+2);
    ALL(1:Nw-1,1:Nw-1)=A(1:Nw-1,1:Nw-1);
    ALL(1:Nw-1,Nw+Nb-1:Nw+Nb)=A(1:Nw-1,Nw:Nw+1);
    ALL(Nw:Nw+Nb-2,Nw:Nw+Nb-2)=B(1:Nb-1,1:Nb-1);
    ALL(Nw:Nw+Nb-2,Nw+Nb+1:Nw+Nb+2)=B(1:Nb-1,Nb:Nb+1);
    %upper boundary
    ALL(Nw+Nb-1,1:Nw-1)        =1.0;
    ALL(Nw+Nb-1,Nw+Nb-1:Nw+Nb) =1.0;
    %lower boundary
    bottom_boundary=(-1.0).^(0:Nb);
    if(Lowerboundary==1)
        bottom_boundary=bottom_boundary*D2;
    end
    ALL(Nw+Nb+2,Nw:Nw+Nb-2)     = bottom_boundary(1:Nb-1);
    ALL(Nw+Nb+2,Nw+Nb+1:Nw+Nb+2)= bottom_boundary(Nb:Nb+1);
    %interface boundary 1
    ALL(Nw+Nb,1:Nw-1)         = (-1.0).^(0:Nw-2);
    ALL(Nw+Nb,Nw:Nw+Nb-2)     = -1.0;
    ALL(Nw+Nb,Nw+Nb-1:Nw+Nb)  = (-1.0).^(Nw-1:Nw);
    ALL(Nw+Nb,Nw+Nb+1:Nw+Nb+2)= -1.0;
    %interface boundary 2
    Pu=1/ rhow(Nw+1)/interface* ((-1.0).^(0:Nw)) * D1;
    Pd=-1/rhob(1) / (bottom-interface) * ones(1,Nb+1)* D2;

    ALL(Nw+Nb+1,1:Nw-1)         = Pu(1:Nw-1);
    ALL(Nw+Nb+1,Nw:Nw+Nb-2)     = Pd(1:Nb-1);
    ALL(Nw+Nb+1,Nw+Nb-1:Nw+Nb)  = Pu(Nw:Nw+1);
    ALL(Nw+Nb+1,Nw+Nb+1:Nw+Nb+2)= Pd(Nb:Nb+1);
    %blocking
    L11=ALL(1:Nw+Nb-2,1:Nw+Nb-2);
    L12=ALL(1:Nw+Nb-2,Nw+Nb-1:Nw+Nb+2);
    L21=ALL(Nw+Nb-1:Nw+Nb+2,1:Nw+Nb-2);
    L22=ALL(Nw+Nb-1:Nw+Nb+2,Nw+Nb-1:Nw+Nb+2);

    L=L11-L12*(L22\L21);
    [v,k2]=eig(L);

    v2=-(L22\L21)*v;

    eigvector_w=[v(1:Nw-1,:);v2(1:2,:)];
    eigvector_b=[v(Nw:Nw+Nb-2,:);v2(3:4,:)];

    k2=sqrt(diag(k2));
    [~,ind]=sort(real(k2),'descend');
    kr=k2(ind);

    eigvector_w=eigvector_w(:,ind);
    eigvector_b=eigvector_b(:,ind);

end

function [n_modes,kr,eigvector_w,eigvector_b]= Modes_number...
    (freq,kr,eigvector_w,eigvector_b,cpmax)

    cp =2*pi*freq./real(kr);
    n_modes  =0;
    for i=1:length(kr)
        if(cp(i)<=cpmax)
            n_modes=i;
        end
    end

    if(n_modes==0)
        error('Incorrect maximum phase speed input!');
    end

    kr=kr(1:n_modes);
    eigvector_w=eigvector_w(:,1:n_modes);
    eigvector_b=eigvector_b(:,1:n_modes);

end

function a = Normalization_2layers(eigvector_w,eigvector_b,n_modes,...
    rhow,rhob,interface,bottom)

    Nw=size(eigvector_w,1)-1;
    Nb=size(eigvector_b,1)-1;
    x1=cos((0:Nw)*pi/Nw)';
    x2=cos((0:Nb)*pi/Nb)';
    Rw=Convolution(Cheb(1./rhow,x1));
    Rb=Convolution(Cheb(1./rhob,x2));

    a=zeros(n_modes,1);

    P1      = zeros(Nw+1,1);
    k       = 0:2:Nw;
    P1(k+1) = -2./(k.^2-1);
    P2      = zeros(Nb+1,1);
    k       = 0:2:Nb;
    P2(k+1) = -2./(k.^2-1);

    for j=1:n_modes

        f1 = Convolution(eigvector_w(:,j))*eigvector_w(:,j);
        f2 = Convolution(eigvector_b(:,j))*eigvector_b(:,j);
        f1 = Rw*f1;        
        f2 = Rb*f2;

        a(j)=sqrt(dot(P1,f1)*interface/2+dot(P2,f2)*(bottom-interface)/2);
    end

end

function [psi,z] = Generate_modes(eigvector_w,eigvector_b,n_modes,dz,...
    rhow,rhob,interface,bottom)

    zt1  = 0:dz:interface;
    zt2  = interface:dz:bottom;

    xt1  = -2/interface*zt1+1;
    xt2  = -2/(bottom-interface)*zt2+(bottom+interface)/(bottom-interface);

    psi1 = invCheb(eigvector_w,xt1);
    psi2 = invCheb(eigvector_b,xt2);
    psi  = [psi1(1:length(xt1)-1,:); psi2];

    a    = Normalization_2layers(eigvector_w,eigvector_b,...
           n_modes,rhow,rhob,interface,bottom);

    for j=1:n_modes
        psi (:,j)=psi(:,j)./a(j);
    end

    z=0:dz:bottom;

end

function [tl,tl_zr] = Calculate_tl_field(n_modes,nr,r,kr,zs,rho_zs,dz,psi,zr)

    bessel=zeros(n_modes,nr);
    for im=1:n_modes
        for ir=1:nr
            bessel(im,ir)=besselh(0,1,r(ir)*kr(im));
        end
    end
    s=round(zs/dz)+1;

    psi = psi*diag(psi(s,:));
    p   = psi*bessel*1i*pi/rho_zs;

    tl   = -20 * log10(abs(p));
    tl_zr=tl(round(zr/dz)+1,:);

end

function Print_wavenumbers(kr,casename)

    disp('plot the modal wavenumbers!');
    figure;
    plot(real(kr),imag(kr),'r*');
    xlabel('Real Wavenumber (1/m)');
    ylabel('Imaginary Wavenumber (1/m)');
    title(casename);
    set(gca,'FontSize',16,'FontName','Times New Roman');

end

function Print_mode(psi,z)

    figure;
    mode_num = input('What mode number do you want to plot?:');
    plot(imag(psi(:,mode_num)),z,'k--','LineWidth',1);hold on;
    plot(real(psi(:,mode_num)),z,'r-','LineWidth',0.5);
    set(gca,'YDir','reverse');ylabel( 'Depth (m)');
    set(gca,'FontSize',16,'FontName','Times New Roman');

end

function Print_tl_zr(r,zr,tl_zr)

    figure;
    disp('plot the transmission loss curve at zr!');
    plot(r,tl_zr,'b-','LineWidth',1.5);
    set(gca,'YDir','reverse');
    xlabel( 'Range (m)'); ylabel('TL (dB)');
    title(['Depth=',num2str(zr),'m']);
    set(gca,'FontSize',16,'FontName','Times New Roman');

end

function Print_tl_field(r,z,tl,tlmin,tlmax,casename,interface)

    figure;
    disp('plot the transmission loss field!');
    pcolor( r, z, tl); hold on;
    plot(r,interface*ones(length(r)),'k--','Linewidth',1.5);
    caxis( [tlmin tlmax] ); colormap( flipud(jet) );
    shading flat; colorbar; view( 0, -90 );
    xlabel( 'Range (m)'); ylabel( 'Depth (m)');
    colorbar( 'YDir', 'Reverse' );title(casename);
    set(gca,'FontSize',16,'FontName','Times New Roman');

end

function Write_to_file(filename,tlmin,tlmax,r,z,tl)

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