function ShowWavenumbers(kr,casename)

    figure;
    disp('plot the modal wavenumbers!');
    plot(real(kr),imag(kr),'r*');title(casename);
    xlabel('Real Wavenumber (1/m)');
    ylabel('Imaginary Wavenumber (1/m)');
    set(gca,'FontSize',16,'FontName','Times New Roman');

end
