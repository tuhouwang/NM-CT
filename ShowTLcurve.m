function ShowTLcurve(r,zr,tl_zr)

    figure;
    disp('plot the transmission loss curve at zr!');
    plot(r,tl_zr,'b-','LineWidth',1.5);
    set(gca,'YDir','reverse');
    xlabel( 'Range (m)'); ylabel('TL (dB)');
    title(['Depth=',num2str(zr),'m']);
    set(gca,'FontSize',16,'FontName','Times New Roman');

end