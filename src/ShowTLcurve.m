function ShowTLcurve(r, zr, tl_zr)

    disp('plot the transmission loss curve at zr!');

    figure;
    plot(r, tl_zr, 'b-', 'LineWidth', 1.5);

    title (['Depth=', num2str(zr), 'm']);
    xlabel( 'Range (m)'); 
    ylabel('TL (dB)');
    set   (gca, 'YDir', 'reverse');
    set   (gca, 'FontSize', 16, 'FontName', 'Times New Roman');

end
