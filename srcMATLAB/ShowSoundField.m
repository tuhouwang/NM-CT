function ShowSoundField(r, z, tl, tlmin, tlmax, casename, interface)

    disp('plot the transmission loss field!');

    figure;
    pcolor( r, z, tl); hold on;
    plot([0, max(r)], [interface, interface], 'k--', 'Linewidth', 1.5);

    title(casename);
    caxis([tlmin tlmax]); 
    colormap(flipud(jet));
    shading flat; 
    view(0, -90);
    xlabel('Range (m)'); 
    ylabel('Depth (m)');
    colorbar('YDir', 'Reverse');
    set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');

end
