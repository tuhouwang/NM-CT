function ShowSoundField(r, z, tl, tlmin, tlmax, casename, interface)

    disp('plot the transmission loss field!');
    r = r ./ 1000;
    figure;     
    pcolor(r, z, tl); 
    shading flat; 
    view(0, -90);
    title(casename);
    caxis([tlmin tlmax]); 
    xlabel('Range (km)'); ylabel('Depth (m)');
    colormap(flipud(jet)); colorbar('YDir', 'Reverse');
    set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
    
    hold on;
    for i = 1 : length(interface) - 1
        plot([0, max(r)], [interface(i), interface(i)], 'k--', 'Linewidth', 1.5);
    end

end
