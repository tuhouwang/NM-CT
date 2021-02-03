function ShowMode(psi, z)

    mode_num = input('What mode number do you want to plot?:');

    figure;
    plot(imag(psi(:, mode_num)), z, 'k--', 'LineWidth', 1) ;
    hold on;
    plot(real(psi(:, mode_num)), z, 'r-', 'LineWidth', 0.5);

    ylabel( 'Depth (m)');
    set(gca, 'YDir', 'reverse');
    set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');

end
