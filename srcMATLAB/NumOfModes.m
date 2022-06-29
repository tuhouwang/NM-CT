function [nmodes, kr, eigvector] = NumOfModes(Layers, freq, kr, eigvector, cpmax)

    cp   =  2 * pi * freq ./ real(kr);
    mode = find( cp <= cpmax );
    kr   = kr(mode);
   
    for i = 1 : Layers
        eigvector(i) = {eigvector{i}(:, mode)};
    end
   
    nmodes = length( kr );
    if(nmodes == 0)
        error('Incorrect maximum phase speed input!');
    end
    
end
