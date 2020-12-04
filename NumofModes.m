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