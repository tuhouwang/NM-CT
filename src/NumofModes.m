function [nmodes,kr,v1,v2] = NumofModes(w,kr,v1,v2,cpmax)

    cp     = w ./ real(kr);
    nmodes = 0;
    for i=1 : length(kr)
        if(cp(i) <= cpmax)
            nmodes = i;
        end
    end

    if(nmodes == 0)
        error('Incorrect maximum phase speed input!');
    end
    kr = kr(  1 : nmodes);
    v1 = v1(:,1 : nmodes);
    v2 = v2(:,1 : nmodes);
    
end
