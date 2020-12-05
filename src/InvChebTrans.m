function fx = InvChebTrans(fk, x)

    N  = size(fk, 1) - 1;
    T  = ChebPoly(N, x);   
    fx = T * fk;

end
