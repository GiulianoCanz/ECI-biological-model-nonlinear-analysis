function dx = ECI_equations(t,x,g,k,c_c,l,n) 

    dx = [g.c0*hill_function(x(1),c_c.CC0,n.cc,l.CCp) - k.c*x(1)*hill_function(x(2),c_c.DC0,n.dc,l.DCp)*hill_function(x(3),c_c.KC0,n.kc,l.KCp);
        g.d0*hill_function(x(2),c_c.DD0,n.dd,l.DDp)*hill_function(x(3),c_c.KD0,n.kd,l.KDp)*hill_function(x(1),c_c.CD10,n.cd1,l.CD1m)*hill_function(x(1),c_c.CD20,n.cd2,l.CD2p) - k.d * x(2);
        g.k0*x(2)*hill_function(x(1),c_c.CK10,n.ck1,l.CK1p)*hill_function(x(1),c_c.CK20,n.ck2,l.CK2m) - k.k*x(3)*hill_function(x(1),c_c.CK30,n.ck3,l.CK3p)
        ];
    
    % il vettore di stato è x = [C,D,K]
    % x(1) = C;
    % x(2) = D;
    % x(3) = K;
    end
    
    