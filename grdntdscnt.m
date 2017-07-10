function [rrrm, ee] = grdntdscnt(rrr, rrrrf, qq, tp, aa, bb, dd, nbndd, s)

    dbnds =  (length(rrr) - 2)/3 - 1;
    gg = zeros(1, dbnds);
        
    ee = zeros(1, s);
    
    m = 1;
    
    while m < s

        for n = 1:dbnds

            rrrg = rttntn(rrr, 1, n);
            
            gg(n) = nrgntn(rrrg, rrr, qq, tp, aa, bb, dd, nbndd);

        end
        
        gg = gg/norm(gg);
        
        for n = 1:dbnds
            
            rrr = rttntn(rrr, -gg(n), n);
           
        end
        
        ee(m) = nrgntn(rrr, rrrrf, qq, tp, aa, bb, dd, nbndd);
        
        m = m + 1;
    
    end
    
    rrrm = rrr;

end