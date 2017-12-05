function [rrr, dhd, ee] = grdntdscnt(rrr, rrrrf, qq, tp, aa, bb, dd, nbndd, s)

    dbnds =  (length(rrr) - 2)/3 - 3;
    gg = zeros(2, dbnds);
    ee = zeros(1, s);
    dhd = zeros(s, (length(rrr)-2)/3 -3);
    dhd(1,:) = plthlndhdrl(rrr);
    
    m = 1;
    
    while m <= s

        gg(1, :) = gg(2, :);
        
        for n = 2:dbnds+1

            rrrg = rttntn(rrr, 1/60, n);
            
            gg(2, n) = (nrgntnw(rrrg, qq, tp, aa, bb, dd, nbndd)-nrgntnw(rrr, qq, tp, aa, bb, dd, nbndd))*60; 

        end
  
        gg(2,:)
        
        gmm = 1/12/norm(gg(2,:));
        
        for n = 2:dbnds+1
            
            rrr = rttntn(rrr, -gmm.*gg(2,n), n);
           
        end
               
        dhd(m, :) = plthlndhdrl(rrr);
        
        ee(m) = nrgntnw(rrr, qq, tp, aa, bb, dd, nbndd) - nrgntnw(rrrrf, qq, tp, aa, bb, dd, nbndd);
        
        m = m + 1;
    
    end
 
end