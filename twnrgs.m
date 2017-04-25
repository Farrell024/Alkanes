function E = twnrgs(rrr, n, m, rsltn, qq, tp, aa, bb, dd, nbndd)

  rrrrf = rrr;
  
  s = 360/(rsltn);
  
  for a = 1:rsltn
  
    rrr = rttntn(rrr, s/2*(2*a-1), n);
    
    for b = 1:rsltn
    
      rrr = rttntn(rrr, s/2*(2*b-1), m);
      
      E(a, b) = nrgntn(rrr, rrrrf, qq, tp, aa, bb, dd, nbndd);
  
    end
    
  end

end