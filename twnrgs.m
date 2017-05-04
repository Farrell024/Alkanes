function E = twnrgs(rrr, n, m, rsltn, qq, tp, aa, bb, dd, nbndd)

  rrrrf = rrr;
  
  s = 360/(rsltn+1);
  
  E = zeros(rsltn);
  
  for a = 1:rsltn
  
    rrr = rttntn(rrrrf, s/2 + s*(2*a-1), n);
    
    for b = 1:rsltn
    
      rrr = rttntn(rrr, s/2 + s, m);
      
      E(a, b) = nrgntn(rrr, rrrrf, qq, tp, aa, bb, dd, nbndd);
  
    end
    
  end

end