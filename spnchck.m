function F = spnchck(rrr, n, m, k, l, rsltn)

  rrrrf = rrr;
  
  s = 360/(rsltn+1);
  
  F = zeros(rsltn);
  
  for a = 1:rsltn
  
    rrr = rttntn(rrrrf, s/2 + s*a, n);
    
    rrrx = rrr;
    
    for b = 1:rsltn
    
      rrr = rttntn(rrrx, s/2 + s*b, m);
      
      F(a, b) = tbchck(rrrrf, rrr, k, l);
  
    end
    
  end

end