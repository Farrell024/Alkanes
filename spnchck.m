function F = spnchck(rrr, n, m, k, l, rsltn)

  rrrrf = rrr;
  
  s = 360/(rsltn);
  
  F = zeros(rsltn);
  
  for a = 1:rsltn
  
    rrr = rttntn(rrrrf, s/2 + s*(2*a-1), n);
    
    for b = 1:rsltn
    
      rrr = rttntn(rrr, s/2 + s, m);
      
      F(a, b) = tbchck(rrrrf, rrr, k, l);
  
    end
    
  end

end