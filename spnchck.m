function F = spnchck(rrr, n, m, k, l)

  rrrrf = rrr;
  
  rsltn = 36;
  
  s = 360/(rsltn + 1);
  
  rrr = rttntn(rrr, s/2, n);
  rrr = rttntn(rrr, s/2, m);
  
  for a = 1:rsltn
  
    rrr = rttntn(rrr, s, n);
    
    for b = 1:rsltn
    
      rrr = rttntn(rrr, s, m);
      
      F(a, b) = tbchck(rrrrf, rrr, k, l);
  
    end
    
  end

end