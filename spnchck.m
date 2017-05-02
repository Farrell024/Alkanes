function [F, Q] = spnchck(rrr, n, m, k, l, rsltn)

  rrrrf = rrr;
  
  s = 360/(rsltn);
  
  for a = 1:rsltn
  
    rrr = rttntn(rrr, s/2*(2*a-1), n);
    
    for b = 1:rsltn
    
      rrr = rttntn(rrr, s/2*(2*b-1), m);
      
      F(a, b) = tbchck(rrrrf, rrr, k, l);
      
      Q(:, a, b) = rrr(:, 150);
  
    end
    
  end

end