function f = tbchck(rrrrf, rrr, n, m)

    R = 2.46/pi*sqrt(n*n + n*m + m*m)/2;
    
    a = ( rrrrf(:,1) + rrrrf(:,4) )/2;
    b = ( rrrrf(:,length(rrrrf) - 6) + rrrrf(:,length(rrrrf) - 3) )/2; 
 
    axs = (b - a)/norm(b - a);
   
    f = 1;
   
    for n = 1:length(rrr)
   
        if norm( rrr(:, n) - axs*dot(axs, rrr(:,n) ) ) > R
   
            f = 0;   
       
        end
        
    end 

end