function rrrtt = rttntn(rrr, p, b)

    % do nothing if alkane is too small to rotate
  
    if ( length(rrr) == 4 )
    
        rrrtt = rrr;
        
        return;
        
    end
    
   % correct polymer ordering
    
    if ( length(rrr) > 14)
        
        rrrtmp = zeros(3, 7);
        
        rrrtmp(:, 1:3) = rrr(:, end-2:end); 
        rrrtmp(:, 4:7) = rrr(:, end-6:end-3);
        
        rrr(:, end-6:end) = rrrtmp;
        
    end
    

    % determine length of alkane and determine the atom number characteristic of the bond of interest
    
    N = length(rrr);
    b = fix(b);
    b = 1 + mod(b-1, (N-2)/3-1); 
    M = 1 + 3*b;
   
    rrrtt = zeros(3,M);
    
    for j = 1:M % for atoms beneath the atom before the bond of interest, do not rotate
   
        rrrtt(1,j) = rrr(1,j); 
        rrrtt(2,j) = rrr(2,j);
        rrrtt(3,j) = rrr(3,j);  
    
    end
    
    kr  = 0;
    
    if ( M == 4 )
    
        kr = -1;
    
    end

    %create two matrices to rotate around a general direction as determined by the bond of interest
    
    zp = rrr(:, M+1) - rrr(:, M-2+kr);    
    p1 = atan2( zp(1), zp(2) );    
    t1 = [cos(p1), -sin(p1), 0; sin(p1), cos(p1), 0; 0, 0, 1];
    zpp = t1*zp;    
    p2 = atan2( zpp(2), zpp(3) );
    t2 = [1, 0, 0; 0, cos(p2), -sin(p2); 0, sin(p2), cos(p2)];
    
    t = [cosd(p), -sind(p), 0; sind(p), cosd(p), 0; 0, 0, 1];
    
    RT = (t1^-1)*(t2^-1)*t*(t2)*(t1);
    
    for j = M+1:N
        
        rrrtt(:,j) = rrr(:, M-2+kr) +  RT*( rrr(:,j) - rrr(:, M-2+kr) ); 
        
    end
    
    if ( length(rrr) > 14)
        
        rrrtmp(:, 1:4) = rrrtt(:, end-3:end); 
        rrrtmp(:, 5:7) = rrrtt(:, end-6:end-4);
        
        rrrtt(:, end-6:end) = rrrtmp;
        
    end    
    
end
