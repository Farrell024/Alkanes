function [dd, crbn] = plthlndhdrl(rrr)

    dbnds =  (length(rrr) - 2)/3 - 3;
    dd = zeros(1, dbnds);
    
    if ( length(rrr) > 23)
        
        rrrtmp = zeros(3, 7);
        
        rrrtmp(:, 1:3) = rrr(:, end-2:end); 

        rrrtmp(:, 4:7) = rrr(:, end-6:end-3);
        
        rrr(:, end-6:end) = rrrtmp;
        
    end
    
    crbn = zeros(1, dbnds+3);
    crbn(1) = 1;
    crbn(2) = 5;
    
    for j = 3:dbnds+3
        
       crbn(j) = crbn(j-1) + 3; 
    
    end
    
    for j = 1:dbnds
        
        v1 = cross(rrr(:, crbn(j+1) ) - rrr(:, crbn(j) ), rrr(:, crbn(j+2) ) - rrr(:, crbn(j+1) ) );
        v1 = -v1/norm(v1);
        
        v2 = cross(rrr(:, crbn(j+2) ) - rrr(:, crbn(j+1) ), rrr(:, crbn(j+2) ) - rrr(:, crbn(j+3) ) );
        v2 = v2/norm(v2);
        
        vx = rrr(:, crbn(j+2) ) - rrr(:, crbn(j+1) );
        vx = vx/norm(vx);
        
        sgn = dot(vx, cross(v1, v2) )/norm(dot(vx, cross(v1, v2) ));
        csn = dot(v1, v2);
        
        dd(1, j) = sgn*acosd(csn);
       
            
    end
    
end
