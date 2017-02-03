function pltntn(rrr)

    xx = rrr(1,:);
    yy = rrr(2,:);
    zz = rrr(3,:);
   
    xxcc(1) = xx(1);
    yycc(1) = yy(1);
    zzcc(1) = zz(1);
 
    xxhh(1) = xx(2);
    yyhh(1) = yy(2); 
    zzhh(1) = zz(2);
    xxhh(2) = xx(3);
    yyhh(2) = yy(3);
    zzhh(2) = zz(3);
    xxhh(3) = xx(4);
    yyhh(3) = yy(4);
    zzhh(3) = zz(4);
    
    N = (length(rrr)-2)/3;
    
    m = 5;  
   
    if N == 1
        
        xxhh(4) = xx(5);
        yyhh(4) = yy(5); 
        zzhh(4) = zz(5);  
        
    end         
    
    for n = 2:N
    
        xxcc(n) = xx(m);
        yycc(n) = yy(m);
        zzcc(n) = zz(m);
        m = m + 1;
        
        if n == N
        
            xxhh(m - n) = xx(m);
            yyhh(m - n) = yy(m); 
            zzhh(m - n) = zz(m);
            xxhh(m - n + 1) = xx(m+1);
            yyhh(m - n + 1) = yy(m+1);
            zzhh(m - n + 1) = zz(m+1);
            xxhh(m - n + 2) = xx(m+2);
            yyhh(m - n + 2) = yy(m+2);
            zzhh(m - n + 2) = zz(m+2);
        
        else

            xxhh(m - n) = xx(m);
            yyhh(m - n) = yy(m); 
            zzhh(m - n) = zz(m);
            xxhh(m - n + 1) = xx(m+1);
            yyhh(m - n + 1) = yy(m+1);
            zzhh(m - n + 1) = zz(m+1);
        
        end
       
       m = m + 2;
    
    end    
       
    plot3(xxcc,yycc,zzcc,'-',xxhh,yyhh,zzhh,'or');
    
    L = norm( zz( length(zz) ) - zz(2) );
    
    L = fix( L + 1 );
    
    axis([-(L+2)/2, (L+2)/2, -2, L, -2, L]);    

end