function [rrr, tp, qq, mm, aa, bb, dd, nbndd, pdb, psf] = bldr()

    mlcl = fopen('25.txt');
    mlclst= fopen('25strctr.txt');     % Opens files used for sizing
 
    mlcl2 = fopen('25.txt');
    mlclst2 = fopen('25strctr.txt');   % Opens files used for reading
    
%    mlcl = fopen('hept.txt');
%    mlclst= fopen('heptstrctr.txt');
 
%    mlcl2 = fopen('hept.txt');
%    mlclst2 = fopen('heptstrctr.txt');  
   
    chk = 0;                          % used to see if end of file has been reached
    
    N = 0;                            % Size of matrix read characters from psf/pdb
    MX = 0;                           % Candidate for above
   
    while (chk < 3)                   
    
        d = fgetl(mlcl2);
        d2 = fgetl(mlclst2);        
        
        N = N + 1;

        % Read through file line by line until end, keeping track of longest line, and number of lines
        
        if ( d ~= -1 | d2 ~= -1 )

            if ( MX < max( length(d), length(d2) ) )
        
                MX = max( length(d), length(d2) );
            
            end
    
            chk = 0;
    
        else 
        
            chk = chk + 1;
        
        end
    
    end
    
   pdb = zeros( max(N, MX) + 2 );
   psf = zeros( max(N, MX) + 2 );   % Preallocated matrices for both files that are sufficiently large
    
    j = 0;
    k = 0;    
    
    chk = 0;
    
    while (chk < 3)
    
        d = fgetl(mlcl);
        d2 = fgetl(mlclst);        
        
        if ( d ~= -1 | d2 ~= -1 )

            if (d ~= -1 ) 
            
                j = j + 1;
        
                for ( n = 1:length(d) )
        
                    pdb(j, n) = d(n);
                
                end
        
            end
            
            if ( d2 ~= -1)
               
                k = k + 1; 
              
                for ( n = 1:length(d2) )
        
                    psf(k, n) = d2(n);
                
                end
               
            end
    
            chk = 0;
    
        else 
        
            chk = chk + 1;
        
        end
    
    end
    
    N = 1;
    
    while ( pdb(N, 1) ~= 0 )
    
    N = N + 1;
    
    end
    
    N = N - 2;
    
    rrr = zeros(3, N);
    
    for ( n = 1:N )
    
        j = 27;
    
        while ( pdb(n + 1, j) == 32 )
        
            j = j + 1;
            
        end
        
        k = 1;
        
        while ( pdb(n + 1, j) ~= 32 )
        
            x(k) = pdb(n + 1, j);
       
            k = k + 1;  
            
            j = j + 1;
    
        end
 
   %         if (n == 152)
                
    %             qrty = str2double(char(x))
                
     %       end
                        
        while ( pdb(n + 1, j) == 32 )
        
            j = j + 1;
            
        end
        
        k = 1;
        
        while ( pdb(n + 1, j) ~= 32 )
        
            y(k) = pdb(n + 1, j);
            
            k = k + 1; 
            
            j = j + 1;
    
        end
  %          if (n == 152)
                
 %               y
          
%           end
            
        while ( pdb(n + 1, j) == 32 )
        
            j = j + 1;
            
        end
        
        k = 1;
        
        while ( pdb(n + 1, j) ~= 32 )
        
            z(k) = pdb(n + 1, j);

            
            k = k + 1;        
            j = j + 1;
    
        end
            
      %      if (n == 152)
                
     %           z
                
    %        end
            
        rrr(1, n) = str2double(char(x));
        rrr(2, n) = str2double(char(y));
        rrr(3, n) = str2double(char(z));
        
 %       if n == 152 
            
  %          rrr(1, n) = qrty;
            
   %     end        
    
    end
    
    ooo = rrr(:, 1);
    
    for ( j = 1:length(rrr) )
        
        rrr(:, j) = rrr(:, j) - ooo;
        
    end
    
    nbndd = zeros(N);
    
    for ( j = 1:N )
    
        for ( k = (j + 1):N )
        
            nbndd(j, k) = 1;
            
        end
    
    end
    
    l = 48;
    
    for ( j = 5:(4 + N) ) 
    
        k = 48;
        
        clear tt;
        
        while ( psf(j, k) ~= 32 )
            
            tt(k - 47) = psf(j, k);
            
            k = k + 1;
            
            if ( k > l )
            
                l = k;
                
            end
           
        end
        
        tp{j - 4} = char(tt);
        
    end
    
    for ( j = 5:(4 + N) ) 
    
        k = l;
        
        while ( psf(j, k) == 32) 
        
            k = k + 1;
            
        end
        
        h = k;
        
        clear tt;
        
        while ( psf(j, k) ~= 32 )
            
            tt(k - h + 1) = psf(j, k);
            
            k = k + 1;
           
        end
        
        qq(j - 4) = str2double(char(tt));
        
    end    
    
    l = k;
    
    for ( j = 5:(4 + N) ) 
    
        k = l;
        
        while ( psf(j, k) == 32) 
        
            k = k + 1;
            
        end
        
        h = k;
        
        clear tt;
        
        while ( psf(j, k) ~= 32 )
            
            tt(k - h + 1) = psf(j, k);
            
            k = k + 1;
           
        end
        
        mm(j - 4) = str2double(char(tt));
        
    end  
    
    j = j + 3;
    
    k = 1;
    
    while ( psf(j, k) == 32 )
    
        k = k + 1;
        
    end   
    
    h = 1;

    while ( psf(j, k) > 32 )

        while ( psf(j, k) > 32 )
        
            l = 1;
        
            clear nn;
        
            while ( psf(j, k + l - 1) > 32 )
        
                nn(l) = psf(j, k + l - 1);
                
                l = l + 1;
            
            end
        
            bb(1, h) = str2double( char(nn) );
            
            k = k + l;
            
            while ( psf(j, k) == 32 )
    
                k = k + 1;
        
            end 
            
            l = 1;
            
            clear nn;
        
            while ( psf(j, k + l - 1) > 32 )
        
                nn(l) = psf(j, k + l - 1);
                
                l = l + 1;
            
            end
            
            bb(2, h) = str2double( char( nn ) );
            
            h = h + 1;
            k = k + l;

            while ( psf(j, k) == 32 )
            
                k = k + 1;
                
            end 
            
        end
        
        k = 1;
        j = j + 1;
        
        while ( psf(j, k) == 32 )
        
            k = k + 1;
            
        end 
        
    end
    
    j = j + 2;
    h = 1;
    k = 1;
    
    while ( psf(j, k) == 32 )
    
        k = k + 1;
        
    end  
    
    while ( psf(j, k) > 32 )

        while ( psf(j, k) > 32 )
        
            l = 1;
        
            clear nn;
        
            while ( psf(j, k + l - 1) > 32 )
        
                nn(l) = psf(j, k + l - 1);
                
                l = l + 1;
            
            end
        
            aa(1, h) = str2double( char(nn) );
            
            k = k + l;
            
            while ( psf(j, k) == 32 )
    
                k = k + 1;
        
            end 

            l = 1;
        
            clear nn;
        
            while ( psf(j, k + l - 1) > 32 )
        
                nn(l) = psf(j, k + l - 1);
                
                l = l + 1;
            
            end
        
            aa(2, h) = str2double( char(nn) );
            
            k = k + l;
            
            while ( psf(j, k) == 32 )
    
                k = k + 1;
        
            end 
            
            l = 1;
            
            clear nn;
        
            while ( psf(j, k + l - 1) > 32 )
        
                nn(l) = psf(j, k + l - 1);
                
                l = l + 1;
            
            end
            
            aa(3, h) = str2double( char( nn ) );
            
            h = h + 1;
            k = k + l;

            while ( psf(j, k) == 32 )
            
                k = k + 1;
                
            end 
            
        end
        
        k = 1;
        j = j + 1;
        
        while ( psf(j, k) == 32 )
        
            k = k + 1;
            
        end 
        
    end    

    j = j + 2;
    h = 1;
    k = 1;
    
    while ( psf(j, k) == 32 )
    
        k = k + 1;
        
    end  
    
    while ( psf(j, k) > 32 )

        while ( psf(j, k) > 32 )
        
            l = 1;
        
            clear nn;
        
            while ( psf(j, k + l - 1) > 32 )
        
                nn(l) = psf(j, k + l - 1);
                
                l = l + 1;
            
            end
        
            dd(1, h) = str2double( char(nn) );
            
            k = k + l;
            
            while ( psf(j, k) == 32 )
    
                k = k + 1;
        
            end 

            l = 1;
        
            clear nn;
        
            while ( psf(j, k + l - 1) > 32 )
        
                nn(l) = psf(j, k + l - 1);
                
                l = l + 1;
            
            end
        
            dd(2, h) = str2double( char(nn) );
            
            k = k + l;
            
            while ( psf(j, k) == 32 )
    
                k = k + 1;
        
            end 

            l = 1;
        
            clear nn;
        
            while ( psf(j, k + l - 1) > 32 )
        
                nn(l) = psf(j, k + l - 1);
                
                l = l + 1;
            
            end
        
            dd(3, h) = str2double( char(nn) );
            
            k = k + l;
            
            while ( psf(j, k) == 32 )
    
                k = k + 1;
        
            end
            
            l = 1;
            
            clear nn;
        
            while ( psf(j, k + l - 1) > 32 )
        
                nn(l) = psf(j, k + l - 1);
                
                l = l + 1;
            
            end
            
            dd(4, h) = str2double( char( nn ) );
            
            h = h + 1;
            k = k + l;

            while ( psf(j, k) == 32 )
            
                k = k + 1;
                
            end 
            
        end
        
        k = 1;
        j = j + 1;
        
        while ( psf(j, k) == 32 )
        
            k = k + 1;
            
        end
        
    end
   
    for n = 1:length(dd)
 
        nbndd( dd(1, n), dd(2, n) ) = 0;
        nbndd( dd(1, n), dd(3, n) ) = 0;
        nbndd( dd(1, n), dd(4, n) ) = 0;
        nbndd( dd(2, n), dd(3, n) ) = 0;
        nbndd( dd(2, n), dd(4, n) ) = 0;
        nbndd( dd(3, n), dd(4, n) ) = 0;
        nbndd( dd(2, n), dd(1, n) ) = 0;
        nbndd( dd(3, n), dd(1, n) ) = 0;
        nbndd( dd(4, n), dd(1, n) ) = 0;
        nbndd( dd(3, n), dd(2, n) ) = 0;
        nbndd( dd(4, n), dd(2, n) ) = 0;
        nbndd( dd(4, n), dd(3, n) ) = 0;
    
    end
    
    for n = 1:length(aa)
 
        nbndd( aa(1, n), aa(2, n) ) = 0;
        nbndd( aa(1, n), aa(3, n) ) = 0;
        nbndd( aa(2, n), aa(3, n) ) = 0;
        nbndd( aa(2, n), aa(1, n) ) = 0;
        nbndd( aa(3, n), aa(1, n) ) = 0;
        nbndd( aa(3, n), aa(2, n) ) = 0;
    
    end
    
end
    
    