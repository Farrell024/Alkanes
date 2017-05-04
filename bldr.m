function [rrr, tp, qq, mm, aa, bb, dd, nbndd, pdb, psf] = bldr()

    mlcl = fopen('hept.txt');
    mlclst= fopen('heptstrctr.txt');     % Opens files used for sizing
 
    mlcl2 = fopen('hept.txt');
    mlclst2 = fopen('heptstrctr.txt');   % Opens files used for reading
    

%    mlcl = fopen('PE25.pdb');
%    mlclst = fopen('PE25.txt');
 
%    mlcl2 = fopen('PE25.pdb');
%    mlclst2 = fopen('PE25.txt');  
   
    chk = 0;                          % used to see if end of file has been reached
    
    N = 0;                            % Size of matrix read characters from psf/pdb
    MX = 0;                           % Candidate for above
   
    while (chk < 3)                   
    
        d = fgetl(mlcl2);
        d2 = fgetl(mlclst2);        
        
        N = N + 1;

        % Read through file line by line until end, keeping track of longest line, and number of lines
        
        if ( sum(d ~= -1) || sum(d2 ~= -1) )

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
    
    while (chk < 3)                 %read through the file but now add the character values to the pdb/psf object
    
        d = fgetl(mlcl);
        d2 = fgetl(mlclst);        
        
        if ( sum(d ~= -1) || sum(d2 ~= -1) )

            if (d ~= -1 ) 
            
                j = j + 1;
        
                for  n = 1:length(d) 
        
                    pdb(j, n) = d(n);
                
                end
        
            end
            
            if ( d2 ~= -1)
               
                k = k + 1; 
              
                for  n = 1:length(d2) 
        
                    psf(k, n) = d2(n);
                
                end
               
            end
    
            chk = 0;
    
        else 
        
            chk = chk + 1;
        
        end
    
    end
    
    N = 1;
    
    while ( pdb(N, 1) ~= 0 ) % Sizes up the position matrix (rrr)
    
    N = N + 1;
    
    end
    
    N = N - 2;      % Elimates the line in pdb which is unrelated to the atoms and the +1 overcount from the loop
    
    rrr = zeros(3, N);
    x = zeros(1, N);
    z = zeros(1, N);
    y = zeros(1, N);
    
    for  n = 1:N 
    
        j = 27;
    
        while ( pdb(n + 1, j) == 32 ) % n+1 because the first line does not have atom information
        
            j = j + 1;      % charge j up to the colum where the x-coordinate values start
            
        end
        
        k = 1;      % index for character values
        
        while ( pdb(n + 1, j) ~= 32 )   % read the nth x coordinate
        
            x(k) = pdb(n + 1, j);   
       
            k = k + 1;  
            
            j = j + 1;
    
        end
                        
        while ( pdb(n + 1, j) == 32 ) % charge j to y column
        
            j = j + 1;
            
        end
        
        k = 1;
        
        while ( pdb(n + 1, j) ~= 32 ) % read the nth y coordinate 
        
            y(k) = pdb(n + 1, j);
            
            k = k + 1; 
            
            j = j + 1;
    
        end
            
        while ( pdb(n + 1, j) == 32 ) % charge j to z column
        
            j = j + 1;
            
        end
        
        k = 1;
        
        while ( pdb(n + 1, j) ~= 32 ) % read the nth z coordinate
        
            z(k) = pdb(n + 1, j);

            k = k + 1;        
            
            j = j + 1;
    
        end
          
        rrr(1, n) = str2double(char(x));          % fill the three coordinates of the nth atom  
        rrr(2, n) = str2double(char(y));
        rrr(3, n) = str2double(char(z));
        
    end
    
    ooo = rrr(:, 1);    % the origin position
    
    for  j = 1:length(rrr)  % shift everything so atom 1 sits at <0, 0, 0>
        
        rrr(:, j) = rrr(:, j) - ooo;
        
    end
    
    nbndd = zeros(N);      % insatiate the nonbonded matrix of atom pair
    
    for  j = 1:N       %fill it in upper triangular
    
        for  k = (j + 1):N 
        
            nbndd(j, k) = 1; 
            
        end
    
    end
    
    l = 48; % charge the column position to the atom type column
    
    bb = zeros(1,1);
    aa = zeros(1,1);
    dd = zeros(1,1);
    nn = zeros(1,1);
    mm = zeros(1,N);
    qq = zeros(1,N);
    tt = zeros(1,1);
    tp = cell(1,N);
    
    for  j = 5:(4 + N)  % Start at line 5 and continue for the length reading in the type names into tt
    
        k = 48;
        
        clear tt;
        
        while ( psf(j, k) ~= 32 ) % read until end of word

                tt(k - 47) = psf(j, k);

                k = k + 1;

                if ( k > l )

                    l = k; % l sees when the column for type name stops, given the longest name

                end

        end
        
        tp{j - 4} = char(tt);
        
    end
    
    for  j = 5:(4 + N)  
    
        k = l; % start at the end of the type name column
        
        while ( psf(j, k) == 32)  % move up to the charge column
        
            k = k + 1; 
            
        end
        
        h = k;
        
        clear tt;
        
        while ( psf(j, k) ~= 32 ) % read the charge values 
            
            tt(k - h + 1) = psf(j, k);
            
            k = k + 1;
           
        end
        
        qq(j - 4) = str2double(char(tt));
        
    end    
    
    l = k;
    
    for  j = 5:(4 + N)  % read in the mass column
    
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
    
    j = j + 3; % move down to the directly bonded group blocked
    
    k = 1;
    
    while ( psf(j, k) == 32 ) %find the first value
    
        k = k + 1;
        
    end   
    
    h = 1;

    while ( psf(j, k) > 32 )

        while ( psf(j, k) > 32 ) %run until it falls of the end of the line
        
            l = 1;
        
            clear nn;
        
            while ( psf(j, k + l - 1) > 32 ) % read the value into nn
        
                nn(l) = psf(j, k + l - 1);
                
                l = l + 1;
            
            end
        
            bb(1, h) = str2double( char(nn) ); % put the value as the first of the ordered pair
            
            k = k + l;
            
            while ( psf(j, k) == 32 ) % move to the next value
    
                k = k + 1;
        
            end 
            
            l = 1;
            
            clear nn;
        
            while ( psf(j, k + l - 1) > 32 ) % read the value 
        
                nn(l) = psf(j, k + l - 1);
                
                l = l + 1;
            
            end
            
            bb(2, h) = str2double( char( nn ) ); % put the value as the second entry of the same ordered pair
            
            h = h + 1; % increase the ordered pair index
            k = k + l; % move past the value last read

            while ( psf(j, k) == 32 ) % move to next equation
            
                k = k + 1;
                
            end 
            
        end
        
        k = 1; % start at the left edge
        j = j + 1; % drop a line down
        
        while ( psf(j, k) == 32 ) % move over to the first value
        
            k = k + 1;
            
        end 
        
    end
    
    j = j + 2; % drop down to the angle bond block and start at the left edge
    h = 1;
    k = 1;
    
    while ( psf(j, k) == 32 ) % move up to the first value
    
        k = k + 1;
        
    end  
    
    while ( psf(j, k) > 32 ) % this code block shoul follow a parallel procedure to the one above

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
    
    while ( psf(j, k) > 32 ) % This block should follow a parallel procedure to the one above

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
   
    for n = 1:length(dd) % eliminate pairs from the non-bonded group belonging to dd
 
%        nbndd( dd(1, n), dd(2, n) ) = 0;
%        nbndd( dd(1, n), dd(3, n) ) = 0;
        nbndd( dd(1, n), dd(4, n) ) = 2;
%        nbndd( dd(2, n), dd(3, n) ) = 0;
%        nbndd( dd(2, n), dd(4, n) ) = 0;
%        nbndd( dd(3, n), dd(4, n) ) = 0;
%        nbndd( dd(2, n), dd(1, n) ) = 0;
%        nbndd( dd(3, n), dd(1, n) ) = 0;
%        nbndd( dd(4, n), dd(1, n) ) = 0;
%        nbndd( dd(3, n), dd(2, n) ) = 0;
%        nbndd( dd(4, n), dd(2, n) ) = 0;
%        nbndd( dd(4, n), dd(3, n) ) = 0;
    
    end
    
    for n = 1:length(aa) % elimate pairs from the non-bonded group belonging to aa this would include bb
 
        nbndd( aa(1, n), aa(2, n) ) = 0;
        nbndd( aa(1, n), aa(3, n) ) = 0;
        nbndd( aa(2, n), aa(3, n) ) = 0;
        nbndd( aa(2, n), aa(1, n) ) = 0;
        nbndd( aa(3, n), aa(1, n) ) = 0;
        nbndd( aa(3, n), aa(2, n) ) = 0;
    
    end
    
end
    
    
