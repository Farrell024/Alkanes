function E = nrgntn(rrr, rrrrf, qq, tp, aa, bb, dd, nbndd)

    E = nrgntnrw(rrr, qq, tp, aa, bb, dd, nbndd); 
% - nrgntnrw(rrrrf, qq, tp, aa, bb, dd, nbndd); % return the energy from all the contributions, optionally with respect to a reference energy level

end

function E = nrgntnrw(rrr, qq, tp, aa, bb, dd, nbndd) % manages the calls to individual energy components

    E = 0;
    
    N = length(rrr);
   
   for j = 1:N

        for k = j+1:N % go through all the pairs and calculate LJ and C given they are nonbonded
   
            r = norm( rrr(:,j) - rrr(:,k) );
   
            if ( nbndd(j, k) )                      
                  
               %E = E + LJ( r, tp(j), tp(k) );              
               %E = E + C( r, qq(j) , qq(k) );
                
            end
                                    
       end
        
    end    
    
    if ( N == 4 ) % if it 4 atom molecule WHO CARES
    
        return;
        
    end    
    
    for n = 1:length(dd) % do all the dihedrals
    
 %       E = E + dhdrl( dd(1, n), dd(2, n), dd(3, n), dd(4, n), tp, rrr);
    
    end

    for n = 1:length(bb) % do all the direct bonds

  %      E = E + bnd( bb(1, n), bb(2, n), tp, rrr);
    
    end    

    for n = 1:length(aa) % do all the angled bonds

        E = E + ngl( aa(1, n), aa(2, n), aa(3, n), tp, rrr)
    
    end
    
end

function U = LJ(r, t1, t2)

    %Identify the types, and assign parameters
			       
    if ( strcmp(t1, 'CC33A') )
        
        nrgj = .078;
        
        rmj = 2.08; 
        
    elseif ( strcmp(t1, 'HCA3') )
    
        nrgj = .024;
        
        rmj = 1.34;
        
    elseif ( strcmp(t1, 'CC32A') )
    
        nrgj = .056;
        
        rmj = 2.01;
        
    else
    
        nrgj = .035;
    
        rmj = 1.34;
        
    end
    
    if ( strcmp(t2, 'CC33A') )
        
        nrgk = .078;
        
        rmk = 2.08; 
        
    elseif ( strcmp(t2, 'HCA3') )
    
        nrgk = .024;
        
        rmk = 1.34;
        
    elseif ( strcmp(t2, 'CC32A') )
    
        nrgk = .056;
        
        rmk = 2.01;
        
    else
    
        nrgk = .035;
    
        rmk = 1.34;
        
    end
    
    nrg = sqrt(nrgk*nrgj);
    rm = (rmk+rmj);

    U = nrg*( (rm/r)^12 - 2*(rm/r)^6 ); % calculate   

end

function U = C(r, q1, q2)

    U = (331.6*q1*q2)/r;

end

function A = ngl(m1, m2, m3, tp, rrr)

    % calculate angle
  
    r1 = rrr(:, m1) - rrr(:, m2);
    r2 = rrr(:, m3) - rrr(:, m2);

    D = dot(r1, r2);
    M = norm(r1)*norm(r2);

    T = abs( acos(D/M) )

    if ( strcmp( tp(m1), 'CC32A' ) ) % determine type and assign paratemers 

        if ( strcmp( tp(m3), 'CC32A' ) )

            A = 58.53*( T - 113.6/180*3.1415926536 )^2;
            
        elseif ( strcmp( tp(m3), 'CC33A' ) )

            A = 58*( T - 115/180*3.1415926536 )^2;
                       
        elseif ( strcmp( tp(m3), 'HCA3' ) ) 

            A = 26*( T - 110.1/180*3.1415926535 )^2;
            
        else
       
            A = 26.5*( T - 110.1/180*3.1415926535 )^2;
                    
        end

    elseif ( strcmp( tp(m1), 'CC33A' ) )

        if ( strcmp( tp(m3), 'CC32A' ) )

            A = 58*( T - 115/180*3.1415926536 )^2;
                        
        elseif ( strcmp( tp(m3), 'CC33A' ) )

            A = 53.35*( T - 114/180*3.1415926535 )^2;
                        
        elseif ( strcmp( tp(m3), 'HCA3' ) ) 

	    A = 37.5*( T - 110.1/180*3.1415926536 )^2;
            
        else
       
            A = 34.6*( T - 110.1/180*3.1415926536 )^2;
                             
        end

    elseif ( strcmp( tp(m1), 'HCA3' ) )

        if ( strcmp( tp(m3), 'CC32A' ) )

            A = 34.6*( T - 110.1/180*3.1415926536 )^2;
                       
        elseif ( strcmp( tp(m3), 'CC33A' ) )

	    A = 37.5*( T - 110.1/180*3.1415926536 )^2;
                                   
        else

            A = 35.5*( T - 108.4/180*3.1415926536 )^2;
                             
        end
        
    else

        if ( strcmp( tp(m3), 'CC32A' ) )

            A = 26.5*( T - 110.1/180*3.1415926536 )^2;
            
        elseif ( strcmp( tp(m3), 'CC33A' ) )

            A = 34.6*( T - 110.1/180*3.1415926536 )^2;
            
        else

	    A = 35.5*( T - 109/180*3.1415926536 )^2;
            
        end
        
    end
        
end

function B = bnd(m1, m2, tp, rrr)

    if ( strcmp( tp(m1), 'CC32A' ) ) % read type and assign parameters
    
        if ( strcmp( tp(m2), 'CC32A' ) )
        
            B = 222.5*( norm( rrr(:, m1) - rrr(:, m2) ) - 1.53 )^2;
            
        elseif ( strcmp( tp(m2), 'CC33A' ) )
        
            B = 222.5*( norm( rrr(:, m1) - rrr(:, m2) ) - 1.528 )^2;
            
        else 
        
            B = 309.0*( norm( rrr(:, m1) - rrr(:, m2) ) - 1.111 )^2;
        
        end
        
    elseif ( strcmp( tp(m1), 'CC33A' ) )
    
        if ( strcmp( tp(m2), 'CC32A' ) )
        
            B = 222.5*( norm( rrr(:, m1) - rrr(:, m2) ) - 1.53 )^2;
                        
        else 
        
            B = 322.0*( norm( rrr(:, m1) - rrr(:, m2) ) - 1.111 )^2;
        
        end

    elseif ( strcmp( tp(m1), 'HCA2' ) )
        
        B = 309*( norm( rrr(:, m1) - rrr(:, m2) ) - 1.111 )^2;
    
    else

        B = 322*( norm( rrr(:, m1) - rrr(:, m2) ) - 1.111 )^2;    
        
    end
        
end

function D = dhdrl (m1, m2, m3, m4, tp, rrr)

      %calculate dihedral angle
      
    v1 = cross( rrr(:, m2) - rrr(:, m1), rrr(:, m3) - rrr(:, m2) );
    v2 = cross( rrr(:, m3) - rrr(:, m2), rrr(:, m4) - rrr(:, m3) );

    c = dot(v1, v2)/( norm(v1) * norm(v2) );    
    p = acosd(c);
    
 %read type and assign parameters

 if ( strcmp( tp(m1), 'CC33A' ) ) 
 
        if ( strcmp( tp(m4), 'CC32A' ) )
        
            k = [.20391, .10824, .08133, .15051];
            n = [5, 4, 3, 2];
            t = [0, 0, 180, 0];
        
        elseif ( strcmp( tp(m4), 'CC33A' ) )
        
            k = [.03179, .03819];
            n = [6, 2];
            t = [180, 0];       
        
        elseif ( strcmp( tp(m4), 'HCA3' ) )
        
            k = .16;
            n = 3;
            t = 0;
        
        else
        
            k = .19;
            n = 3;
            t = 0;        
        
        end
        
    elseif ( strcmp( tp(m1), 'CC32A' ) )
        
        if ( strcmp( tp(m4), 'CC32A' ) )
        
            k = [.11251, .09458, .14975, .06450];
            n = [5, 4, 3, 2];
            t = [0, 0, 180, 0];        
        
        elseif ( strcmp( tp(m4), 'CC33A' ) )
        
             k = [.20391, .10824, .08133, .15051];
             n = [5, 4, 3, 2];
             t = [0, 0, 180, 0];       
        
        elseif ( strcmp( tp(m4), 'HCA3' ) )
        
             k = .16;
             n = 3;
             t = 0;       
        
        else
        
             k = .19;
             n = 3;
             t = 0;         
        
        end
        
    elseif ( strcmp( tp(m1), 'HCA3' ) )

        k = .16;
        n = 3;
        t = 0;        

    else

        k = .19;
        n = 3;
        t = 0;       
        
        if ( strcmp( tp(m4), 'HCA3' ) )
        
            k = .16;
            n = 3;
            t = 0;         
        
        end
        
    end
    
    D = sum( k.*( 1 + cosd(n.*p - t) ) ); % calculate dihedral for all terms

end
