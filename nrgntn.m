function E = nrgntnr(rrr, rrrrf, qq, tp, aa, bb, dd, nbndd)

    E = nrgntnrw(rrr, qq, tp, aa, bb, dd, nbndd);
   % - nrgntnrw(rrrrf, qq, tp, aa, bb, dd, nbndd);

end

function E = nrgntnrw(rrr, qq, tp, aa, bb, dd, nbndd)

    E = 0;
    
    N = length(rrr);
   
   for j = 1:N

        for k = j+1:N
   
            r = norm( rrr(:,j) - rrr(:,k) );
   
            if ( nbndd(j, k) )                      
                  
               %E = E + LJ( r, tp(j), tp(k) );              
               % E(j, k) = C( r, qq(j) , qq(k) );
               %E = E + C( r, qq(j) , qq(k) );
                
            end
                                    
             %E = E + C( r, qq(j) , qq(k) );
            
        end
        
    end

%    r = norm( rrr(:, 7) - rrr(:, 15) );
%    E = C( r, qq(7), qq(15) );   
    
    
    if ( N == 4 )
    
        return;
        
    end    
    
    for n = 1:length(dd)
    
 %       E = E + dhdrl( dd(1, n), dd(2, n), dd(3, n), dd(4, n), tp, rrr);
    
    end

    for n = 1:length(bb)

  %      E = E + bnd( bb(1, n), bb(2, n), tp, rrr);
    
    end    

    for n = 1:length(aa)
    
       aa(1, n) 
       aa(2, n)
       aa(3, n)
       
       E = E + ngl( aa(1, n), aa(2, n), aa(3, n), tp, rrr)
    
    end
    
end

function U = LJ(r, t1, t2)

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

    U = nrg*( (rm/r)^12 - 2*(rm/r)^6 );    


end

function U = C(r, q1, q2)

    U = (331.6*q1*q2)/r;

end

function A = ngl(m1, m2, m3, tp, rrr)

    r1 = rrr(:, m1) - rrr(:, m2);
    r2 = rrr(:, m3) - rrr(:, m2);

    D = dot(r1, r2);
    M = norm(r1)*norm(r2);

    T = abs( acosd(D/M) )

    if ( strcmp( tp(m1), 'CC32A' ) )

        if ( strcmp( tp(m3), 'CC32A' ) )

            A = 58.53*( T - 113.6 )^2;
            
            T - 113.6
            
        elseif ( strcmp( tp(m3), 'CC33A' ) )

            A = 58*( T - 115 )^2;
            
            T - 115
            
        elseif ( strcmp( tp(m3), 'HCA3' ) ) 

            A = 26*( T - 110.1 )^2;

            T - 110.1
            
        else
       
            A = 26.5*( T - 110.1 )^2;
        
            T - 110.1
            
        end

    elseif ( strcmp( tp(m1), 'CC33A' ) )

        if ( strcmp( tp(m3), 'CC32A' ) )

            A = 58*( T - 115 )^2;
            
            T - 115
            
        elseif ( strcmp( tp(m3), 'CC33A' ) )

            A = 53.35*( T - 114 )^2;
            
            T - 114
            
        elseif ( strcmp( tp(m3), 'HCA3' ) ) 

            A = 37.5*( T - 110.1 )^2;

            T - 110.1
            
        else
       
            A = 34.6*( T - 110.1 )^2;
            
            T - 110.1
        
        end

    elseif ( strcmp( tp(m1), 'HCA3' ) )

        if ( strcmp( tp(m3), 'CC32A' ) )

            A = 34.6*( T - 110.1 )^2;
            
            T - 110.1
            
        elseif ( strcmp( tp(m3), 'CC33A' ) )

            A = 37.5*( T - 110.1 )^2;
            
            T - 110.1
            
        else

            A = 35.5*( T - 108.4 )^2;
        
            T - 108.4
            
        end
        
    else

        if ( strcmp( tp(m3), 'CC32A' ) )

            A = 26.5*( T - 110.1 )^2;
            
            T - 110.1
            
        elseif ( strcmp( tp(m3), 'CC33A' ) )

            A = 34.6*( T - 110.1 )^2;
            
            T - 110.1
            
        else

            A = 35.5*( T - 109 )^2;
        
            T - 109
            
        end
        
    end
        
end

function B = bnd(m1, m2, tp, rrr)

    if ( strcmp( tp(m1), 'CC32A' ) )
    
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

    v1 = cross( rrr(:, m2) - rrr(:, m1), rrr(:, m3) - rrr(:, m2) );
    v2 = cross( rrr(:, m3) - rrr(:, m2), rrr(:, m4) - rrr(:, m3) );

    c = dot(v1, v2)/( norm(v1) * norm(v2) );    
    p = acosd(c);
    
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
    
    D = sum( k.*( 1 + cosd(n.*p - t) ) );
    
end
