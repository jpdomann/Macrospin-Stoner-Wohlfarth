function N = DemagFactor(shape,Dims)
%DEMAGFACTOR calculates the demagnetization factors for ellipsoidal and
%rectangular shapes.
%   To do: enable the user to rotate the ellipse with a set of Euler angles

[Dims_sorted, Sort_ind] = sort(Dims,'descend');

a = Dims_sorted(1);
b = Dims_sorted(2);
c = Dims_sorted(3);
if (a <0 || b<0 || c<0) 
    warning('A negative dimension was entered')
    N = [0 0 0];
    return
end

switch shape
    case 'Ellipse'
        if (a == b) && (b == c) %sphere
            N1 = 1/3;
            N2 = 1/3;
            N3 = 1/3;
        elseif (a==b) && a>c     %oblate spheroid
            e = sqrt((a^2-c^2)/(a^2));
            N3 =1/e^2 .* (1 - sqrt(1-e^2)/e *asin(e));
            N1 = 0.5 .* (1-N3);
            N2 = N1;
        elseif (a == b) && (a < c)  %prolate spheroid
            e = sqrt((c^2-a^2)/(c^2));
            N3 =(1-e^2)/e^2 .* (1/(2*e).*log((1+e)/(1-e)) -1 );
            N1 = 0.5 .* (1-N3);
            N2 = N1;
        elseif (a > b) && (b == c)  %prolate spheroid
            m = a/c;
            n = m^2-1;
            N1 = 1/n * ( m /(2*sqrt(n)) * log((m+sqrt(n))/(m-sqrt(n))) - 1 );
            N2 = m/(2*n) * (m - 1/(2*sqrt(n))* log((m+sqrt(n))/(m-sqrt(n))) );
            N3 = N2; 
        elseif (a >= b) && (b >= c)                        
            %"Depolartization Tensor Method - 2nd Ed."
            B = b/a;    B2 = 1-B^2;
            G = c/a;    G2 = 1-G^2;                         
            k = sqrt(B2 / G2 )^2; 
            phi = acos(G);
            N1 = B*G /( B2*sqrt(G2) )*( ellipticF(phi,k) - ellipticE(phi,k) );
            N2 = - G^2/(B^2 - G^2) +...
                B*G*sqrt(G2) / ( B2*(B^2-G^2) ) * ellipticE(phi,k) -...
                B*G/( B2*sqrt(G2) ) * ellipticF(phi,k);
            N3 = B^2/(B^2-G^2) -...
                B*G/( (B^2-G^2)*sqrt(G2) ) * ellipticE(phi,k);
            
%             phi = asin((a^2-c^2)/a^2);
%             k = sqrt((a^2 - b^2) / (a^2 - c^2))^2;
%             N1 = a*b*c / (sqrt(a^2 - c^2) * (a^2 - b^2) ) * (-ellipticE(phi,k) + ellipticF(phi,k));
%             N2 = -a*b*c^2 / ( a*b*(b^2-c^2) ) ...
%                  +a*b*c*sqrt(a^2-c^2) / ( (a^2-b^2)*(b^2-c^2) )* ellipticE(phi,k) ...
%                  -a*b*c / ( sqrt(a^2-c^2)*(a^2-b^2) ) * ellipticF(phi,k);
%             N3 =  a*b^2*c/(a*c*(b^2-c^2)) ...
%                  -a*b*c/(sqrt(a^2-c^2)*(b^2-c^2)) * ellipticE(phi,k) ;
        else 
            error('Demag error: For a general ellipsoid, dimensions must be ordered [a,b,c], where a>=b>=c')
            
        end
                        
    case 'Rectangle'
        m = @(a,b,c) sqrt(a^2 + b^2 + c^2);
        D = @(a,b,c) ...
            (b^2-c^2)/(2*b*c)*log( (m(a,b,c) - a) / ( m(a,b,c) + a ) ) + ...
            (a^2-c^2)/(2*a*c)*log( (m(a,b,c) - b) / ( m(a,b,c) + b ) ) + ...
            b/(2*c)*log( (m(a,b,0) +a )/( m(a,b,0)-a ) ) + ...
            a/(2*c)*log( (m(a,b,0) +b )/( m(a,b,0)-b ) ) + ...
            c/(2*a)*log( (m(0,b,c) -b )/( m(0,b,c)+b ) ) + ...
            c/(2*b)*log( (m(a,0,c) -a )/( m(a,0,c)+a ) ) + ...
            2*atan( (a*b) / ( c * m(a,b,c) ) )+ ...
            (a^3 + b^3 -2*c^3)/(3*a*b*c) + ...
            (a^2 + b^2 -2*c^2)/(3*a*b*c)*m(a,b,c) + ...
            c/(a*b)*( m(a,0,c) + m(0,b,c) ) - ...
            (m(a,b,0)^3  +m(0,b,c)^3 + m(a,0,c)^3)/(3*a*b*c);
        
        N3 = D(a,b,c)/(pi);
        N1 = D(b,c,a)/(pi);
        N2 = D(c,a,b)/(pi);
end

check = N1+N2+N3 -1 < 1e-6;
if ~check 
    error('Demag factors don''t add to 1!')
end

%return to original order 
N(Sort_ind(1))= N1;
N(Sort_ind(2))= N2;
N(Sort_ind(3))= N3;

end