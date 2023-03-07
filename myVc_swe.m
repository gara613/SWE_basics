%% q = myVc_swe(F,J,varargin)
% Vectorial field expansion
% Series expansion of a given field in terms of Vector spherical harmonics
%
% Inputs
%   * F: Vector Field to expand. (Struct with fields theta and phi)
%   * J: Number of harmonics, J = 2N(N+2)
%   * varargin: theta, phi (rad)
%
% Outputs
%   * q: Vector with odd (s=1) and even (s=2) expansion coefficients in alternate fashion
% 
% Germán A. Ramírez
% EPFL - MAG
% Sept 2022

function q = myVc_swe(F,J,varargin)
    if exist('varargin','var') & ~isempty(varargin)        
        el = unique(varargin{1}(:));
        az = unique(varargin{2}(:));        
        [phi,theta] = meshgrid(az,el);
    else
        el = linspace(0,pi,size(F.theta,1)).';
        az = linspace(0,2*pi,size(F.theta,2)).';
        [phi, theta] = meshgrid(az,el);
    end
       
    N = -1 + sqrt(J/2+1);  
    q = zeros(J,1);
    j = 0;    

    cont_j = 1;
	for n = 1:N                         % DEGREE, n = floor(sqrt((j-s)/2+1))
         j = j(end) + (1:(2*n+1)*2);    % Indices of j for a given 'n'

        Pnm = legendre(n,cos(el)).'; % legendre(n,X) computes the associated Legendre functions of degree n and order m = 0, 1, ..., n, for each x in X
        
        qj_o = zeros(1,2*n+1);
        qj_e = zeros(1,2*n+1);
        cont = 1;
        for m = -n:n    % ORDER, m = (j-s)/2+1-n*(n+1)
            % Normalization must include the (-1)^m because of the definition of Legendre function in Matlab uses the Condon-Shortley phase
            A = (-1)^m*sqrt((2*n+1)/(2*prod(n-abs(m)+(1:2*abs(m)))) ) * sqrt(2/(n*(n+1)))*(-m/abs(m))^m *exp(1i*m*az).'; 
            if m == 0
            	A = sqrt((2*n+1)/(n*(n+1)));                 
                mP = zeros(length(el),1); 
                dP = Pnm(:,2); % Pn1

            elseif abs(m)<n    % Recurrence relations differ from those in Stratton and Hansen due to the phase term, see derivation on documentation 
                mP = (m/abs(m))*(-0.5*cos(el).*( (n-abs(m)+1)*(n+abs(m))*Pnm(:,abs(m)) + Pnm(:,abs(m)+2) ) + abs(m)*sin(el).*Pnm(:,abs(m)+1));                 
                dP = -0.5*( (n-abs(m)+1)*(n+abs(m))*Pnm(:,abs(m)) - Pnm(:,abs(m)+2) );
                
            elseif abs(m)==n
                mP = (m/abs(m))*(-0.5*cos(el).*( (n-abs(m)+1)*(n+abs(m))*Pnm(:,abs(m))) + abs(m)*sin(el).*Pnm(:,abs(m)+1));
                dP = -0.5*( (n-abs(m)+1)*(n+abs(m))*Pnm(:,abs(m)) );                
            end
 
            F1.theta = (-1i)^(n+1)*1i*mP*A;    % Outer product is used as theta and phi variations are independent 
            F1.phi = -(-1i)^(n+1)*dP*A;
            F2.theta = (-1i)^(n)*dP*A;
            F2.phi = (-1i)^(n)*1i*mP*A; 

            In.odd = (F.theta.*conj(F1.theta) + F.phi.*conj(F1.phi)).*sin(theta); 
            In.even = (F.theta.*conj(F2.theta) + F.phi.*conj(F2.phi)).*sin(theta); 
            
            qj_o(cont) = 1/(4*pi)*trapz(az,trapz(el,In.odd,1));  % nested integrations on an array of numeric data -> (Least squares approach should also work)
            qj_e(cont) = 1/(4*pi)*trapz(az,trapz(el,In.even,1)); % NOTE: Advanced algorithms such as Lebedev quadrature should be used here 

            cont = cont+1;
            cont_j = cont_j+1;
        end         
        q(j(1:2:end)) = qj_o; 
        q(j(2:2:end)) = qj_e;        
	end
end