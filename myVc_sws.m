%% F = myVc_sws(q,varargin)
% Vectorial farfield reconstruction from a Spherical wave expansion
% Synthesis, given the coefficients of a series expansion in spherical harmonics
% inputs
%   * q: Vector of coefficients
%   * varargin: theta, phi
% 
% Germán A. Ramírez
% EPFL - MAG
% Sept 2022

function F = myVc_sws(q,varargin)
    if exist('varargin','var') & ~isempty(varargin)        
        el = unique(varargin{1}(:));
        az = unique(varargin{2}(:));
        
        [phi,theta] = meshgrid(az,el);
    else
        el = linspace(0,pi,91).';
        az = linspace(0,2*pi,181).';
        [phi, theta] = meshgrid(az,el);
    end   
    
    F.theta = zeros(size(theta));
    F.phi = zeros(size(theta));
    J = length(q);              % Number of harmonics
        
    N = -1 + sqrt(J/2+1);       % J = 2N(N+2)
    j = 1;
    
    for n = 1:N         % DEGREE
        Pnm = legendre(n,cos(el)).';   % legendre(n,X) computes the associated Legendre functions of degree n and order m = 0, 1, ..., n, for each x in X

        for m = -n:n    % ORDER
            A = (-1)^m*sqrt((2*n+1)/(2*prod(n-abs(m)+(1:2*abs(m)))) ) * sqrt(2/(n*(n+1)))*(-m/abs(m))^m *exp(1i*m*az).'; 
                
            if m == 0
            	A = sqrt((2*n+1)/2) * sqrt(2/(n*(n+1))); 
                mP = zeros(size(el));
                dP = Pnm(:,2); % Pn1

            elseif abs(m)<n    % Recursive relations for the Legendre functions, derived following Stratton and Hansen 
            	mP = (m/abs(m))*( -0.5*cos(el).*( (n-abs(m)+1)*(n+abs(m))*Pnm(:,abs(m)) + Pnm(:,abs(m)+2) ) + abs(m)*sin(el).*Pnm(:,abs(m)+1) ); 
                dP = -0.5*( (n-abs(m)+1)*(n+abs(m))*Pnm(:,abs(m)) - Pnm(:,abs(m)+2) );
                    
            elseif abs(m)==n
            	mP = (m/abs(m))*( -0.5*cos(el).*( (n-abs(m)+1)*(n+abs(m))*Pnm(:,abs(m)) ) + abs(m)*sin(el).*Pnm(:,abs(m)+1) );
                dP = -0.5*( (n-abs(m)+1)*(n+abs(m))*Pnm(:,abs(m)) );
            end
                          
            F1.theta = (-1i)^(n+1)*1i*mP*A;
            F1.phi = -(-1i)^(n+1)*dP*A;                
            F2.theta = (-1i)^(n)*dP*A;
            F2.phi = (-1i)^(n)*1i*mP*A;
            
            F.theta = F.theta + q(j)*F1.theta + q(j+1)*F2.theta;    
            F.phi = F.phi + q(j)*F1.phi + q(j+1)*F2.phi;
            
            j=j+2;
        end    
    end       
end