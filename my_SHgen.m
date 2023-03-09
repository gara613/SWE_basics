%% [F,Pnm,mP,dP] = my_SHgen(m,n,s,varargin)
% Vectorial field expansion
% Series expansion of a given field in terms of Scalar spherical harmonics
%
% inputs
%   * n: (scalar) Degree 
%   * m: (scalar) Order
%   * s: (scalar) Function (1: odd <-> TE, 2: even <-> TM)
%   * varargin: 
%       - theta: vector
%       - phi: vector
%       - test: boolean, enables the comparison between the calculated
%       Functions and direct and (selected) theoretical expressions taken
%       from Hansen book p.322-324, for n == 2 | n == 3 | n == 4 
%
% outputs 
%   * F: Vectorial Spherical Harmonic (struct with fields theta and phi)
%   * Pnm: normalized associated Legendre function
%   * mP: normalized associated Legendre function, times m, over sin\theta
%   * dP: derivative of normalized associated Legendre function
%
% Germán A. Ramírez
% EPFL - MAG
% Feb 2023

function [F,Pnm,mP,dP] = my_SHgen(n,m,s,varargin)
    if exist('varargin','var') & ~isempty(varargin)        
        el = unique(varargin{1}(:));
        az = unique(varargin{2}(:));
        dx = mean(diff(el));
        test = varargin{3};
    else
        dx = 2*pi/180;     % 2 deg sampling
        el = 0:dx:pi;
        az = 0:dx:2*pi;
        test = false;
    end  

    Pnm = (-1)^m*sqrt( (2*n+1)/(2*prod(n-abs(m)+(1:2*abs(m)))) )*legendre(n,cos(el)).';

    % NOTE: singularity extraction in Legendre function divided by sin\theta,
    % mP = m/sin\theta P_n^m(cos\theta), uses Recurrence relations different to those 
    % in Stratton and Hansen due to the phsae term, see documentation on derivation
    % dP = d/d\theta P_n^m(cos\theta) is the Derivative of P_n^m(cos\theta) 

    if m == 0
        mP = zeros(size(el)); 
        dP = squeeze(Pnm(:,2)); % = Pn1
    elseif abs(m) < n
        mP = (m/abs(m))*(-0.5*cos(el).*( (n-abs(m)+1)*(n+abs(m))*Pnm(:,abs(m)) + Pnm(:,abs(m)+2) ) + abs(m)*sin(el).*Pnm(:,abs(m)+1));
        dP = -0.5*( (n-abs(m)+1)*(n+abs(m))*Pnm(:,abs(m)) - Pnm(:,abs(m)+2) );
    elseif abs(m) == n
        mP = (m/abs(m))*(-0.5*cos(el).*( (n-abs(m)+1)*(n+abs(m))*Pnm(:,abs(m)) ) + abs(m)*sin(el).*Pnm(:,abs(m)+1));
        dP = -0.5*( (n-abs(m)+1)*(n+abs(m))*Pnm(:,abs(m)) );
    end
    
    A = sqrt(2/(n*(n+1)))*(-m/abs(m))^m *exp(1i*m*az).';
    if m==0 
        A = sqrt(2/(n*(n+1)));
    end
    if s ==1        % F1 functions
        F.theta = (-1i)^(n+1)*1i*mP*A;
        F.phi = -(-1i)^(n+1)*dP*A;
    elseif s ==2	% F2 functions
        F.theta = (-1i)^(n)*dP*A;
        F.phi = (-1i)^(n)*1i*mP*A;
    end

%% Comparison to direct and (selected) theoretical expressions
	if test
        % Direct expressions
        mPdir = m*Pnm(:,abs(m)+1)./sin(el);
        dPdir = gradient(Pnm(:,abs(m)+1),dx);
        mPerr = sum(abs(mPdir(:) - mP(:)).^2)/length(mP(:));
        dPerr = sum(abs(dPdir(:) - dP(:)).^2)/length(dP(:));

        % Theoretical expressions taken from Hansen book p.322-324
        if n == 2 | n == 3 | n == 4 
            if n == 2 & abs(m) == 0
                P_th  = sqrt(10)/8*(3*cos(2*el)+1);
                mP_th = 0*el;
                dP_th = -3*sqrt(10)/4*sin(2*el);        
            elseif n == 2 & abs(m) == 1
                P_th  = sqrt(15)/4*sin(2*el);
                mP_th = (m/abs(m))*sqrt(15)/2*cos(el);
                dP_th = sqrt(15)/2*cos(2*el);
            elseif n == 2 & abs(m) == 2
                P_th  = -sqrt(15)/8*(cos(2*el)-1);
                mP_th = (m/abs(m))*sqrt(15)/2*sin(el);
                dP_th = sqrt(15)/4*sin(2*el);
            elseif n == 3 & abs(m) == 0
                P_th  = sqrt(14)/16*(5*cos(3*el)+3*cos(el));
                mP_th = 0*el;
                dP_th = -3*sqrt(14)/16*(5*sin(3*el)+sin(el));
            elseif n == 3 & abs(m) == 1
                P_th  = sqrt(42)/32*(5*sin(3*el)+sin(el));
                mP_th = (m/abs(m))*sqrt(42)/16*(5*cos(2*el)+3);
                dP_th = sqrt(42)/32*(15*cos(3*el)+cos(el));
            elseif n == 3 & abs(m) == 2
                P_th  = -sqrt(105)/16*(cos(3*el)-cos(el));
                mP_th = (m/abs(m))*sqrt(105)/4*sin(2*el);
                dP_th = sqrt(105)/16*(3*sin(3*el)-sin(el));
            elseif n == 3 & abs(m) == 3
                P_th  = -sqrt(70)/32*(sin(3*el)-3*sin(el));
                mP_th = (m/abs(m))*-3*sqrt(70)/16*(cos(2*el)-1);
                dP_th = -3*sqrt(70)/32*(cos(3*el)-cos(el));
            elseif n == 4 & abs(m) == 0
                P_th  = 3*sqrt(2)/128*(35*cos(4*el)+20*cos(2*el)+9);
                mP_th = 0*el;
                dP_th = -15*sqrt(2)/32*(7*sin(4*el)+2*sin(2*el));
            elseif n == 4 & abs(m) == 1
                P_th  = 3*sqrt(10)/64*(7*sin(4*el)+2*sin(2*el));
                mP_th = (m/abs(m))*3*sqrt(10)/32*(7*cos(3*el)+9*cos(el));
                dP_th = 3*sqrt(10)/16*(7*cos(4*el)+cos(2*el));
            elseif n == 4 & abs(m) == 2
                P_th  = -3*sqrt(5)/64*(7*cos(4*el)-4*cos(2*el)-3);
                mP_th = (m/abs(m))*3*sqrt(5)/16*(7*sin(3*el)+3*sin(el));
                dP_th = 3*sqrt(5)/16*(7*sin(4*el)-2*sin(2*el));
            elseif n == 4 & abs(m) == 3
                P_th  = -3*sqrt(70)/64*(sin(4*el)-2*sin(2*el));
                mP_th = (m/abs(m))*-9*sqrt(70)/32*(cos(3*el)-cos(el));
                dP_th = -3*sqrt(70)/16*(cos(4*el)-cos(2*el));
            elseif n == 4 & abs(m) == 4
                P_th  = 3*sqrt(35)/128*(cos(4*el)-4*cos(2*el)+3);
                mP_th = (m/abs(m))*-3*sqrt(35)/16*(sin(3*el)-3*sin(el));
                dP_th = -3*sqrt(35)/32*(sin(4*el)-2*sin(2*el));
            end
%%
            figure, 
            subplot(2,2,[1,2]), plot(el,squeeze(Pnm(abs(m)+1,:,1)),	el,P_th,'--', 'linewidth',2);  legend('def','th'); title(['Pmn (m=',num2str(m),'n=',num2str(n),')']);
            subplot 223, plot(el,mP(:,1),	el,mPdir(:,1),'--',     el,mP_th,'-.','linewidth',2);  legend('rec','dir','th'), title({['mP (m=', num2str(m),',n=',num2str(n),')'], ['err=',num2str(mPerr)]});
            subplot 224, plot(el,dP(:,1),   el,squeeze(dPdir(1,:,1)),'--',	el,dP_th,'-.','linewidth',2);  legend('rec','dir','th'); title({['dP (m=', num2str(m),',n=',num2str(n),')'], ['err=',num2str(dPerr)]});
       
        else            
            figure,
            subplot 121, plot(el,mP(:,1), el,mPdir(:,1), '--','linewidth',2); legend('rec','dir'); title({['mP (m=', num2str(m),',n=',num2str(n),')'], ['err=',num2str(mPerr)]});
            subplot 122, plot(el,dP(:,1), el,squeeze(dPdir(1,:,1)), '--','linewidth',2); legend('rec','dir'); title({['dP (m=', num2str(m),',n=',num2str(n),')'], ['err=',num2str(dPerr)]});
        end
	end
end