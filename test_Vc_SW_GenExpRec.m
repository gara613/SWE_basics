% Testing of the my_SHgen, myVc_swe, myVc_sws functions
% Generates a linear combination of the specified vector spherical
% harmonics, then calculates the expansion coefficients and the
% corresponding function reconstruction. 
% Relative RMS error between original and reconstructed functions is displayed as well as the and computation time

close all; 
clc; clearvars; 

dx = 2*pi/180;     % 2 deg sampling
el = 0:dx:pi;
az = 0:dx:2*pi;
[phi,theta] = meshgrid(az,el);
[N_Theta, N_Phi] = size(theta);
L = N_Theta*N_Phi;

%% ---------------- Pnm: RECURSIVE VS DIRECT VS THEORETICAL ----------------
n = [1, 2, 4, 5, 2, 4, 3, 3, 4, 5, 3];  
m = [0, -1, 3, -4, 0, 3, -2, 2, -1, 2, -2];
s = [1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 2]; 
a = [1-1i, 2+1i, 1+0.5i, 0.5-0.5i, 1+1i, 1-1i, 1.5-0.5i, 2-1i, 1+0.5i, 2+1i, 0.5+2i];
%  L = length(n);
%  a = (1-rand(1,L)).*exp(1i*pi*(2-rand(1,L)));

j = n.*(n+1) + m,    % NOTE: Actual positions within the q vector: j = 2*(n*(n+1) +m-1) + [1:2]


%% Target function
%[F,Pnm,mP,dP] = my_SHgen(n(1),m(1),s(1),el,az,false);
F.theta = zeros(size(theta));
F.phi = zeros(size(theta));
for cont = 1:length(n) % cont = 2:length(n)
    tF = my_SHgen(n(cont),m(cont),s(cont),el,az,false);
    F.theta = F.theta + a(cont)*tF.theta;
    F.phi = F.phi+ a(cont)*tF.phi;
end

%% Decomposition 
N = 10;  
J = 2*N*(N+2);
Jr = J;                         % Number of coefficients for reconstruction 

tic
    q_F = myVc_swe(F,J,theta,phi);  
toc

figure,
jj = 1:J/2;
subplot 211
stem(jj.', [real(q_F(1:2:J)), imag(q_F(1:2:J))], 'linewidth',2); grid on; 
title('Odd Coeffs'), legend('real','imag'); axis([1,50,-2,2]),
subplot 212
stem(jj.',[real(q_F(2:2:J)), imag(q_F(2:2:J))], 'linewidth',2); grid on; 
title('Even Coeffs'), legend('real','imag'); axis([1,50,-2,2]),


%% Reconstruction 
tic
    F_rec = myVc_sws(q_F(1:Jr),theta,phi);  
toc

%% Validation 
targ_F = sqrt(abs(F.theta).^2 + abs(F.phi).^2);
targ_F_rec = sqrt(abs(F_rec.theta).^2 + abs(F_rec.phi).^2);

plotPattern(targ_F,theta,phi,'norma',false,'scale','dB','plStyle','pwr','coords','pol3D','limsScale',[-30,30],'cut2Plot',[-1,-1],'multPatts',false,'aggrPatts',false)
title({'Normalized Power Pattern (dB)', 'Original'});
plotPattern(targ_F_rec,theta,phi,'norma',false,'scale','dB','plStyle','pwr','coords','pol3D','limsScale',[-30,30],'cut2Plot',[-1,-1],'multPatts',false,'aggrPatts',false)
title({'Normalized Power Pattern (dB)', 'Reconstructed'});

erFd_th = F.theta - F_rec.theta;
erFd_ph = F.phi - F_rec.phi;
erFd = targ_F - targ_F_rec;

rel_rmse_th = sqrt( 1/L* sum(sum( abs(erFd_th).^2 )) )/(max(abs(F.theta(:)))-min(abs(F.theta(:)))),     
rel_rmse_ph = sqrt( 1/L* sum(sum( abs(erFd_ph).^2 )) )/(max(abs(F.phi(:)))-min(abs(F.phi(:)))),     
rel_rmse_tot = sqrt( 1/L* sum(sum( abs(erFd).^2 )) )/(max(targ_F(:))-min(targ_F(:))),               % Normalized relative error"