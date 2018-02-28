%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%              Magnon's frequency conversion to wavevector's              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matlab script to calculate the wavevector's value for each magnon's
% frequency. The calculation is based on Vladimir V. Vlasov calculation.

clearvars,
% close all,

%% Constant
mu0=4*pi*1e-7;  % Vacuum permeability (V.s/A.m)
h_bar=1.054e-34;    % reduced Planck constant or Dirac constant (J.s)
gamma=1.76e11;    % gyromagnetic ratio (rad-1.s-1.T-1 ?)

HAll = 6.5425/mu0; % DC magnetic field (T), double
thickness = 30e-9; % thickness of the film (m), double
xi = (45*pi/2)/90; % DC magnetic field out of plane angle in rad
L = 30e-9; % thickness of the thin film (m)

D2=430; % exchange constant (meV*A^2) %550 for cobalt
D = D2*1.602e-42/(gamma*mu0*h_bar); % exchange constant for fields

M0=484.1e3; % saturation magnetization %1.8 for cobalt

load('Data.mat','tet');    % equilibrum magnetization out-of-plane angle in rad
load('Data2.mat','fqcy','spectraY');
% load('Data2.mat','fqcy');


%% Function
% definition of spinwaves resonnant frequencies (in Hz)
fun1 = @(k) sqrt((mu0*gamma)^2 * ((HAll * sin(xi) + (D * k^2 - M0) * sin(tet))^2 + (HAll * cos(xi) + (D * k^2 + M0) * cos(tet)) * (HAll * cos(xi) + D * k^2 * sin(tet))))/(2*pi);



%% Determination of magnon's k
Nk=length(fqcy); % number of points for dispersion relation calculations
n_max = 20;
kmax=pi*n_max/L;    % maximal value of the wave vector
k_ini = linspace(0,kmax,Nk);
k_est = zeros(1,length(fqcy));
opt = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off','TolFun',1e-6,'TolX',1e-6);   % selection of options

for nn = 2:length(fqcy)
% for nn = 2    
    if fqcy(nn) < fun1(0)
        k_est(nn) = 0;

    else
        fun = @(k) (fqcy(nn)*(2*pi))^2 - ((mu0*gamma)^2 * ((HAll * sin(xi) + (D * k^2 - M0) * sin(tet))^2 + (HAll * cos(xi) + (D * k^2 + M0) * cos(tet)) * (HAll * cos(xi) + D * k^2 * sin(tet))));
        k_est(nn) = fsolve(fun,k_ini(nn),opt);
    end
end

%% Figure

figure(1)
% semilogx(spectraY,k_est./1e8)
semilogy(k_est,spectraY,'LineWidth',2)
xlim([0 13e8])
ax1 = gca;
fig1 = gcf;
% set(ax1,'XTickMode','Manual','XTick',[0 0.5 1],'XTickLabel',[0 0.5 1])%x,'XDir','reverse'
set(ax1,'XTickMode','Manual','XTick',[0 2.5e8 5e8 7.5e8 10e8 12.5e8],'XTickLabel',[0 2.5 5 7.5 10 12.5])
% set(ax1,'YTickMode','Manual','YTick',[0 200 400 600 800 1000],'YTickLabel',[0 0.2 0.4 0.6 0.8 1]');
set(ax1,'YTickMode','Manual','YTick',[10^(-4) 10^(-2)],'YTickLabel',[{'10^{-4}' '10^{-2}'}]);
set(ax1,'FontSize',25,'FontName','CMU Serif','FontWeight','bold');
xlabel('Wavevector (10^{-1} nm^{-1})','FontName','CMU Serif','FontWeight','bold')
ylabel('Amplitude (A.U.)','FontName','CMU Serif','FontWeight','bold')
set(fig1,'units','centimeters','outerposition',[0 0 30 15])
filename = sprintf('spectrumAsAFunctionOfWavevector');
print(strcat(filename,'.svg'),'-dsvg');
print(strcat(filename,'.png'),'-dpng');