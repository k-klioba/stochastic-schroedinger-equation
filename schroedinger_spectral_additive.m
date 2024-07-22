%% numerical simulation of stochastic Schrödinger equation
%% authors: Katharina Klioba, Hamburg University of Technology,
%% Mark Veraar, Delft University of Technology,
%% for the paper 'Pathwise uniform convergence of 
%% Time Discretisation Schemes for SPDEs'
%% 12.07.2024
%% 
%% space: spectral Galerkin, periodic BC
%% time: exponential Euler, implicit Euler, Crank-Nicolson 
%% initial condition: Gaussian pulse centred at pi
%% additive noise g(t)=-iQ^(1/2), Q has eigenvalues 1/(1+n^5.1)
%% expected rate 1|2/3|1/2 for pathwise uniform errors for EXP|CN|IE

clear all; close all;
rng(42)

fprintf('schroedinger.m\n'); 

% kind of spatial error: measured in H^{sigma}
sigma = 0; % 0 for L2 (default), 1 for H^1
barsigma = 1 - (1 - 1/2) * (sigma == 0); % prefactor 1/2 for L2, 1 else

% number of samples
N_samples = 100;

% space discretisation
Il = 0; Ir = 2*pi;
M = 2^10; % number of Fourier modes 

% wave numbers
kk = [0:M/2 -M/2+1:-1]';
ksquared = kk.^2;

% initial value in frequency space
Fu_0 = 1./( ones(size(kk)) + abs(kk).^6 );

% eigenvalues for the noise
regularity_exponent = 5.1; % default: 5.1 for sigma = 0
Lambda = (1 ./ ( ones(M,1) + (abs(kk)).^regularity_exponent ));
sqrtLambda = (1 ./ sqrt( ones(M,1) + (abs(kk)).^regularity_exponent ));

% time discretisation
t_0 = 0; t_end = 0.5;
dt_num = [2.^(-5:-1:-9),2^(-12)]; % vector of different time steps
N_dt = length(dt_num);
dt_exact = dt_num(end); % time step for 'exact' solution
Nt_num = floor((t_end - t_0)./ dt_num);
Nt_exact = floor((t_end - t_0) / dt_exact);

% exact solution via exponential Euler with small time step
FU_exact = zeros(N_samples,M,Nt_exact);

% errors
spaceErrEXP = zeros(N_dt-1,N_samples,Nt_exact);
spaceErrIE = zeros(N_dt-1,N_samples,Nt_exact);
spaceErrCN = zeros(N_dt-1,N_samples,Nt_exact);
ErrEXP = zeros(N_dt-1,1);
ErrIE = zeros(N_dt-1,1);
ErrCN = zeros(N_dt-1,1);


% simulate pathwise
for m = 1:N_samples
 
  %% noise
  Wtt = sqrt(dt_exact) * randn(M,Nt_exact);

  %% calculate exact solution
  FU1_exact = Fu_0; % EXP for exact solution
  for t = 1:Nt_exact
      FU1_exact = exp(1i * dt_exact * ksquared) .* (FU1_exact - 1i * sqrtLambda .* Wtt(:,t) );
      FU_exact(m,:,t) = FU1_exact; 
  end  
     
  %% approximate solutions for different time step sizes
  Wn = zeros(M,1);  
  for j = 1:N_dt-1 % for each time step size dt except the exact one

    Nt = (t_end-t_0) / dt_num(j);
    dt = dt_num(j);
    ratio = dt / dt_exact; % multiple of exact time step
    
    %% initialise
    FU1EXP = Fu_0; 
    FU1CN = Fu_0; 
    FU1IE = Fu_0;
   
    %% iteration: approximate solutions
    for t = 1:Nt

      % sum Brownian motion over corresponding time interval ((t-1)*dt,t*dt]
      Wtn = sum(Wtt(:,ratio * (t-1) + 1 : ratio * t),2);

      %% Exponential Euler 
      FU1EXP = exp(1i * dt * ksquared) .* (FU1EXP - 1i * sqrtLambda .* Wtn );
      %% Implicit Euler
      FU1IE = (ones(M,1) - 1i * dt * ksquared).^(-1) .* (FU1IE - 1i * sqrtLambda .* Wtn); 
      %% Crank-Nicolson
      FU1CN = (ones(M,1) - 1i * dt/2 * ksquared).^(-1) .* (ones(M,1) + 1i * dt/2 * ksquared) .* (FU1CN - 1i * sqrtLambda .* Wtn);  
     
      % calculate spatial H^{sigma} norm of error vector at this time step
      index_t_exact = t * ratio;
      spaceErrEXP(j,m,t) = norm( barsigma * ( abs(kk).^(2*sigma) + ones(size(kk)) ) .* (FU1EXP - FU_exact(m,:,index_t_exact).') , 2)^2;
      spaceErrIE(j,m,t) = norm( barsigma * ( abs(kk).^(2*sigma) + ones(size(kk)) ) .* (FU1IE - FU_exact(m,:,index_t_exact).') , 2)^2;
      spaceErrCN(j,m,t) = norm( barsigma * ( abs(kk).^(2*sigma) + ones(size(kk)) ) .* (FU1CN - FU_exact(m,:,index_t_exact).') , 2)^2;

    end % t
  end  % j   
  m
end % m 

% pathwise uniform errors
ErrEXP = sqrt(mean( max(spaceErrEXP,[],3), 2 ));
ErrIE = sqrt(mean( max(spaceErrIE,[],3), 2 ));
ErrCN = sqrt(mean( max(spaceErrCN,[],3), 2 ));

set(0,'DefaultTextFontSize',12)
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultLineMarkerSize',10)

figure(),
loglog(dt_num(1:end-1),dt_num(1:end-1).^(.5) * ErrEXP(1)/dt_num(1)^(.5)*1.2,'k--', ...
       dt_num(1:end-1),dt_num(1:end-1).^1 * ErrEXP(1)/dt_num(1)*0.9,'k-.', ...
       dt_num(1:end-1),ErrEXP,'bs-', ... 
       dt_num(1:end-1),ErrIE,'md-', ...
       dt_num(1:end-1),ErrCN,'c*-');
xlabel('$k$','Interpreter','latex','FontSize',14);
ylabel('Error','Rotation',90,'HorizontalAlignment','right','FontSize',14);
legend('Slope 1/2','Slope 1','EXP', ... 
     'IE','CN','Location','SouthEast');

% numerical convergence rates
numericalConvRateEXP = log2(ErrEXP(1:end-1)).' - log2(ErrEXP(2:end)).';
numericalConvRateIE = log2(ErrIE(1:end-1)).' - log2(ErrIE(2:end)).';
numericalConvRateCN = log2(ErrCN(1:end-1)).' - log2(ErrCN(2:end)).';

averageNumericalConvRateExponentialEuler = 1/(N_dt-2) * ( log2(ErrEXP(1)) - log2(ErrEXP(end)) )
averageNumericalConvRateImplicitEuler = 1/(N_dt-2) * ( log2(ErrIE(1)) - log2(ErrIE(end)) )
averageNumericalConvRateCrankNicolson = 1/(N_dt-2) * ( log2(ErrCN(1)) - log2(ErrCN(end)) )

print -depsc2 -r0 schroedinger.eps

save('schroedinger.mat','dt_num','ErrEXP', 'ErrIE', ...
'ErrCN')
