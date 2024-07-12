%% numerical simulation of stochastic Schrödinger equation
%% authors: Katharina Klioba, Hamburg University of Technology,
%% Mark Veraar, Delft University of Technology,
%% for the paper 'Pathwise uniform convergence of 
%% Time Discretisation Schemes for SPDEs'
%% 12.07.2024
%% 
%% Create numerical convergence rate plots
%% requires 'ErrEXP_additive_samples100_M1024_dt25912.mat'
%% as well as 'ErrEXP_additive_samples100_M1024_dt25912.mat'

load('ErrEXP_additive_samples100_M1024_dt25912.mat', 'ErrEXP', 'ErrIE', 'ErrCN')
ErrEXPadd = ErrEXP;
ErrIEadd = ErrIE;
ErrCNadd = ErrCN;

load('ErrEXP_multiplicative_samples100_M1024_dt25912.mat', 'ErrEXP', 'ErrIE', 'ErrCN')
ErrEXPmult = ErrEXP;
ErrIEmult = ErrIE;
ErrCNmult = ErrCN;


set(0,'DefaultTextFontSize',12)
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultLineMarkerSize',10)

figure;

% additive noise
subplot(1,2,1);
loglog(dt_num(1:end-1),dt_num(1:end-1).^(.5) * ErrIEadd(1)/dt_num(1)^(.5)*1.2,'k--', ...
       dt_num(1:end-1),dt_num(1:end-1).^1 * ErrEXPadd(1)/dt_num(1)*0.8,'k-.', ...
       dt_num(1:end-1),ErrEXPadd,'bs-', ... 
       dt_num(1:end-1),ErrIEadd,'md-', ...
       dt_num(1:end-1),ErrCNadd,'c*-');
xlabel('$k$','Interpreter','latex','FontSize',14);
ylabel('Error','Rotation',90,'HorizontalAlignment','right','FontSize',14);
legend('Slope 1/2','Slope 1','EXP', ... 
     'IE','CN','Location','SouthEast');
axis([(dt_num(end-1)) (dt_num(1)) 10^-3 10^-1])

% multiplicative noise
subplot(1,2,2);
loglog(dt_num(1:end-1),dt_num(1:end-1).^(.25) * ErrIEmult(1)/dt_num(1)^(.25)*1.2,'k:', ...
       dt_num(1:end-1),dt_num(1:end-1).^(.5) * ErrEXPmult(1)/dt_num(1)^(.5)*0.7,'k--', ...
       dt_num(1:end-1),ErrEXPmult,'bs-', ... 
       dt_num(1:end-1),ErrIEmult,'md-', ...
       dt_num(1:end-1),ErrCNmult,'c*-');
xlabel('$k$','Interpreter','latex','FontSize',14);
ylabel('Error','Rotation',90,'HorizontalAlignment','right','FontSize',14);
legend('Slope 1/4','Slope1/2','EXP', ... 
     'IE','CN','Location','SouthEast');
axis([(dt_num(end-1)) (dt_num(1)) 10^-2.5 10^-0.5])