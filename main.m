%MAIN Runs scripts mainoutdist and mainoutnoise to simulate the scenarios
%described in Gleizer and Mazo Jr. (2019). The results in output figures 1
%and 2 should match those in Figs. 4 and 5 of the paper, respectively.
%
%   Requirements: Optimization Toolbox, Ellipsoidal Toolbox 
%   The Ellipsoidal Toolbox can be found in
%       https://nl.mathworks.com/matlabcentral/fileexchange/21936-ellipsoidal-toolbox-et
% 
%   Author: Gabriel de A. Gleizer, 2018--2019 (g.gleizer@tudelft.nl)
%
% REFERENCES:
% G. A. Gleizer and M. Mazo Jr. Self-triggered output-feedback control of
%     LTI systems subject to disturbances and noise. Submitted to
%     Automatica, 2019.
% G. A. Gleizer and M. Mazo Jr. Self-triggered output-feedback control for
%     perturbed linear systems. IFAC-PapersOnLine, 51(23):248--253
%     (Presented in IFAC NecSys 2018).

%% Get NecSys data for comparison
clearvars;

% Run NecSys script
mainoutdist;

% Rename log variables
% 1. STC
xlognecsys = xlog;
klognecsys = klog;
kslognecsys = kslog;
dklognecsys = dklog;
% 2. PETC
xlogpetc = xloge;
klogpetc = kloge;
kslogpetc = ksloge;
dklogpetc = dkloge;

% Put states back in original form (it is the canonical form for the
% noiseless case)
xlognecsys(:,1:np) = xlognecsys(:,1:np)*(inv(T)');
xlogpetc(:,1:np) = xlogpetc(:,1:np)*(inv(T)');

% The log of instants for the petc is shifted because kappa is computed a
% posteriori. Shift back for comparison purposes.
kslogpetc = kslogpetc - dklogpetc;

clearvars -except *necsys *petc

%% Now run the Automatica version without noise

% Change verbosity of the Ellipsoidal Toolbox.
global ellOptions; 
intersection_ea(ell_unitball(1),ell_unitball(1));  % Dummy command to 
                                                   % initialize it.
if isstruct(ellOptions)
    ellOptions.verbose = 0;
end

% Set noise to zero for the first experiment.
V_GLOBAL = 0;

% Run script
mainnoise;
% Put states back to original coordinates
xlog(:,1:np) = xlog(:,1:np)*(inv(T)');

%% Plots (Figure 4 in the paper)
try
    close(1);
end

LINEWIDTH = 1;
MARKERSIZE = 12^2;

figure(1);  
hdle(1) = subplot(211);
plot(klog*h,sqrt(sum(xlog.^2,2)),'-','LineWidth',LINEWIDTH); hold all;
plot(klognecsys*h,sqrt(sum(xlognecsys.^2,2)),'--','LineWidth',LINEWIDTH); 
plot(klogpetc*h,sqrt(sum(xlogpetc.^2,2)),':','LineWidth',LINEWIDTH)
legend('PSTC', 'STC from \cite{gleizer2018selftriggered}', 'PETC');
ylabel('State norm $\|\xi\|$');
grid on;

hdle(2) = subplot(212);
scatter(kslog*h,dklog,MARKERSIZE,'o'); hold all;
scatter(kslognecsys*h,dklognecsys,MARKERSIZE,'.');
scatter(kslogpetc*h,dklogpetc,MARKERSIZE,'x');
legend('PSTC', 'STC from \cite{gleizer2018selftriggered}', 'PETC')
ylabel('$\Delta k$ (samples)');
ylim([0 kfinal]);
grid on;

linkaxes(hdle,'x');
grid on;

%% Now run the PSTC with noise
V_GLOBAL = 0.01;

% Experiment one: epsilon=0
mainnoise;
% Rename log variables to avoid rewriting from the second experiment
xlog0 = xlog;
klog0 = klog;
kslog0 = kslog;
dklog0 = dklog;
dkzlog0 = dkzlog;

% Experiment two: epsilon^2=0.01
TRIG_LEVEL_GLOBAL = 0.01;
mainnoise;

%% Plots
try
    close(2);
end

LINEWIDTH = 1;
MARKERSIZE = 12^2;

figure(2);  
hdle(1) = subplot(211);
scatter(kslog0*h,dklog0,MARKERSIZE,'o'); hold all;
plot(kslog0*h,dkzlog0,'LineWidth',LINEWIDTH);
legend('Our STC', 'PETC')
ylabel('$\kappa$ (samples)');
ylim([0 kfinal]);
grid on;

hdle(2) = subplot(212);
scatter(kslog*h,dklog,MARKERSIZE,'o'); hold all;
plot(kslog*h,dkzlog,'LineWidth',LINEWIDTH);
legend('Our STC', 'PETC')
ylabel('$\kappa$ (samples)');
ylim([0 kfinal]);
grid on;

linkaxes(hdle,'x');
%ylim(yl);
grid on;
   
%% Get some other metrics

% state norm vectors (xn0 for epsilon=0, the other for epsilon=0.1)
xn0 = sqrt(sum(xlog0.^2,2));
xn = sqrt(sum(xlog.^2,2));

% Get data after 6 time units (transients gone)
i0after6 = kslog0 > 6/h;
iafter6 = kslog > 6/h;
meanxn0 = mean(xn0(i0after6));
meanxn = mean(xn(iafter6));

fprintf('Average state norms after t=6: %g (eta > 0) vs %g (eta > %g)\n',...
    meanxn0, meanxn, TRIG_LEVEL_GLOBAL);

% How more frequently does PSTC trigger when compared to ETC? Get values
% after state norm reaches 10% of initial norm and 1% of initial norm 
i10 = xn > xn(1)/10;
i100 = xn > xn(1)/100;
k10 = klog(find(i10,1,'last'));
k100 = klog(find(i100,1,'last'));
i10 = kslog <= k10;
i100 = kslog <= k100;
% "Worseness" metrics
worseness10 = mean(dkzlog(i10))/mean(dklog(i10)) - 1;
worseness100 = mean(dkzlog(i100))/mean(dklog(i100)) - 1;


%% For the reviewer
figure(3); 
plot(klog0*h,sqrt(sum(xn0.^2,2)),'-','LineWidth',LINEWIDTH); hold all;
plot(klog*h,sqrt(sum(xn.^2,2)),'--','LineWidth',LINEWIDTH);
legend('PSTC $\epsilon=0$', 'PSTC $\epsilon = 0.1$');
ylabel('State norm $\|\xiv\|$');
xlabel('$t$')
grid on;

matlab2tikz('with_noise.tex','height','30mm','width','75mm',...
    'interpretTickLabelsAsTex',true,'standalone',false,...
    'parseStrings',false,'showInfo', false);