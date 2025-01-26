clear all
close all
clc
%% EXERCISE 1
tab_NaN = readtable (" OR.PA.csv ");
tab = rmmissing ( tab_NaN );
P = tab. Close ; % close price
log_p = log(P); % log prices
y = 100* diff ( log_p ); % logarithmic returns series
n = length (y);
dates = tab. Date (2: end);
vt = 1:n; % time index , t
% plot of returns
figure (1) ;
plot (dates ,y);
ylabel ("L' O r a l logarithmic returns ");
xlabel (" years ");
T = 3500; % size of the rolling window
m = n - T - 1; % observations
% Gaussian GARCH
Mdl1 = garch (1 ,1);
mPar1 = NaN (3,m);
vh_frcst1 = NaN(m ,1);
% Gaussian EGARCH
Mdl2 = egarch (1 ,1);
mPar2 = NaN (4,m);
vh_frcst2 = NaN(m ,1);
% Student -t GARCH
Mdl3 = garch (1 ,1); Mdl3 . Distribution = "t";
mPar3 = NaN (4,m);
vh_frcst3 = NaN(m ,1);
% Student -t GJR
Mdl4 = gjr (1 ,1); Mdl4 . Distribution = "t";
mPar4 = NaN (5,m);
vh_frcst4 = NaN(m ,1);
9;
% Riskmetrics
MdlES = arima ('Constant ', 0, 'D', 1,'MA ' ,-.94, 'Variance ', 1);
vh_frcst_es = NaN(m ,1);
vVar1 = [];
vVar2 = [];
vVar3 = [];
vVar4 = [];
vVarES = [];
dp = 0.01;
dt = norminv (0.01);
for i = 0:m -1
yr = y(i +1:T+i); % rolling sample
if any(i ==0:20:m -1)
[ EstMdl1 ,~ ,~ ,~] = estimate (Mdl1 , yr); % garch (1,1) gaussian
[ EstMdl2 ,~ ,~ ,~] = estimate (Mdl2 , yr); % egarch
gaussian
[ EstMdl3 ,~ ,~ ,~] = estimate (Mdl3 , yr); % garch (1 ,1) t
[ EstMdl4 ,~ ,~ ,~] = estimate (Mdl4 , yr); % gjr garch t
end
% Gaussian GARCH
mPar1 (:,i+1) = [ EstMdl1 . Constant ; EstMdl1 . ARCH
{1}; EstMdl1 . GARCH {1}];
dh_frcst1 = forecast ( EstMdl1 , 1,'Y0 ',yr);
vh_frcst1 (i+1) = dh_frcst1 ; % forecast about volatility
dmu1 = EstMdl1 . Offset ;
vVar1 = - ( dmu1 + dt* sqrt ( vh_frcst1 ));
% Gaussian EGARCH
mPar2 (:,i+1) = [ EstMdl2 . Constant ; EstMdl2 . ARCH
{1}; EstMdl2 . GARCH {1}; EstMdl2 . Leverage
{1}];
dh_frcst2 = forecast ( EstMdl2 , 1,'Y0 ',yr);
vh_frcst2 (i+1) = dh_frcst2 ;
dmu2 = EstMdl2 . Offset ;
vVar2 = - ( dmu2 + dt* sqrt ( vh_frcst2 ));
% Student -t GARCH
mPar3 (:,i+1) = [ EstMdl3 . Constant ; EstMdl3 . ARCH {1}; EstMdl3 . GARCH {1}; EstMdl3 . Distribution10 .DoF ];
dh_frcst3 = forecast ( EstMdl3 , 1,'Y0 ',yr);
vh_frcst3 (i+1) = dh_frcst3 ;
dmu3 = EstMdl3 . Offset ;
dnu3 = EstMdl3 . Distribution .DoF ;
dt_3 = tinv (dp , dnu3 )/ sqrt ( dnu3 /( dnu3 - 2));
vVar3 = - ( dmu3 + dt_3 * sqrt ( vh_frcst3 ));
% Student -t GJR
mPar4 (:,i+1) = [ EstMdl4 . Constant ; EstMdl4 . ARCH
{1}; EstMdl4 . GARCH {1}; EstMdl4 . Leverage
{1}; EstMdl4 . Distribution .DoF ];
dh_frcst4 = forecast ( EstMdl4 , 1,'Y0 ',yr);
vh_frcst4 (i+1) = dh_frcst4 ;
dmu4 = EstMdl4 . Offset ;
dnu4 = EstMdl4 . Distribution .DoF ;
dt_4 = tinv (dp , dnu4 )/ sqrt ( dnu4 /( dnu4 - 2));
vVar4 = - ( dmu4 + dt_4 * sqrt ( vh_frcst4 ));
% Riskmetrics
dh_frcstES = forecast (MdlES , 1,'Y0 ',yr .^2) ;
vh_frcst_es (i +1) = dh_frcstES ; % forecast about
volatility
vVarES = -( sqrt ( vh_frcst_es )*dt);
end
ve1 = y(T+1: end -1) .^2 - vh_frcst1 ; % forecast error Gaussian GARCH
ve2 = y(T+1: end -1) .^2 - vh_frcst2 ; % forecast error Gaussian EGARCH
ve3 = y(T+1: end -1) .^2 - vh_frcst3 ; % forecast error Student -t GARCH
ve4 = y(T+1: end -1) .^2 - vh_frcst4 ; % forecast error Student -t GJR
veES = y(T+1: end -1) .^2 - vh_frcst_es ; % forecast error Riskmetrics
dates_forecast = dates (T+1: end -1);
figure ('Name ', 'Comparison Volatility Prediction with Benchmark')
plot ( dates_forecast ,[y(T+1: end -1) .^2 , vh_frcst1 , vh_frcst2 , vh_frcst3 , vh_frcst4 ]);
hold on
plot ( dates_forecast , vh_frcst_es ,'LineWidth ' ,1,'color ', 'r')
xlabel ('years ')
ylabel ('squared returns ')
11;
xlim ([ dates_forecast (1) , dates_forecast (end )])
legend ('Squared Returns ','Gaussian GARCH ','Gaussian EGARCH', 't- GARCH ', 't-GJR - GARCH ', 'Benchmark ')
hold off
figure ('Name ','Comparison Forecast Errors ')
plot ( dates_forecast ,[ ve1 ,ve2 , ve3 , ve4 , veES ])
xlabel ('years ')
ylabel ('forecast errors ')
xlim ([ dates_forecast (1) , dates_forecast (end )])
legend ('Gaussian GARCH ', 'Gaussian EGARCH ', 't- GARCH ','t-GJR - GARCH '; 'Benchmark ')
mean ([ ve1 .^2 , ve2 .^2 , ve3 .^2 , ve4 .^2])
figure ('Name ','Comparison VaR ')
plot ( dates_forecast ,[y(T+1: end -1) ,- vVar1 , - vVar2 , - vVar3 , - vVar4 ])
hold on
plot ( dates_forecast ,- vVarES ,'LineWidth ' ,2,'color ','r')
xlabel ('years ')
ylabel ('returns ')
xlim ([ dates_forecast (1) , dates_forecast (end )])
legend ('Returns ','Gaussian GARCH ', 'Gaussian EGARCH ','t- GARCH '; 't-GJR - GARCH '; 'Benchmark ')
hold off
expected_exceed = m *0.01;
actual_exceed_1 = sum(y(T+1: end -1) < - vVar1 );
actual_exceed_2 = sum(y(T+1: end -1) < - vVar2 );
actual_exceed_3 = sum(y(T+1: end -1) < - vVar3 ); % less exceeds since t- student
actual_exceed_4 = sum(y(T+1: end -1) < - vVar4 ); % less exceeds since t- student
actual_exceed_ES = sum (y(T+1: end -1) < - vVarES );
idx_1 = y(T+1: end -1) < - vVar1 ;
returns_1 = y(T+1: end -1) .* idx_1 ;
returns_1 ( returns_1 == 0) = nan; % We find the returns that exceed the VAR
idx_2 = y(T+1: end -1) < - vVar2 ;
returns_2 = y(T+1: end -1) .* idx_2 ;
returns_2 ( returns_2 == 0) = nan; % We find the returns that exceed the VAR
12;
idx_3 = y(T+1: end -1) < - vVar3 ;
returns_3 = y(T+1: end -1) .* idx_3 ;
returns_3 ( returns_3 == 0) = nan; % We find the returns that exceed the VAR
idx_4 = y(T+1: end -1) < - vVar4 ;
returns_4 = y(T+1: end -1) .* idx_4 ;
returns_4 ( returns_4 == 0) = nan; % We find the returns that exceed the VAR
idx_ES = y(T+1: end -1) < - vVarES ;
returns_es = y(T+1: end -1) .* idx_ES ;
returns_es ( returns_es == 0) = nan; % We find the returns that exceed the VAR
subplot (5 ,1 ,1)
plot ( dates_forecast , y(T+1: end -1) )
hold on
plot ( dates_forecast , - vVarES ),
hold on
plot ( dates_forecast , returns_es ,'r*')
legend ('Returns ' , 'Var Riskmetrics ', 'Returns that exceed Var ')
xlabel ('years ')
ylabel ('returns ')
hold off
subplot (5 ,1 ,2)
plot ( dates_forecast , y(T+1: end -1) )
hold on
plot ( dates_forecast , - vVar1 )
hold on
plot ( dates_forecast , returns_1 ,'r*')
legend ('Returns ' , 'VaR Gaussian GARCH ','Returns that exceed Var')
xlabel ('years ')
ylabel ('returns ')
hold off
subplot (5 ,1 ,3)
plot ( dates_forecast , y(T+1: end -1) )
hold on
plot ( dates_forecast , - vVar2 )
hold on
plot ( dates_forecast , returns_2 ,'r*')
13;
legend ('Returns ' , 'VaR Gaussian EGARCH ','Returns that exceed Var')
xlabel ('years ')
ylabel ('returns ')
hold off
subplot (5 ,1 ,4)
plot ( dates_forecast , y(T+1: end -1) )
hold on
plot ( dates_forecast , - vVar3 )
hold on
plot ( dates_forecast , returns_3 ,'r*')
legend ('Returns ' , 'VaR T- Student GARCH ','Returns that exceed Var')
xlabel ('years ')
ylabel ('returns ')
hold off
subplot (5 ,1 ,5)
plot ( dates_forecast , y(T+1: end -1) )
hold on
plot ( dates_forecast , - vVar4 )
hold on
plot ( dates_forecast , returns_4 ,'r*')
legend ('Returns ' , 'VaR T- Student GJR - GARCH ','Returns that exceed Var')
xlabel ('years ')
ylabel ('returns ')
hold off
%% LOSS
y1 = y(T +1: end -1) ;
% Gaussian GARCH
loss_1 = fCheckFunction (y1 + vVar1 , dp);
total_loss_1 = sum ( loss_1 );
% Gaussian EGARCH
loss_2 = fCheckFunction (y1 + vVar2 , dp);
total_loss_2 = sum ( loss_2 );
% T- Student GARCH
loss_3 = fCheckFunction (y1 + vVar3 , dp);
total_loss_3 = sum ( loss_3 );
% T- Student GJR - GARCH
14;
loss_4 = fCheckFunction (y1 + vVar4 , dp);
total_loss_4 = sum ( loss_4 );
% Riskmetrics
loss_ES = fCheckFunction (y1 + vVarES , dp);
total_loss_ES = sum( loss_ES );
%% Diebold - Mariano test
m_DM = m ^(1/3) ;
omega = 0;
% Gaussian GARCH
d_1 = loss_1 - loss_ES ;
D_1 = mean (d_1);
spec_dens_1 = 2* pi* fBartlettSpectralDensityEst (d_1 , m_DM , omega );
DM_1 = D_1 ./( sqrt ( spec_dens_1 ./m)); % Result of Diebold - Mariano test
pvalue_1 = normcdf (DM_1 ,dp); % pvalue of Diebold - Mariano test
% Gaussian EGARCH
d_2 = loss_2 - loss_ES ;
D_2 = mean (d_2);
spec_dens_2 = 2* pi* fBartlettSpectralDensityEst (d_2 , m_DM , omega );
DM_2 = D_2 ./( sqrt ( spec_dens_2 ./m)); % Result of Diebold - Mariano test
pvalue_2 = normcdf (DM_2 ,dp); % pvalue of Diebold - Mariano test
% T- Student GARCH
d_3 = loss_3 - loss_ES ;
D_3 = mean (d_3);
spec_dens_3 = 2* pi* fBartlettSpectralDensityEst (d_3 , m_DM , omega );
DM_3 = D_3 ./( sqrt ( spec_dens_3 ./m)); % Result of Diebold - Mariano test
pvalue_3 = normcdf (DM_3 ,dp); % pvalue of Diebold - Mariano test
% T- Student GJR - GARCH
d_4 = loss_4 - loss_ES ;
D_4 = mean (d_4);
spec_dens_4 = 2* pi* fBartlettSpectralDensityEst (d_4 , m_DM , omega );
15;
DM_4 = D_4 ./( sqrt ( spec_dens_4 ./m)); % Result of Diebold - Mariano test
pvalue_4 = normcdf (DM_4 ,dp); % pvalue of Diebold - Mariano test
%% MSFE
sqrd_returns = y1 .^2;
MSFE_ES = 1/m*sum (( sqrd_returns - vh_frcst_es ) .^2) ; % Mean squared forecast error of Riskmetrics
MSFE_1 = 1/m*sum (( sqrd_returns - vh_frcst1 ) .^2) ; % Mean squared forecast error of Gaussian Garch
MSFE_2 = 1/m*sum (( sqrd_returns -vh_frcst2 ) .^2) ; % Mean squared forecast error of Gaussian EGARCH
MSFE_3 = 1/m*sum (( sqrd_returns -vh_frcst3 ) .^2) ; % Mean squared forecast error of T- Student Garch
MSFE_4 = 1/m*sum (( sqrd_returns -vh_frcst4 ) .^2) ; % Mean squared forecast error of T- Student GJR - GARCH