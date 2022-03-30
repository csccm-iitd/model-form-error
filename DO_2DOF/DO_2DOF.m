clc
clear
close all
rng('default')

%% Part I

global dt;
load ChiChi_X.txt; eq  = ChiChi_X; 
dt_new = 0.005; t = eq(:,1); ug = eq(:,2); dt = t(2)-t(1);
tsin = timeseries(t,t(1):dt:t(end));
tsout = resample(tsin,t(1):dt_new:t(end));
t = tsout.Data;
tsin = timeseries(ug,t(1):dt:t(end));
tsout = resample(tsin,t(1):dt_new:t(end));
f_e = tsout.Data'; dt = t(2)-t(1);

for i = 1:3
    fs = 1/dt; f = randn(1,length(f_e));
    f = (max(abs(f_e))/max(abs(f_e)))*f;
    f = bp(detrend(f),fs,0.5,4); f_ug(i,:) = f;
end

%% Part II

ug = f_ug(3,:)*9.81; figure; plot(t,ug);
plot_properties('','\textbf{Time (t)}','\textbf{Gr. Acc.}','Ground Acceleration - Data Simultaion',1,1);

m1 = 30; c1 = 10; k1 = 1000; f1 = -m1*ug;
m2 = 15; c2 = 5; k2 = 1000; f2 = -m2*ug;
alpha_do = 100; % alpha_do_1 = alpha_do_2 = 100;

xdv(:,1) = [0; 0; 0; 0];
for i = 2:length(t)
    y1 = xdv(1,i-1); y2 = xdv(2,i-1); y3 = xdv(3,i-1); y4 = xdv(4,i-1);
    xdv(:,i) = disp_vel_IT(y1, y2, y3, y4, f1(i-1), f2(i-1), f_ug(1,i-1), f_ug(2,i-1), k1, k2, c1, c2, m1, m2);
end

xa(:,1) = [0; 0];
for i = 1:length(t)
    xa(:,i) = [(1/m1)*(-(c1*xdv(3,i)+k1*xdv(1,i)+alpha_do*(xdv(1,i))^3+...
        c2*(xdv(3,i)-xdv(4,i))+k2*(xdv(1,i)-xdv(2,i))+alpha_do*(xdv(1,i)-xdv(2,i))^3));
        (1/m2)*(-(c2*(xdv(4,i)-xdv(3,i))+k2*(xdv(2,i)-xdv(1,i))+alpha_do*(xdv(2,i)-xdv(1,i))^3))];
end

c1 = c1*1.2; k1 = k1*0.90; c2 = c2*0.9; k2 = k2*0.85;

M = [m1, 0; 0, m2]; K = [k1+k2, -k2; -k2, k2]; C = [c1+c2, -c2; -c2, c2];
[~, wnsq] = eig(K,M); wn = diag(wnsq).^0.5;

Ac = [zeros(2), eye(2); -inv(M)*K, -inv(M)*C]; Bc = [zeros(2); inv(M)];
Ad = expm(Ac*dt); Bd = (Ad-eye(4))*inv(Ac)*Bc;

xdv_lin(:,1) = [0; 0; 0; 0];
for i = 2:length(t)
    xdv_lin(:,i) = Ad*xdv_lin(:,i-1)+Bd*[f1(i-1); f2(i-1)];
end

snr = 1e10;
z(1,:) = xa(1,:)+sqrt((1/snr)*(std(xa(1,:))^2))*randn(1,length(xa(1,:)));
z(2,:) = xa(2,:)+sqrt((1/snr)*(std(xa(2,:))^2))*randn(1,length(xa(2,:)));

figure;
subplot(2,2,1); hold on; plot(t,xdv_lin(1,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(1,:),'r');
plot_properties({'Linear-Linear','Duffing Oscillator-Duffing Oscillator'},'\textbf{Time (t)}','\textbf{Disp.}','Displacement Mass-1',1,1);
subplot(2,2,2); hold on; plot(t,xdv_lin(2,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(2,:),'r');
plot_properties('','\textbf{Time (t)}','\textbf{Disp.}','Displacement Mass-2',1,1);
subplot(2,2,3); hold on; plot(t,xdv_lin(3,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(3,:),'r');
plot_properties('','\textbf{Time (t)}','\textbf{Vel.}','Velocity Mass-1',1,1);
subplot(2,2,4); hold on; plot(t,xdv_lin(4,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(4,:),'r');
plot_properties('','\textbf{Time (t)}','\textbf{Vel.}','Velocity Mass-2',1,1);
sub_title('Data Simulation',0.515,0.96);

%% Part III

M = [m1, 0; 0, m2]; K = [k1+k2, -k2; -k2, k2];
C = [c1+c2, -c2; -c2, c2];
Ac = [zeros(2), eye(2); -inv(M)*K, -inv(M)*C]; Bc = [zeros(2); inv(M)]; Cc = [zeros(2); -inv(M)];
Ad = expm(Ac*dt); Bd = (Ad-eye(4))*inv(Ac)*Bc; Cd = (Ad-eye(4))*inv(Ac)*Cc;

Am = [-inv(M)*K, -inv(M)*C]; Bm = inv(M); Cm = -inv(M);

ps = [0; 0; 0; 0]; pps = diag([1 1 1 1]); pss = eye(4).*1e-4;
pi = [0; 0]; ppi = diag([1 1]); pii = 100*[1 0; 0 1];

R = eye(2).*1e-10;
s = zeros(4,length(t)); nl = zeros(2,length(t));

for i = 1:length(t)
    ei = z(:,i) - Am*ps - Cm*pi;
    kgi = (inv(Cm*ppi*Cm'+R)*Cm*ppi)';
    
    ci = pi+kgi*ei;
    pci = ppi - kgi*Cm*ppi;
    
    es = z(:,i) - Am*ps - Cm*ci;
    kgs = (inv(Am*pps*Am'+R)*Am*pps)';
    
    cs = ps+kgs*es;
    pcs = pps - kgs*Am*pps;
    
    pi = ci;
    ppi = pci+pii;
    
    ps = Ad*cs+Bd*[f1(i); f2(i)]+Cd*ci;
    pps = Ad*pcs*Ad'+pss;
    
    s(:,i) = cs; nl(:,i) = ci;
end

figure;
subplot(2,3,1); plot(t,xdv(1,:),'r'); hold on; plot(t,s(1,:),'b--');
plot_properties({'True','Filtered'},'\textbf{Time (t)}','\textbf{Disp.}','Displacement Mass-1',1,1);
subplot(2,3,2); plot(t,xdv(2,:),'r'); hold on; plot(t,s(2,:),'b--');
plot_properties('','\textbf{Time (t)}','\textbf{Disp.}','Displacement Mass-2',1,1);
subplot(2,3,3); plot(t,xdv(3,:),'r'); hold on; plot(t,s(3,:),'b--');
plot_properties('','\textbf{Time (t)}','\textbf{Vel.}','Velocity Mass-1',1,1);
subplot(2,3,4); plot(t,xdv(4,:),'r'); hold on; plot(t,s(4,:),'b--');
plot_properties('','\textbf{Time (t)}','\textbf{Vel.}','Velocity Mass-2',1,1);
subplot(2,3,5); plot(t,(-(m1*xa(1,:)+c1*xdv(3,:)+c2*xdv(3,:)-c2*xdv(4,:)+k1*xdv(1,:)+k2*xdv(1,:)-k2*xdv(2,:))),'r'); hold on; plot(t,detrend(nl(1,:)),'b--');
plot_properties('','\textbf{Time (t)}','$\mathbf{R_1}$','Non-linearity Mass-1',1,1);
subplot(2,3,6); plot(t,(-(m2*xa(2,:)-c2*xdv(3,:)+c2*xdv(4,:)-k2*xdv(1,:)+k2*xdv(2,:))),'r'); hold on; plot(t,detrend(nl(2,:)),'b--');
plot_properties('','\textbf{Time (t)}','$\mathbf{R_2}$','Non-linearity Mass-2',1,1);
sub_title('Filter Results',0.515,0.96);

%% Part IV a

tr_lim = 15;

mdl1  = fitrgp(s(:,1:tr_lim/dt)',nl(1,1:tr_lim/dt)',...
    'FitMethod','sr','PredictMethod','exact','Standardize',1);
[a1,~,ci1] = predict(mdl1,s');

mdl2  = fitrgp(s(:,1:tr_lim/dt)',nl(2,1:tr_lim/dt)',...
    'FitMethod','sr','PredictMethod','exact','Standardize',1);
[a2,~,ci2] = predict(mdl2,s');

figure;
subplot(2,1,1); plot(t,(-(m1*xa(1,:)+c1*xdv(3,:)+c2*xdv(3,:)-c2*xdv(4,:)+k1*xdv(1,:)+k2*xdv(1,:)-k2*xdv(2,:))),'b','linewidth',1.5);
hold on; plot(t,a1,'m--','linewidth',1.5);
xline(t(tr_lim/dt),'k-.','linewidth',2.5); xlim([0 t(end)]);
plot_properties({'Ground Truth','GP Estimates','Last Training Data'},'\textbf{Time (t)}','$\mathbf{R_1}$','',1,1);

subplot(2,1,2); plot(t,(-(m2*xa(2,:)-c2*xdv(3,:)+c2*xdv(4,:)-k2*xdv(1,:)+k2*xdv(2,:))),'b','linewidth',1.5);
hold on; plot(t,a2,'m--','linewidth',1.5);
xline(t(tr_lim/dt),'k-.','linewidth',2.5); xlim([0 t(end)]);
plot_properties({'Ground Truth','GP Estimates','Last Training Data'},'\textbf{Time (t)}','$\mathbf{R_2}$','',1,1);

actual_nl_1 = (-(m1*xa(1,:)+c1*xdv(3,:)+c2*xdv(3,:)-c2*xdv(4,:)+k1*xdv(1,:)+k2*xdv(1,:)-k2*xdv(2,:)));
actual_nl_2 = (-(m2*xa(2,:)-c2*xdv(3,:)+c2*xdv(4,:)-k2*xdv(1,:)+k2*xdv(2,:)));
rmse(1,1) = sqrt(mean((actual_nl_1-a1').^2));
rmse(2,1) = sqrt(mean((actual_nl_2-a2').^2));

clear mdl1 mdl2 a1 a2

%% Part IV b

tr_lim = 25;

mdl1  = fitrgp(s(:,1:tr_lim/dt)',nl(1,1:tr_lim/dt)',...
    'FitMethod','sr','PredictMethod','exact','Standardize',1);
a1 = predict(mdl1,s');

mdl2  = fitrgp(s(:,1:tr_lim/dt)',nl(2,1:tr_lim/dt)',...
    'FitMethod','sr','PredictMethod','exact','Standardize',1);
a2 = predict(mdl2,s');

figure;
subplot(2,1,1); plot(t,(-(m1*xa(1,:)+c1*xdv(3,:)+c2*xdv(3,:)-c2*xdv(4,:)+k1*xdv(1,:)+k2*xdv(1,:)-k2*xdv(2,:))),'b','linewidth',1.5);
hold on; plot(t,a1,'m--','linewidth',1.5);
xline(t(tr_lim/dt),'k-.','linewidth',2.5); xlim([0 t(end)]);
plot_properties({'Ground Truth','GP Estimates','Last Training Data'},'\textbf{Time (t)}','$\mathbf{R_1}$','',1,1);

subplot(2,1,2); plot(t,(-(m2*xa(2,:)-c2*xdv(3,:)+c2*xdv(4,:)-k2*xdv(1,:)+k2*xdv(2,:))),'b','linewidth',1.5);
hold on; plot(t,a2,'m--','linewidth',1.5);
xline(t(tr_lim/dt),'k-.','linewidth',2.5); xlim([0 t(end)]);
plot_properties({'Ground Truth','GP Estimates','Last Training Data'},'\textbf{Time (t)}','$\mathbf{R_2}$','',1,1);

actual_nl_1 = (-(m1*xa(1,:)+c1*xdv(3,:)+c2*xdv(3,:)-c2*xdv(4,:)+k1*xdv(1,:)+k2*xdv(1,:)-k2*xdv(2,:)));
actual_nl_2 = (-(m2*xa(2,:)-c2*xdv(3,:)+c2*xdv(4,:)-k2*xdv(1,:)+k2*xdv(2,:)));
rmse(1,2) = sqrt(mean((actual_nl_1-a1').^2));
rmse(2,2) = sqrt(mean((actual_nl_2-a2').^2));

clear mdl1 mdl2 a1 a2

%% Part IV c

tr_lim = 35;

mdl1  = fitrgp(s(:,1:tr_lim/dt)',nl(1,1:tr_lim/dt)',...
    'FitMethod','sr','PredictMethod','exact','Standardize',1);
a1 = predict(mdl1,s');

mdl2  = fitrgp(s(:,1:tr_lim/dt)',nl(2,1:tr_lim/dt)',...
    'FitMethod','sr','PredictMethod','exact','Standardize',1);
a2 = predict(mdl2,s');

figure;
subplot(2,1,1); plot(t,(-(m1*xa(1,:)+c1*xdv(3,:)+c2*xdv(3,:)-c2*xdv(4,:)+k1*xdv(1,:)+k2*xdv(1,:)-k2*xdv(2,:))),'b','linewidth',1.5);
hold on; plot(t,a1,'m--','linewidth',1.5);
xline(t(tr_lim/dt),'k-.','linewidth',2.5); xlim([0 t(end)]);
plot_properties({'Ground Truth','GP Estimates','Last Training Data'},'\textbf{Time (t)}','$\mathbf{R_1}$','',1,1);

subplot(2,1,2); plot(t,(-(m2*xa(2,:)-c2*xdv(3,:)+c2*xdv(4,:)-k2*xdv(1,:)+k2*xdv(2,:))),'b','linewidth',1.5);
hold on; plot(t,a2,'m--','linewidth',1.5);
xline(t(tr_lim/dt),'k-.','linewidth',2.5); xlim([0 t(end)]);
plot_properties({'Ground Truth','GP Estimates','Last Training Data'},'\textbf{Time (t)}','$\mathbf{R_2}$','',1,1);

actual_nl_1 = (-(m1*xa(1,:)+c1*xdv(3,:)+c2*xdv(3,:)-c2*xdv(4,:)+k1*xdv(1,:)+k2*xdv(1,:)-k2*xdv(2,:)));
actual_nl_2 = (-(m2*xa(2,:)-c2*xdv(3,:)+c2*xdv(4,:)-k2*xdv(1,:)+k2*xdv(2,:)));
rmse(1,3) = sqrt(mean((actual_nl_1-a1').^2));
rmse(2,3) = sqrt(mean((actual_nl_2-a2').^2));

%% Part V

xdv_s_gp(:,1) = [0; 0; 0; 0];
for i = 2:length(t)
    disp(['step ',num2str(i)]);
    av1 = predict(mdl1,[xdv_s_gp(1,i-1), xdv_s_gp(2,i-1), xdv_s_gp(3,i-1), xdv_s_gp(4,i-1)]);
    av2 = predict(mdl2,[xdv_s_gp(1,i-1), xdv_s_gp(2,i-1), xdv_s_gp(3,i-1), xdv_s_gp(4,i-1)]);
    
    xdv_s_gp(:,i) = Ad*xdv_s_gp(:,i-1)+Bd*[f1(i-1); f2(i-1)]+Cd*[av1; av2];
end

figure;
subplot(2,1,1); hold on; plot(t,xdv(1,:),'r','linewidth',1.5);
plot(t,xdv_s_gp(1,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties({'Ground Truth','Projected States'},'\textbf{Time (t)}','$\mathbf{x_1}$','',1,1);

subplot(2,1,2); hold on; plot(t,xdv(2,:),'r','linewidth',1.5);
plot(t,xdv_s_gp(2,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{x_2}$','',1,1);

figure;
subplot(2,1,1); hold on; plot(t,xdv(3,:),'r','linewidth',1.5);
plot(t,xdv_s_gp(3,:),'b--','linewidth',1.5);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}_1}$','',1,1);

subplot(2,1,2); hold on; plot(t,xdv(4,:),'r','linewidth',1.5);
plot(t,xdv_s_gp(4,:),'b--','linewidth',1.5);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}_2}$','',1,1);

%% Part VI a

load ChiChi_X.txt; eq  = ChiChi_X; 
dt_new = 0.005; t = eq(:,1); ug = eq(:,2); dt = t(2)-t(1);
tsin = timeseries(t,t(1):dt:t(end));
tsout = resample(tsin,t(1):dt_new:t(end));
t = tsout.Data;
tsin = timeseries(ug,t(1):dt:t(end));
tsout = resample(tsin,t(1):dt_new:t(end));
f_e = tsout.Data'; dt = t(2)-t(1);

clear f_ug
for i = 1:2
    fs = 1/dt; f = randn(1,length(f_e));
    f = (max(abs(f_e))/max(abs(f_e)))*f;
    f = bp(detrend(f),fs,0.5,4);
    L = length(f); f_ug(i,:) = f;
end

ug = f_e*9.81; figure; plot(t,ug);
plot_properties('','\textbf{Time (t)}','\textbf{Gr. Acc.}','Ground Acceleration - Testing',1,1);

f1 = -m1*ug; f2 = -m2*ug;

xdv_d_gp(:,1) = [0; 0; 0; 0];
for i = 2:length(t)
    disp(['step ',num2str(i)]);
    av1 = predict(mdl1,[xdv_d_gp(1,i-1), xdv_d_gp(2,i-1), xdv_d_gp(3,i-1), xdv_d_gp(4,i-1)]);
    av2 = predict(mdl2,[xdv_d_gp(1,i-1), xdv_d_gp(2,i-1), xdv_d_gp(3,i-1), xdv_d_gp(4,i-1)]);
    
    xdv_d_gp(:,i) = Ad*xdv_d_gp(:,i-1)+Bd*[f1(i-1); f2(i-1)]+Cd*[av1; av2];
end

m1 = 30; c1 = 10; k1 = 1000; f1 = -m1*ug;
m2 = 15; c2 = 5; k2 = 1000; f2 = -m2*ug; alpha_do = 100;

xdv_d_tr(:,1) = [0; 0; 0; 0];
for i = 2:length(t)
    y1 = xdv_d_tr(1,i-1); y2 = xdv_d_tr(2,i-1); y3 = xdv_d_tr(3,i-1); y4 = xdv_d_tr(4,i-1);
    xdv_d_tr(:,i) = disp_vel_IT(y1, y2, y3, y4, f1(i-1), f2(i-1), f_ug(1,i-1), f_ug(2,i-1), k1, k2, c1, c2, m1, m2);
end

figure;
subplot(2,1,1); hold on; plot(t,xdv_d_tr(1,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(1,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties({'Ground Truth','Projected States'},'\textbf{Time (t)}','$\mathbf{x_1}$','',1,1);

subplot(2,1,2); hold on; plot(t,xdv_d_tr(2,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(2,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{x_2}$','',1,1);

figure;
subplot(2,1,1); hold on; plot(t,xdv_d_tr(3,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(3,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}_1}$','',1,1);

subplot(2,1,2); hold on; plot(t,xdv_d_tr(4,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(4,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}_2}$','',1,1);

%% Part VI b

m1 = 30; c1 = 10; k1 = 1000; m2 = 15; c2 = 5; k2 = 1000;
c1 = c1*1.2; k1 = k1*0.90; c2 = c2*0.9; k2 = k2*0.85;

M = [m1, 0; 0, m2]; K = [k1+k2, -k2; -k2, k2]; C = [c1+c2, -c2; -c2, c2];
[ph, wnsq] = eig(K,M); wn = diag(wnsq).^0.5;

Ac = [zeros(2), eye(2); -inv(M)*K, -inv(M)*C]; Bc = [zeros(2); inv(M)];
Ad = expm(Ac*dt); Bd = (Ad-eye(4))*inv(Ac)*Bc;

load ChiChi_X.txt; eq  = ChiChi_X; 
dt_new = 0.005; t = eq(:,1); ug = eq(:,2); dt = t(2)-t(1);
tsin = timeseries(t,t(1):dt:t(end));
tsout = resample(tsin,t(1):dt_new:t(end));
t = tsout.Data;
tsin = timeseries(ug,t(1):dt:t(end));
tsout = resample(tsin,t(1):dt_new:t(end));
f_e = tsout.Data'; dt = t(2)-t(1);

clear f_ug
for i = 1:3
    fs = 1/dt; f = randn(1,length(f_e));
    f = (max(abs(f_e))/max(abs(f_e)))*f;
    f = bp(detrend(f),fs,0.5,4);
    L = length(f); f_ug(i,:) = f;
end

ug = 0.75*f_ug(3,:)*9.81; figure; plot(t,ug);
plot_properties('','\textbf{Time (t)}','\textbf{Gr. Acc.}','Ground Acceleration - Testing',1,1);

f1 = -m1*ug; f2 = -m2*ug;

xdv_d_gp(:,1) = [0; 0; 0; 0];
for i = 2:length(t)
    disp(['step ',num2str(i)]);
    av1 = predict(mdl1,[xdv_d_gp(1,i-1), xdv_d_gp(2,i-1), xdv_d_gp(3,i-1), xdv_d_gp(4,i-1)]);
    av2 = predict(mdl2,[xdv_d_gp(1,i-1), xdv_d_gp(2,i-1), xdv_d_gp(3,i-1), xdv_d_gp(4,i-1)]);
    
    xdv_d_gp(:,i) = Ad*xdv_d_gp(:,i-1)+Bd*[f1(i-1); f2(i-1)]+Cd*[av1; av2];
end

m1 = 30; c1 = 10; k1 = 1000; f1 = -m1*ug;
m2 = 15; c2 = 5; k2 = 1000; f2 = -m2*ug; alpha_do = 100;

xdv_d_tr(:,1) = [0; 0; 0; 0];
for i = 2:length(t)
    y1 = xdv_d_tr(1,i-1); y2 = xdv_d_tr(2,i-1); y3 = xdv_d_tr(3,i-1); y4 = xdv_d_tr(4,i-1);
    xdv_d_tr(:,i) = disp_vel_IT(y1, y2, y3, y4, f1(i-1), f2(i-1), f_ug(1,i-1), f_ug(2,i-1), k1, k2, c1, c2, m1, m2);
end

figure;
subplot(2,1,1); hold on; plot(t,xdv_d_tr(1,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(1,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties({'Ground Truth','Projected States'},'\textbf{Time (t)}','$\mathbf{x_1}$','',1,1);

subplot(2,1,2); hold on; plot(t,xdv_d_tr(2,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(2,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{x_2}$','',1,1);

figure;
subplot(2,1,1); hold on; plot(t,xdv_d_tr(3,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(3,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}_1}$','',1,1);

subplot(2,1,2); hold on; plot(t,xdv_d_tr(4,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(4,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}_2}$','',1,1);

%% Part VI c

m1 = 30; c1 = 10; k1 = 1000; m2 = 15; c2 = 5; k2 = 1000;
c1 = c1*1.2; k1 = k1*0.90; c2 = c2*0.9; k2 = k2*0.85;

M = [m1, 0; 0, m2]; K = [k1+k2, -k2; -k2, k2]; C = [c1+c2, -c2; -c2, c2];
[ph, wnsq] = eig(K,M); wn = diag(wnsq).^0.5;

Ac = [zeros(2), eye(2); -inv(M)*K, -inv(M)*C]; Bc = [zeros(2); inv(M)];
Ad = expm(Ac*dt); Bd = (Ad-eye(4))*inv(Ac)*Bc;

load ChiChi_X.txt; eq  = ChiChi_X; 
dt_new = 0.005; t = eq(:,1); ug = eq(:,2); dt = t(2)-t(1);
tsin = timeseries(t,t(1):dt:t(end));
tsout = resample(tsin,t(1):dt_new:t(end));
t = tsout.Data;
tsin = timeseries(ug,t(1):dt:t(end));
tsout = resample(tsin,t(1):dt_new:t(end));
f_e = tsout.Data'; dt = t(2)-t(1);

clear f_ug
for i = 1:3
    fs = 1/dt; f = randn(1,length(f_e));
    f = (max(abs(f_e))/max(abs(f_e)))*f;
    f = bp(detrend(f),fs,0.5,4);
    L = length(f); f_ug(i,:) = f;
end

ug = 1.25*f_ug(3,:)*9.81; figure; plot(t,ug);
plot_properties('','\textbf{Time (t)}','\textbf{Gr. Acc.}','Ground Acceleration - Testing',1,1);

f1 = -m1*ug; f2 = -m2*ug;

xdv_d_gp(:,1) = [0; 0; 0; 0];
for i = 2:length(t)
    disp(['step ',num2str(i)]);
    av1 = predict(mdl1,[xdv_d_gp(1,i-1), xdv_d_gp(2,i-1), xdv_d_gp(3,i-1), xdv_d_gp(4,i-1)]);
    av2 = predict(mdl2,[xdv_d_gp(1,i-1), xdv_d_gp(2,i-1), xdv_d_gp(3,i-1), xdv_d_gp(4,i-1)]);
    
    xdv_d_gp(:,i) = Ad*xdv_d_gp(:,i-1)+Bd*[f1(i-1); f2(i-1)]+Cd*[av1; av2];
end

m1 = 30; c1 = 10; k1 = 1000; f1 = -m1*ug;
m2 = 15; c2 = 5; k2 = 1000; f2 = -m2*ug; alpha_do = 100;

xdv_d_tr(:,1) = [0; 0; 0; 0];
for i = 2:length(t)
    y1 = xdv_d_tr(1,i-1); y2 = xdv_d_tr(2,i-1); y3 = xdv_d_tr(3,i-1); y4 = xdv_d_tr(4,i-1);
    xdv_d_tr(:,i) = disp_vel_IT(y1, y2, y3, y4, f1(i-1), f2(i-1), f_ug(1,i-1), f_ug(2,i-1), k1, k2, c1, c2, m1, m2);
end

figure;
subplot(2,1,1); hold on; plot(t,xdv_d_tr(1,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(1,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties({'Ground Truth','Projected States'},'\textbf{Time (t)}','$\mathbf{x_1}$','',1,1);

subplot(2,1,2); hold on; plot(t,xdv_d_tr(2,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(2,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{x_2}$','',1,1);

figure;
subplot(2,1,1); hold on; plot(t,xdv_d_tr(3,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(3,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}_1}$','',1,1);

subplot(2,1,2); hold on; plot(t,xdv_d_tr(4,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(4,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}_2}$','',1,1);

%%

figure;
hold on; plot([15, 20, 35],rmse(1,:),'k-o','linewidth',2);
plot([15, 20, 35],rmse(2,:),'m--*','linewidth',2);
plot_properties({'I-DOF','II-DOF'},'\textbf{Training length (sec)}','RMSE','',1,1);
xlim([10,45]); ylim([min(min(rmse))-1, 1+max(max(rmse))]);

%%  Part VII - Functions

function y = disp_vel_IT(y1, y2, y3, y4, f1, f2, r1, r2, k1, k2, c1, c2, m1, m2)
global dt

alpha_do = 100; sigma1 = 0.05; sigma2 = 0.05;

y = [y1; y2; y3; y4];
a = [y3;
    y4;
    (1/m1)*(f1-(c1*y3+k1*y1+alpha_do*(y1)^3+c2*(y3-y4)+k2*(y1-y2)+alpha_do*(y1-y2)^3));
    (1/m2)*(f2-(c2*(y4-y3)+k2*(y2-y1)+alpha_do*(y2-y1)^3))];
L0a = [-(c1*y3 - f1 + k1*y1 + c2*(y3 - y4) + k2*(y1 - y2) + alpha_do*y1^3 + alpha_do*(y1 - y2)^3)/m1;
    (f2 + c2*(y3 - y4) + k2*(y1 - y2) + alpha_do*(y1 - y2)^3)/m2;
    (y4*(k2 + 3*alpha_do*(y1 - y2)^2))/m1 - (y3*(k1 + k2 + 3*alpha_do*y1^2 + 3*alpha_do*(y1 - y2)^2))/m1 + ((c1 + c2)*(c1*y3 - f1 + k1*y1 + c2*(y3 - y4) + k2*(y1 - y2) + alpha_do*y1^3 + alpha_do*(y1 - y2)^3))/m1^2 + (c2*(f2 + c2*(y3 - y4) + k2*(y1 - y2) + alpha_do*(y1 - y2)^3))/(m1*m2);
    (y3*(k2 + 3*alpha_do*(y1 - y2)^2))/m2 - (c2*(f2 + c2*(y3 - y4) + k2*(y1 - y2) + alpha_do*(y1 - y2)^3))/m2^2 - (y4*(k2 + 3*alpha_do*(y1 - y2)^2))/m2 - (c2*(c1*y3 - f1 + k1*y1 + c2*(y3 - y4) + k2*(y1 - y2) + alpha_do*y1^3 + alpha_do*(y1 - y2)^3))/(m1*m2)];
L1a = [sigma1/m1;
    sigma2/m2;
    (c2*sigma2)/(m1*m2) - (sigma1*(c1 + c2))/m1^2;
    (c2*sigma1)/(m1*m2) - (c2*sigma2)/m2^2];
b = [0;0;sigma1/m1;sigma2/m2];
dz1 = (dt^(3/2)*r1)/2 + (3^(1/2)*dt^(3/2)*r2)/6; dw = dt^(1/2)*r1;

y = y+a*dt+b*dw+L1a*dz1+0.5*L0a*dt^2;
end

function [bp_data] = bp(data, fsample, lfc, hfc)
fnyq = fsample/2;
wl = lfc/fnyq;
wh = hfc/fnyq;
P = 2*ceil(8*4/(wh-wl));
Coeff = fir1(P,[wl, wh],'noscale');
bp_data = filter(Coeff, 1, data);
end

function plot_properties(legends,x_label,y_label,title_fig,index_lt,index_ht,legend_font,fig_font)

if nargin == 6
    legend_font = 16;
    fig_font = 20;
end
if nargin == 7
    fig_font = 20;
end

if isempty(legends)
    set(gca,'LineWidth',2,'FontSize',fig_font,'FontWeight','bold','FontName','Times');
    set(gcf,'Position',[1 1 index_lt*round(2000) index_ht*round(2000)]);
    set(get(gca,'xlabel'),'String',x_label,'FontSize',fig_font','Interpreter','Latex');
    set(get(gca,'ylabel'),'String',y_label,'FontSize',fig_font','Interpreter','Latex');
    set(gcf,'color','w'); box on; title(title_fig);
else
    legend(legends,'FontSize',legend_font,'FontWeight','bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','bold','FontName','Times');
    set(gcf,'Position',[1 1 index_lt*round(2000) index_ht*round(2000)]);
    set(get(gca,'xlabel'),'String',x_label,'FontSize',fig_font','Interpreter','Latex');
    set(get(gca,'ylabel'),'String',y_label,'FontSize',fig_font','Interpreter','Latex');
    set(gcf,'color','w'); box on; title(title_fig);
end
end

function sub_title(plot_title,x,y)
axes('Position',[x, y, 0.01, 0.01],'Xlim',[0 0.01],'Ylim',[0  0.01],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
text(0, 0, plot_title, 'FontSize', 20', 'FontWeight', 'Bold','HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
end
