clc
clear
close all
rng('default')

%% Part I

global dt;
load elcentro.txt; eq  = elcentro;
dt_new = 0.001; t = eq(:,1); ug = eq(:,2); dt = t(2)-t(1);
tsin = timeseries(t,t(1):dt:t(end));
tsout = resample(tsin,t(1):dt_new:t(end));
t = tsout.Data;
tsin = timeseries(ug,t(1):dt:t(end));
tsout = resample(tsin,t(1):dt_new:t(end));
f_e = tsout.Data'; dt = t(2)-t(1);

f_ug = zeros(3,length(t));
for i = 1:3
    fs = 1/dt; f = randn(1,length(t));
    f = (max(abs(f_e))/max(abs(f_e)))*f;
    f = bp(detrend(f),fs,0.5,4); f_ug(i,:) = f;
end

%% Part II

ug = f_ug(3,:); figure; plot(t,ug);
plot_properties('','\textbf{Time (t)}','Gr. Acc.','Ground Acceleration - Data Simultaion',1,1);

m = 20*[20; 19; 18; 17; 16];
c = 10*[10; 20; 19; 18; 17];
k = 10000*[10; 20; 19; 18; 17];

m1 = m(1); m2 = m(2); m3 = m(3); m4 = m(4); m5 = m(5);
c1 = c(1); c2 = c(2); c3 = c(3); c4 = c(4); c5 = c(5);
k1 = k(1); k2 = k(2); k3 = k(3); k4 = k(4); k5 = k(5);

g = 9.81; Qy = 0.05*(m1+m2+m3+m4+m5)*g; alpha = 1;
beta = 0.5; gamma = 0.5; eta = 1; Dy = 0.013; kr = 1/6; dofs = 5;

xdv(:,1) = [zeros(2*dofs,1); 1e-5];
for i = 2:length(t)
    y = xdv(:,i-1);
    xdv(:,i) = disp_vel_IT_5(y, ug(i-1), f_ug(1,i-1), f_ug(2,i-1), m, c, k);
end

xa(:,1) = zeros(dofs,1);
for i = 1:length(t)
    y1 = xdv(1,i); y2 = xdv(2,i); y3 = xdv(3,i); y4 = xdv(4,i); y5 = xdv(5,i);
    y6 = xdv(6,i); y7 = xdv(7,i); y8 = xdv(8,i); y9 = xdv(9,i); y10 = xdv(10,i);
    y11 = xdv(11,i);
    
    xa(:,i) = [(1/m1)*(-(c1*y6+k1*y1+c2*(y6-y7)+k2*(y1-y2)+(1-kr)*Qy*y11));
        (1/m2)*(-(c2*(y7-y6)+k2*(y2-y1)+c3*(y7-y8)+k3*(y2-y3)));
        (1/m3)*(-(c3*(y8-y7)+k3*(y3-y2)+c4*(y8-y9)+k4*(y3-y4)));
        (1/m4)*(-(c4*(y9-y8)+k4*(y4-y3)+c5*(y9-y10)+k5*(y4-y5)));
        (1/m5)*(-(c5*(y10-y9)+k5*(y5-y4)))];
end

snr = 1e10; z = zeros(size(xa));
for i = 1:dofs
    z(i,:) = xa(i,:)+sqrt((1/snr)*(std(xa(i,:))^2))*randn(1,length(xa(i,:)));
end

%%%%%%%%%% %%%%%%%%%% Linear System %%%%%%%%%% %%%%%%%%%%

c = [1.10*c1; 1.05*c2; 0.90*c3; 1.10*c4; 0.95*c5];
k = [1.05*k1; 1.05*k2; 0.95*k3; 0.95*k4; 0.95*k5];

c1 = c(1); c2 = c(2); c3 = c(3); c4 = c(4); c5 = c(5);
k1 = k(1); k2 = k(2); k3 = k(3); k4 = k(4); k5 = k(5);

M = diag(m);
K = diag(k)+diag([k(2:dofs);0]);
for i = 1:dofs-1
    for j = 1:dofs
        if j == i+1
            K(i,j) = -k(i+1);
            K(j,i) = -k(i+1);
        end
    end
end
C = diag(c)+diag([c(2:dofs);0]);
for i = 1:dofs-1
    for j = 1:dofs
        if j == i+1
            C(i,j) = -c(i+1);
            C(j,i) = -c(i+1);
        end
    end
end

[ph, wnsq] = eig(K,M); wn_htz = (diag(wnsq).^0.5)/(2*pi);

Ac = [zeros(dofs), eye(dofs); -inv(M)*K, -inv(M)*C]; Bc = [zeros(dofs); inv(M)]*(-M*g*ones(dofs,1));
Ad = expm(Ac*dt); Bd = (Ad-eye(2*dofs))/(Ac)*Bc;

xdv_lin(:,1) = zeros(2*dofs,1);
for i = 2:length(t)
    xdv_lin(:,i) = Ad*xdv_lin(:,i-1)+Bd*ug(i-1);
end

figure;
subplot(2,3,1); hold on; plot(t,xdv_lin(1,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(1,:),'r');
plot_properties({'Linear-Linear','BW-L'},'\textbf{Time (t)}','Disp.','Displacement Mass-1',1,1);
subplot(2,3,2); hold on; plot(t,xdv_lin(2,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(2,:),'r');
plot_properties('','\textbf{Time (t)}','Disp.','Displacement Mass-2',1,1);
subplot(2,3,3); hold on; plot(t,xdv_lin(3,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(3,:),'r');
plot_properties('','\textbf{Time (t)}','Disp.','Displacement Mass-3',1,1);
subplot(2,3,4); hold on; plot(t,xdv_lin(4,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(4,:),'r');
plot_properties('','\textbf{Time (t)}','Disp.','Displacement Mass-4',1,1);
subplot(2,3,5); hold on; plot(t,xdv_lin(5,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(5,:),'r');
plot_properties({'Linear-Linear','BW-L'},'\textbf{Time (t)}','Disp.','Displacement Mass-5',1,1);
subplot(2,3,6); hold on; plot(t,xdv_lin(6,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(6,:),'r');
plot_properties('','\textbf{Time (t)}','Vel.','Velocity Mass-1',1,1);
sub_title('Data Simulation (1/2)',0.515,0.96);

figure;
subplot(2,3,1); hold on; plot(t,xdv_lin(7,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(7,:),'r');
plot_properties('','\textbf{Time (t)}','Vel.','Velocity Mass-2',1,1);
subplot(2,3,2); hold on; plot(t,xdv_lin(8,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(8,:),'r');
plot_properties('','\textbf{Time (t)}','Vel.','Velocity Mass-3',1,1);
subplot(2,3,3); hold on; plot(t,xdv_lin(9,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(9,:),'r');
plot_properties('','\textbf{Time (t)}','Vel.','Velocity Mass-4',1,1);
subplot(2,3,4); hold on; plot(t,xdv_lin(10,:),'color',[0.5 0.5 0.5],'linewidth',1.5); plot(t,xdv(10,:),'r');
plot_properties('','\textbf{Time (t)}','Vel.','Velocity Mass-5',1,1);
subplot(2,3,[5,6]); hold on; plot(t,xdv(11,:),'r');
plot_properties('','\textbf{Time (t)}','z','Evolutionary Parameter',1,1);
sub_title('Data Simulation (2/2)',0.515,0.96);

%% Part III

Ac = [zeros(dofs), eye(dofs); -inv(M)*K, -inv(M)*C]; Bc = [zeros(dofs); inv(M)]*(-M*g*ones(dofs,1)); Cc = [zeros(dofs); -inv(M)];
Ad = expm(Ac*dt); Bd = (Ad-eye(2*dofs))/(Ac)*Bc; Cd = (Ad-eye(2*dofs))/(Ac)*Cc;

Am = [-inv(M)*K, -inv(M)*C]; Bm = (M)\(-M*g*ones(dofs,1)); Cm = -inv(M);

ps = zeros(2*dofs,1); pps = eye(2*dofs); pss = eye(2*dofs).*1e-10;
pi = zeros(dofs,1); ppi = eye(dofs); pii = 100*eye(dofs);

R = eye(dofs).*1e-10;
s = zeros(2*dofs,length(t)); nl = zeros(dofs,length(t));

for i = 1:length(t)
    ei = z(:,i) - Am*ps - Cm*pi;
    kgi = ((Cm*ppi*Cm'+R)\Cm*ppi)';
    
    ci = pi+kgi*ei;
    pci = ppi - kgi*Cm*ppi;
    
    es = z(:,i) - Am*ps - Cm*ci;
    kgs = ((Am*pps*Am'+R)\Am*pps)';
    
    cs = ps+kgs*es;
    pcs = pps - kgs*Am*pps;
    
    pi = ci;
    ppi = pci+pii;
    
    ps = Ad*cs+Bd*ug(i)+Cd*ci;
    pps = Ad*pcs*Ad'+pss;
    
    s(:,i) = cs; nl(:,i) = ci;
end

figure;
subplot(2,4,1); hold on; plot(t,xdv(1,:),'r'); plot(t,s(1,:),'b--');
plot_properties({'True','Filtered'},'\textbf{Time (t)}','Disp.','Displacement Mass-1',1,1);
subplot(2,4,2); hold on; plot(t,xdv(2,:),'r'); plot(t,s(2,:),'b--');
plot_properties('','\textbf{Time (t)}','Disp.','Displacement Mass-2',1,1);
subplot(2,4,3); hold on; plot(t,xdv(3,:),'r'); plot(t,s(3,:),'b--');
plot_properties('','\textbf{Time (t)}','Disp.','Displacement Mass-3',1,1);
subplot(2,4,4); hold on; plot(t,xdv(4,:),'r'); plot(t,s(4,:),'b--');
plot_properties('','\textbf{Time (t)}','Disp.','Displacement Mass-4',1,1);
subplot(2,4,5); hold on; plot(t,xdv(5,:),'r'); plot(t,s(5,:),'b--');
plot_properties({'True','Filtered'},'\textbf{Time (t)}','Disp.','Displacement Mass-5',1,1);
subplot(2,4,6); hold on; plot(t,xdv(6,:),'r'); plot(t,s(6,:),'b--');
plot_properties('','\textbf{Time (t)}','Vel.','Velocity Mass-1',1,1);
subplot(2,4,7); hold on; plot(t,xdv(7,:),'r'); plot(t,s(7,:),'b--');
plot_properties('','\textbf{Time (t)}','Vel.','Velocity Mass-2',1,1);
subplot(2,4,8); hold on; plot(t,xdv(8,:),'r'); plot(t,s(8,:),'b--');
plot_properties('','\textbf{Time (t)}','Vel.','Velocity Mass-3',1,1);
sub_title('Filteration (1/2)',0.515,0.96);

figure;
subplot(2,4,1); hold on; plot(t,xdv(9,:),'r'); plot(t,s(9,:),'b--');
plot_properties({'True','Filtered'},'\textbf{Time (t)}','Vel.','Velocity Mass-4',1,1);
subplot(2,4,2); hold on; plot(t,xdv(10,:),'r'); plot(t,s(10,:),'b--');
plot_properties('','\textbf{Time (t)}','Vel.','Velocity Mass-5',1,1);
subplot(2,4,[3,4]); plot(t,-(m1*xa(1,:)+c1*xdv(6,:)+c2*(xdv(6,:)-xdv(7,:))+k1*xdv(1,:)+k2*(xdv(1,:)-xdv(2,:))),'r');
hold on; plot(t,nl(1,:),'b--'); plot_properties('','\textbf{Time (t)}','$\mathbf{R_1}$','Non-linearity 1',1,1);
subplot(2,4,5); plot(t,-(m2*xa(2,:)+c2*(xdv(7,:)-xdv(6,:))+c3*(xdv(7,:)-xdv(8,:))+k2*(xdv(2,:)-xdv(1,:))+k3*(xdv(2,:)-xdv(3,:))),'r');
hold on; plot(t,nl(2,:),'b--'); plot_properties('','\textbf{Time (t)}','$\mathbf{R_2}$','Non-linearity 2',1,1);
subplot(2,4,6); plot(t,-(m3*xa(3,:)+c3*(xdv(8,:)-xdv(7,:))+c4*(xdv(8,:)-xdv(9,:))+k3*(xdv(3,:)-xdv(2,:))+k4*(xdv(3,:)-xdv(4,:))),'r');
hold on; plot(t,nl(3,:),'b--'); plot_properties({'True','Filtered'},'\textbf{Time (t)}','$\mathbf{R_3}$','Non-linearity 3',1,1);
subplot(2,4,7); plot(t,-(m4*xa(4,:)+c4*(xdv(9,:)-xdv(8,:))+c5*(xdv(9,:)-xdv(10,:))+k4*(xdv(4,:)-xdv(3,:))+k5*(xdv(4,:)-xdv(5,:))),'r');
hold on; plot(t,nl(4,:),'b--'); plot_properties('','\textbf{Time (t)}','$\mathbf{R_4}$','Non-linearity 4',1,1);
subplot(2,4,8); plot(t,-(m5*xa(5,:)+c5*(xdv(10,:)-xdv(9,:))+k5*(xdv(5,:)-xdv(4,:))),'r');
hold on; plot(t,nl(5,:),'b--'); plot_properties('','\textbf{Time (t)}','$\mathbf{R_5}$','Non-linearity 5',1,1);
sub_title('Filteration (2/2)',0.515,0.96);

%% Part IV

tr_lim = 30; ts = dt;

% parfor i = 1:dofs
for i = 1:dofs
    switch i
        case 1
            disp(['model_start ',num2str(i)]);
            mdl{i}  = fitrgp(s([1,2,6,7],1:tr_lim/ts)',nl(1,1:tr_lim/ts)',... 
                'OptimizeHyperparameters','none',...
                'HyperparameterOptimizationOptions',struct('UseParallel',true));%#ok<*PFBNS>
            a{i} = predict(mdl{i},s([1,2,6,7],:)');
            disp(['model_complete ',num2str(i)]);
        case 2
            disp(['model_start ',num2str(i)]);
            mdl{i}  = fitrgp(s([1,2,3,6,7,8],1:tr_lim/ts)',nl(2,1:tr_lim/ts)',...
                'OptimizeHyperparameters','none',...
                'HyperparameterOptimizationOptions',struct('UseParallel',true));
            a{i} = predict(mdl{i},s([1,2,3,6,7,8],:)');
            disp(['model_complete ',num2str(i)]);
        case 3
            disp(['model_start ',num2str(i)]);
            mdl{i}  = fitrgp(s([2,3,4,7,8,9],1:tr_lim/ts)',nl(3,1:tr_lim/ts)',...
                'OptimizeHyperparameters','none',...
                'HyperparameterOptimizationOptions',struct('UseParallel',true));
            a{i} = predict(mdl{i},s([2,3,4,7,8,9],:)');
            disp(['model_complete ',num2str(i)]);
        case 4
            disp(['model_start ',num2str(i)]);
            mdl{i}  = fitrgp(s([3,4,5,8,9,10],1:tr_lim/ts)',nl(4,1:tr_lim/ts)',...
                'OptimizeHyperparameters','none',...
                'HyperparameterOptimizationOptions',struct('UseParallel',true));
            a{i} = predict(mdl{i},s([3,4,5,8,9,10],:)');
            disp(['model_complete ',num2str(i)]);
        otherwise
            disp(['model_start ',num2str(i)]);
            mdl{i}  = fitrgp(s([4,5,9,10],1:tr_lim/ts)',nl(5,1:tr_lim/ts)',...
                'OptimizeHyperparameters','none',...
                'HyperparameterOptimizationOptions',struct('UseParallel',true));
            a{i} = predict(mdl{i},s([4,5,9,10],:)');
            disp(['model_complete ',num2str(i)]);
    end
end

a1 = a{1}; a2 = a{2}; a3 = a{3}; a4 = a{4}; a5 = a{5};
mdl1 = mdl{1}; mdl2 = mdl{2}; mdl3 = mdl{3}; mdl4 = mdl{4}; mdl5 = mdl{5};

figure;
subplot(2,3,1); plot(t,-(m1*xa(1,:)+c1*xdv(6,:)+c2*(xdv(6,:)-xdv(7,:))+k1*xdv(1,:)+k2*(xdv(1,:)-xdv(2,:))),'b--','linewidth',1.5);
hold on; plot(t,a1,'m','linewidth',1.5);
xline(t(tr_lim/dt),'linewidth',2.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{R_1}$','',1,1);

subplot(2,3,2); plot(t,-(m2*xa(2,:)+c2*(xdv(7,:)-xdv(6,:))+c3*(xdv(7,:)-xdv(8,:))+k2*(xdv(2,:)-xdv(1,:))+k3*(xdv(2,:)-xdv(3,:))),'b--','linewidth',1.5);
hold on; plot(t,a2,'m','linewidth',1.5); 
xline(t(tr_lim/dt),'linewidth',2.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{R_2}$','',1,1);

subplot(2,3,3); plot(t,-(m3*xa(3,:)+c3*(xdv(8,:)-xdv(7,:))+c4*(xdv(8,:)-xdv(9,:))+k3*(xdv(3,:)-xdv(2,:))+k4*(xdv(3,:)-xdv(4,:))),'b--','linewidth',1.5);
hold on; plot(t,a3,'m','linewidth',1.5); 
xline(t(tr_lim/dt),'linewidth',2.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{R_3}$','',1,1);

subplot(2,3,4); plot(t,-(m4*xa(4,:)+c4*(xdv(9,:)-xdv(8,:))+c5*(xdv(9,:)-xdv(10,:))+k4*(xdv(4,:)-xdv(3,:))+k5*(xdv(4,:)-xdv(5,:))),'b--','linewidth',1.5);
hold on; plot(t,a4,'m','linewidth',1.5); 
xline(t(tr_lim/dt),'linewidth',2.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{R_4}$','',1,1);

subplot(2,3,[5,6]); plot(t,-(m5*xa(5,:)+c5*(xdv(10,:)-xdv(9,:))+k5*(xdv(5,:)-xdv(4,:))),'b--','linewidth',1.5);
hold on; plot(t,a5,'m','linewidth',1.5); 
xline(t(tr_lim/dt),'linewidth',2.5); xlim([0 t(end)]);
plot_properties({'Ground Truth','GP Estimates','Last Training Data'},'\textbf{Time (t)}','$\mathbf{R_5}$','',1,1);

%% Part V

xdv_s_gp(:,1) = zeros(2*dofs,1);
for i = 2:length(t)
    disp(['step ',num2str(i)]);

    av(1,:) = predict(mdl1,[xdv_s_gp(1,i-1), xdv_s_gp(2,i-1), xdv_s_gp(6,i-1), xdv_s_gp(7,i-1)]);
    for j = 2:dofs-1
        av(j,:) = predict(mdl{j},[xdv_s_gp(j-1,i-1), xdv_s_gp(j,i-1), xdv_s_gp(j+1,i-1), xdv_s_gp(j+4,i-1), xdv_s_gp(j+5,i-1), xdv_s_gp(j+6,i-1)]);
    end
    av(5,:) = predict(mdl5,[xdv_s_gp(4,i-1), xdv_s_gp(5,i-1), xdv_s_gp(9,i-1), xdv_s_gp(10,i-1)]);
    
    xdv_s_gp(:,i) = Ad*xdv_s_gp(:,i-1)+Bd*ug(i-1)+Cd*av;
end

figure;
subplot(2,1,1); hold on; plot(t,xdv(1,:),'r','linewidth',1.5);
plot(t,xdv_s_gp(1,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties({'Ground Truth','Projected States'},'\textbf{Time (t)}','$\mathbf{x_1}$','',1,0.5);

subplot(2,1,2); hold on; plot(t,xdv(2,:),'r','linewidth',1.5);
plot(t,xdv_s_gp(2,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{x_2}$','',1,0.5);

figure;
subplot(3,1,1); hold on; plot(t,xdv(3,:),'r','linewidth',1.5);
plot(t,xdv_s_gp(3,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{x_3}$','',1,0.5);

subplot(3,1,2); hold on; plot(t,xdv(4,:),'r','linewidth',1.5);
plot(t,xdv_s_gp(4,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{x_4}$','',1,0.5);

subplot(3,1,3); hold on; plot(t,xdv(5,:),'r','linewidth',1.5);
plot(t,xdv_s_gp(5,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{x_5}$','',1,0.5);

%% Part VI

load elcentro.txt; eq  = elcentro;
dt_new = 0.001; t = eq(:,1); ug = eq(:,2); dt = t(2)-t(1);
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
    f_ug(i,:) = f;
end

ug = 0.9*f_ug(3,:); figure; plot(t,ug);
plot_properties('','\textbf{Time (t)}','Gr. Acc.','Ground Acceleration - Testing',1,1);

xdv_d_gp(:,1) = zeros(2*dofs,1);
for i = 2:length(t)
    disp(['step ',num2str(i)]);
    
    av(1,:) = predict(mdl1,[xdv_d_gp(1,i-1), xdv_d_gp(2,i-1), xdv_d_gp(6,i-1), xdv_d_gp(7,i-1)]);
    for j = 2:dofs-1
        av(j,:) = predict(mdl{j},[xdv_d_gp(j-1,i-1), xdv_d_gp(j,i-1), xdv_d_gp(j+1,i-1), xdv_d_gp(j+4,i-1), xdv_d_gp(j+5,i-1), xdv_d_gp(j+6,i-1)]);
    end
    av(5,:) = predict(mdl5,[xdv_d_gp(4,i-1), xdv_d_gp(5,i-1), xdv_d_gp(9,i-1), xdv_d_gp(10,i-1)]);
    
    xdv_d_gp(:,i) = Ad*xdv_d_gp(:,i-1)+Bd*ug(i-1)+Cd*av;
end

m = 20*[20; 19; 18; 17; 16];
c = 10*[10; 20; 19; 18; 17];
k = 10000*[10; 20; 19; 18; 17];

m1 = m(1); m2 = m(2); m3 = m(3); m4 = m(4); m5 = m(5);
c1 = c(1); c2 = c(2); c3 = c(3); c4 = c(4); c5 = c(5);
k1 = k(1); k2 = k(2); k3 = k(3); k4 = k(4); k5 = k(5);

g = 9.81; Qy = 0.05*(m1+m2+m3+m4+m5)*g; alpha = 1;
beta = 0.5; gamma = 0.5; eta = 1; Dy = 0.013; kr = 1/6; dofs = 5;

xdv_d_tr(:,1) = [zeros(2*dofs,1); 1e-5];
for i = 2:length(t)
    y = xdv_d_tr(:,i-1);
    xdv_d_tr(:,i) = disp_vel_IT_5(y, ug(i-1), f_ug(1,i-1), f_ug(2,i-1), m, c, k);
end

figure;
subplot(2,1,1); hold on; plot(t,xdv_d_tr(1,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(1,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties({'Ground Truth','Projected States'},'\textbf{Time (t)}','$\mathbf{x_1}$','',1,1);

subplot(2,1,2); hold on; plot(t,xdv_d_tr(2,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(2,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{x_2}$','',1,1);

figure;
subplot(3,1,1); hold on; plot(t,xdv_d_tr(3,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(3,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{x_3}$','',1,1);

subplot(3,1,2); hold on; plot(t,xdv_d_tr(4,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(4,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{x_4}$','',1,1);

subplot(3,1,3); hold on; plot(t,xdv_d_tr(5,:),'r','linewidth',1.5);
plot(t,xdv_d_gp(5,:),'b--','linewidth',1.5); xlim([0 t(end)]);
plot_properties('','\textbf{Time (t)}','$\mathbf{x_5}$','',1,1);