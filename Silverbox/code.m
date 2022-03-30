clc
clear
close all

load Schroeder80mV.mat

ll = 1; ul = length(V1);
sf = 610.35; dt = 1/sf; samples = ll:ul;
Time_total = samples*dt;

figure;
subplot(2,1,1); plot(Time_total,V1,'color',[0.5,0.5,0.5],'linewidth',1.5);
subplot(2,1,2); plot(Time_total,V2,'color',[0.5,0.5,0.5],'linewidth',1.5);

nt = 5000/2; ll = 1000;
ul = nt+ll-1; samples = ll:ul;
t = samples*dt;

subplot(2,1,1); hold on; plot(t,V1(samples),'r','linewidth',1.5); 
xline(t(1),'r--','linewidth',1.5); xline(t(end),'r--','linewidth',1.5);
subplot(2,1,2); hold on; plot(t,V2(samples),'r','linewidth',1.5);
xline(t(1),'r--','linewidth',1.5); xline(t(end),'r--','linewidth',1.5);

ns = 5000; ll = 10500;
ul = ns+ll-1; samples = ll:ul;

t = samples*dt;

subplot(2,1,1); hold on; plot(t,V1(samples),'b','linewidth',1.5);
xline(t(1),'b-.','linewidth',1.5); xline(t(end),'b-.','linewidth',1.5);
plot_properties({'Force','','GP Training Window','','','Testing Window'},'\textbf{Time (t)}','$\mathbf{f}$','',1,0.5);
xlim([Time_total(1), Time_total(end)]);

subplot(2,1,2); plot(t,V2(samples),'b','linewidth',1.5);
xline(t(1),'b-.','linewidth',1.5); xline(t(end),'b-.','linewidth',1.5);
plot_properties({'Displacement','','GP Training Window','','','Testing Window'},'\textbf{Time (t)}','$\mathbf{x}$','',1,0.5);
xlim([Time_total(1), Time_total(end)]);

%%

V1 = V1-mean(V1);
V2 = V2-mean(V2);

ns = 5000; ll = 1000;
ul = ns+ll-1; samples = ll:ul;

sf = 610.35;
dt = 1/sf; t = samples*dt;

force = V1(ll:ul);
displacement = V2(ll:ul);

velocity = [diff(displacement)./diff(t) 0];
acceleration = [diff(velocity)./diff(t) 0];

figure; subplot(3,2,1); plot(t,force,'r'); title('Force'); xlim([t(1), t(end)]);
subplot(3,2,2); plot(t, displacement,'b'); title('Displacement'); xlim([t(1), t(end)]);

amp = 2e5;
M = 1; K = 50000; C = 2.5; force = amp*force;

Ac = [0, 1; -inv(M)*K, -inv(M)*C];
Bc = [0; inv(M)]; Cc = [zeros(1); -inv(M)];

Ad = expm(Ac*dt); Bd = (Ad-eye(2))*inv(Ac)*Bc;
Cd = (Ad-eye(2))*inv(Ac)*Cc;

z = acceleration;
Am = [-inv(M)*K, -inv(M)*C]; Bm = inv(M); Cm = -inv(M);

ps = [0; 0]; pps = diag([1 1]); pss = eye(2).*1e-4;
pi = 0; ppi = 1; pii = 100; R = 1.*1e-10;

s = zeros(2,length(t)); nl = zeros(1,length(t));

for i = 1:length(t)
    ei = z(:,i) - (Am*ps + Bm*force(i) + Cm*pi);
    kgi = (inv(Cm*ppi*Cm'+R)*Cm*ppi)';
    
    ci = pi+kgi*ei;
    pci = ppi - kgi*Cm*ppi;
    
    es = z(:,i) - (Am*ps + Bm*force(i) + Cm*ci);
    kgs = (inv(Am*pps*Am'+R)*Am*pps)';
    
    cs = ps+kgs*es;
    pcs = pps - kgs*Am*pps;
    
    pi = ci;
    ppi = pci+pii;
    
    ps = Ad*cs+Bd*force(i)+Cd*ci;
    pps = Ad*pcs*Ad'+pss;
    
    s(:,i) = cs; nl(:,i) = ci;
end

subplot(3,2,[3, 4]); plot(t,s(1,:)-mean(s(1,:)),'b',t,displacement,'r')
title('Displacement Estimate'); xlim([t(1), t(end)]);
legend('disp\_est-mean(disp\_est)','ground truth',fontsize = 16);

subplot(3,2,5); plot(t,s(2,:),'b',t,velocity,'r')
title('Velocity Estimate'); xlim([t(1), t(end)]);
legend('velocity\_estimate','ground truth',fontsize = 16);

subplot(3,2,6); plot(t,nl(1,:),'k')
title('Residual Force'); xlim([t(1), t(end)]);

%%

s(1,:) = s(1,:)-mean(s(1,:)); % Displacement Estimate
s(2,:) = s(2,:)-mean(s(2,:)); % Velocity Estimate

nl = nl-mean(nl); % Residual Force Estimate

nte = [750, 1500, 2500];
nt = 750;

mdl  = fitrgp(s(:,1:nt)',nl(1,1:nt)','Basis','pureQuadratic',...
    'KernelFunction','ardmatern32',...
    'FitMethod','sr','PredictMethod','exact','Standardize',1);

pr = predict(mdl,s');

figure; plot(t,nl(1,:),'b:','linewidth',1.5);
hold on; plot(t,pr,'m','linewidth',1.5); xlim([t(1), t(end)]);
xline(t(nt),'k--','linewidth',2.5);
plot_properties({'Estimated (R)','GP Predictions','Last Training Data'},...
    '\textbf{Time (t)}','$\mathbf{R}$','',0.75,0.25);

rmse(1) = sqrt(mean((nl(1,:)-pr').^2));

%%

s(1,:) = s(1,:)-mean(s(1,:)); % Displacement Estimate
s(2,:) = s(2,:)-mean(s(2,:)); % Velocity Estimate

nl = nl-mean(nl); % Residual Force Estimate

nt = 1500;

mdl  = fitrgp(s(:,1:nt)',nl(1,1:nt)','Basis','pureQuadratic',...
    'KernelFunction','ardmatern32',...
    'FitMethod','sr','PredictMethod','exact','Standardize',1);

pr = predict(mdl,s');

figure; plot(t,nl(1,:),'b:','linewidth',1.5);
hold on; plot(t,pr,'m','linewidth',1.5); xlim([t(1), t(end)]);
xline(t(nt),'k--','linewidth',2.5);
plot_properties({'Estimated (R)','GP Predictions','Last Training Data'},...
    '\textbf{Time (t)}','$\mathbf{R}$','',0.75,0.25);

rmse(2) = sqrt(mean((nl(1,:)-pr').^2));

%%

s(1,:) = s(1,:)-mean(s(1,:)); % Displacement Estimate
s(2,:) = s(2,:)-mean(s(2,:)); % Velocity Estimate

nl = nl-mean(nl); % Residual Force Estimate

nt = 2500;

mdl  = fitrgp(s(:,1:nt)',nl(1,1:nt)','Basis','pureQuadratic',...
    'KernelFunction','ardmatern32',...
    'FitMethod','sr','PredictMethod','exact','Standardize',1);

pr = predict(mdl,s');

figure; plot(t,nl(1,:),'b:','linewidth',1.5);
hold on; plot(t,pr,'m','linewidth',1.5); xlim([t(1), t(end)]);
xline(t(nt),'k--','linewidth',2.5);
plot_properties({'Estimated (R)','GP Predictions','Last Training Data'},...
    '\textbf{Time (t)}','$\mathbf{R}$','',0.75,0.25);

rmse(3) = sqrt(mean((nl(1,:)-pr').^2));

%%

ns = 5000; ll = 10500;
ul = ns+ll-1; samples = ll:ul;

t = samples*dt;

figure;
plot(t(nte)-t(1), rmse/max(abs(nl(1,:))), 'r--o', linewidth = 2);
ylim padded
plot_properties('','\textbf{Training length}','\textbf{NRMSE}','',1,0.5);

%%

xdv_s_gp(:,1) = [0; 0];
for i = 2:length(t)
    av = predict(mdl,[xdv_s_gp(1,i-1), xdv_s_gp(2,i-1)]);

    xdv_s_gp(:,i) = Ad*[xdv_s_gp(1,i-1); xdv_s_gp(2,i-1)] +...
        Bd*force(i-1) + Cd*av;
end

figure;
subplot(2,1,1); plot(t(4000:4500),displacement(4000:4500),'r:','linewidth',1.5);
hold on; plot(t(4000:4500),xdv_s_gp(1,4000:4500),'b','linewidth',1.5);
plot_properties({'Ground truth','Estimated'},'\textbf{Time (t)}','$\mathbf{x}$','',1,0.5);

subplot(2,1,2); plot(t(4000:4500),velocity(4000:4500),'r:','linewidth',1.5);
hold on; plot(t(4000:4500),xdv_s_gp(2,4000:4500),'b','linewidth',1.5);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}}$','',1,0.5);

%%

ns = 5000; ll = 10500;
ul = ns+ll-1; samples = ll:ul;

t = samples*dt;
force = amp*V1(ll:ul);
displacement = V2(ll:ul);

velocity = [diff(displacement)./diff(t) 0];

xdv_d_gp(:,1) = [0; 0];
for i = 2:length(t)
    av = predict(mdl,[xdv_d_gp(1,i-1), xdv_d_gp(2,i-1)]);

    xdv_d_gp(:,i) = Ad*[xdv_d_gp(1,i-1); xdv_d_gp(2,i-1)] +...
        Bd*force(i-1) + Cd*av;
end

xdv_lin(:,1) = [0; 0];
for i = 2:length(t)
    xdv_lin(:,i) = Ad*[xdv_lin(1,i-1); xdv_lin(2,i-1)] + Bd*force(i-1);
end

figure;
subplot(2,1,1); plot(t,xdv_lin(1,:),'-.','color',[0.5 0.5 0.5],'linewidth',1.5);
hold on; plot(t,displacement,'r:','linewidth',1.5);
plot(t,xdv_d_gp(1,:),'b','linewidth',1.5);
plot_properties({'Linear System','Ground truth','Estimated'},'\textbf{Time (t)}','$\mathbf{x}$','',1,0.5);

subplot(2,1,2); plot(t,xdv_lin(2,:),'-.','color',[0.5 0.5 0.5],'linewidth',1.5);
hold on; plot(t,velocity,'r:','linewidth',1.5);
plot(t,xdv_d_gp(2,:),'b','linewidth',1.5);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}}$','',1,0.5);

%%

function plot_properties(legends,x_label,y_label,title_fig,index_lt,index_ht,legend_font,fig_font)

if nargin == 6
    legend_font = 16;
    fig_font = 22;
end
if nargin == 7
    fig_font = 22;
end

if isempty(legends)
    set(gca,'LineWidth',2,'FontSize',fig_font,'FontWeight','bold','FontName','Times');
    set(gcf,'Position',[1 1 index_lt*round(2000) index_ht*round(2000)]);
    set(get(gca,'xlabel'),'String',x_label,'FontSize',fig_font','FontWeight','bold','FontName','Times','Interpreter','Latex');
    set(get(gca,'ylabel'),'String',y_label,'FontSize',fig_font','FontWeight','bold','FontName','Times','Interpreter','Latex');
    set(gcf,'color','w'); box on; title(title_fig);
else
    legend(legends,'FontSize',legend_font,'FontWeight','bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','bold','FontName','Times');
    set(gcf,'Position',[1 1 index_lt*round(2000) index_ht*round(2000)]);
    set(get(gca,'xlabel'),'String',x_label,'FontSize',fig_font','FontWeight','bold','FontName','Times','Interpreter','Latex');
    set(get(gca,'ylabel'),'String',y_label,'FontSize',fig_font','FontWeight','bold','FontName','Times','Interpreter','Latex');
    set(gcf,'color','w'); box on; title(title_fig);
end
end