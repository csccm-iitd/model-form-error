clc
clear
close all

global dt alpha_do alpha_dvp

dt = 0.001; t = 0:dt:50;
ug = sin(10*t)*9.81;

%%

m1 = 10; c1 = 2.5; k1 = 100;
f1 = -m1*ug; alpha_dvp = 10;

xdv(:,1) = [0; 0];
for i = 2:length(t)
    y1 = xdv(1,i-1); y2 = xdv(2,i-1);
    y = [y1; y2]; r = [randn; randn];
    xdv(:,i) = disp_vel_IT(y, f1(i-1), r, k1, c1, m1);
end

figure;
subplot(2,1,1); hold on; plot(t,xdv(1,:),'r','linewidth',1.5);
plot_properties({'Ground Truth'},'\textbf{Time (t)}','$\mathbf{x}$','',1,1);
subplot(2,1,2); hold on; plot(t,xdv(2,:),'r','linewidth',1.5);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}}$','',1,1);

xa(:,1) = 0;
for i = 1:length(t)
    y1 = xdv(1,i); y2 = xdv(2,i);
    
    xa(:,i) = (- alpha_dvp*y1^3 + k1*y1 - c1*y2)/m1;
end

snr = 1e10;
z(1,:) = xdv(1,:)+sqrt((1/snr)*(std(xdv(1,:))^2))*randn(1,length(xdv(1,:)));
z(2,:) = xa(1,:)+sqrt((1/snr)*(std(xa(1,:))^2))*randn(1,length(xa(1,:)));

%% Part III

alpha_do = 11; k1 = 50;
 
M = m1; Cc = [0; -inv(M)];
Cd = dt*Cc; Cm = -inv(M);

m02 = 0; n2 = size(m02,1); al2 = 1; beta2 = 2; kappa2 = 0;
la2 = al2^2*(n2+kappa2)-n2; wm2 = zeros(2*n2+1,1); wc2 = zeros(2*n2+1,1);
for j=1:2*n2+1
    if j==1
        wm2(j) = la2/(n2+la2);
        wc2(j) = la2/(n2+la2)+(1-al2^2+beta2);
    else
        wm2(j) = 1/(2*(n2+la2));
        wc2(j) = wm2(j);
    end
end
p02 = eye(1); mk2 = m02; pk2 = p02;

m03 = zeros(2,1); n3 = size(m03,1); al3 = 1; beta3 = 2; kappa3 = 0;
la3 = al3^2*(n3+kappa3)-n3; wm3 = zeros(2*n3+1,1); wc3 = zeros(2*n3+1,1);
for j=1:2*n3+1
    if j==1
        wm3(j) = la3/(n3+la3);
        wc3(j) = la3/(n3+la3)+(1-al3^2+beta3);
    else
        wm3(j) = 1/(2*(n3+la3));
        wc3(j) = wm3(j);
    end
end
p03 = eye(2); mk3 = m03; pk3 = p03;
nl = 0; s = [0; 0];

dofs = 1; R = cov(z');

for i = 1:length(z)
    mkm = mk2; pkm = pk2;
    xhikm = [mkm mkm+sqrt(n2+la2)*chol(pkm,'lower') mkm-sqrt(n2+la2)*chol(pkm,'lower')];
    zk = zeros(dofs+1, length(xhikm));
    for j=1:length(xhikm)
        y1 = mk3(1); y2 = mk3(2);
        
        misc = -(alpha_do*y1^3 + k1*y1 + c1*y2)/m1;

        zk(:,j) = [y1; (misc + Cm*xhikm(:,j))];
        
    end
    muk = 0;
    for j=1:length(xhikm)
        muk = muk + wm2(j)*zk(:,j);
    end
    sk = 0;
    for j=1:length(xhikm)
        sk = sk + wc2(j)*(zk(:,j)-muk)*(zk(:,j)-muk)';
    end
    sk = sk+R;
    ck = 0;
    for j=1:length(xhikm)
        ck = ck + wc2(j)*(xhikm(:,j)-mkm)*(zk(:,j)-muk)';
    end
    kk = ck*inv(sk); %#ok<*MINV>
    mk2 = mkm + kk*(z(:,i)-muk);
    pk2 = pkm - kk*sk*kk';
    nl(:,i) = mk2;
    clear mkm pkm xhikm zk sk ck kk qc Q
        
    %%%%%%%%%%
    
    mkm = mk3; pkm = pk3;
    xhikm = [mkm mkm+sqrt(n3+la3)*chol(pkm,'lower') mkm-sqrt(n3+la3)*chol(pkm,'lower')];
    zk = zeros(dofs+1, length(xhikm));
    for j=1:length(xhikm)
        y1 = xhikm(1,j); y2 = xhikm(2,j);
        
        misc = -(alpha_do*y1^3 + k1*y1 + c1*y2)/m1;

        zk(:,j) = [y1; (misc + Cm*mk2)];
    end
    muk = 0;
    for j=1:length(xhikm)
        muk = muk + wm3(j)*zk(:,j);
    end
    sk = 0;
    for j=1:length(xhikm)
        sk = sk + wc3(j)*(zk(:,j)-muk)*(zk(:,j)-muk)';
    end
    sk = sk + R;
    ck = 0;
    for j=1:length(xhikm)
        ck = ck + wc3(j)*(xhikm(:,j)-mkm)*(zk(:,j)-muk)';
    end
    kk = ck*inv(sk); %#ok<*MINV>
    mk3 = mkm + kk*(z(:,i)-muk);
    pk3 = pkm - kk*sk*kk';
    s(:,i) = mk3;
    clear mkm pkm xhikm zk sk ck kk qc Q
    
    %%%%%%%%%%
    
%     qc = 1e2; Q = qc*qc';
%     pk2 = pk2+Q;    
    
    yhikmi = [mk2 mk2+sqrt(n2+la2)*chol(pk2,'lower') mk2-sqrt(n2+la2)*chol(pk2,'lower')];
    yhik = zeros(size(yhikmi));
    for j=1:length(yhikmi)
        yhik(:,j) = yhikmi(:,j);
    end
    mkm = 0;
    for j=1:length(yhik)
        mkm = mkm + wm2(j)*yhik(:,j);
    end
    pkm = 0;
    for j=1:length(yhik)
        pkm = pkm + wc2(j)*((yhik(:,j)-mkm)*(yhik(:,j)-mkm)');
    end
    qc = 1e2; Q = qc*qc';
    pkm = pkm+Q;
    mk2 = mkm; pk2 = pkm;
    
    %%%%%%%%%%
    
    fprintf('step (input-state) %d\n',i);
    yhikmi = [mk3 mk3+sqrt(n3+la3)*chol(pk3,'lower') mk3-sqrt(n3+la3)*chol(pk3,'lower')];
    yhik = zeros(size(yhikmi));
    for j=1:length(yhikmi)
        y1 = yhikmi(1,j); y2 = yhikmi(2,j);
        
        misc = disp_vel_filter(y1, y2, f1(i), k1, c1, m1,t(i));
        
        yhik(:,j) = misc+Cd*mk2;
    end
    mkm = 0;
    for j=1:length(yhik)
        mkm = mkm + wm3(j)*yhik(:,j);
    end
    pkm = 0;
    for j=1:length(yhik)
        pkm = pkm + wc3(j)*((yhik(:,j)-mkm)*(yhik(:,j)-mkm)');
    end
    qc = 1e-3.*diag([0 1]); Q = qc*qc';
    pkm = pkm+Q;
    mk3 = mkm; pk3 = pkm;
end

figure;
subplot(2,2,1); plot(t,xdv(1,:),'r','linewidth',1.5); hold on; plot(t,s(1,:),'b--','linewidth',1.5);
plot_properties({'Ground Truth','Filtered Estimate'},'\textbf{Time (t)}','$\mathbf{x}$','',1,1);
subplot(2,2,2); plot(t,xdv(2,:),'r','linewidth',1.5); hold on; plot(t,s(2,:),'b--','linewidth',1.5);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}}$','',1,1);
subplot(2,2,[3,4]); plot(t,-(m1*xa(1,:)+alpha_do*xdv(1,:).^3+k1*xdv(1,:)+c1*xdv(2,:)),'r','linewidth',1.5); hold on; plot(t,nl(1,:),'b--','linewidth',1.5);
plot_properties('','\textbf{Time (t)}','R','',1,1);

%%

tr_lim = 20;

mdl  = fitrgp(s(:,1:tr_lim/dt)',nl(1,1:tr_lim/dt)','PredictMethod','sr');
[pr,~,ci] = predict(mdl,s');

figure; plot(t,-(m1*xa(1,:)+alpha_do*xdv(1,:).^3+k1*xdv(1,:)+c1*xdv(2,:)),'b','linewidth',1.5);
hold on; plot(t,pr,'m--','linewidth',1.5); xline(t(tr_lim/dt),'linewidth',1.5);
plot_properties({'Ground Truth','GP Predictions','Last Training Data'},'\textbf{Time (t)}','$\mathbf{R}$','',1,0.25);

%%

gp_a = 25000;
gp_b = 27500;

figure; plot(t,-(m1*xa(1,:)+alpha_do*xdv(1,:).^3+k1*xdv(1,:)+c1*xdv(2,:)),'b','linewidth',1.5);
hold on; plot(t,pr,'m--','linewidth',1.5); xline(t(tr_lim/dt),'linewidth',1.5);
hold on; patch([t';flipud(t')],[ci(:,1);flipud(ci(:,2))],'m','FaceAlpha',0.1);
plot_properties({'Ground Truth','GP Predictions','CI','Last Training Data'},'\textbf{Time (t)}','$\mathbf{R}$','',1,0.25);

figure; plot(t(gp_a : gp_b),-(m1*xa(1,(gp_a : gp_b))+alpha_do*xdv(1,(gp_a : gp_b)).^3+k1*xdv(1,(gp_a : gp_b))+c1*xdv(2,(gp_a : gp_b))),'b','linewidth',1.5);
hold on; plot(t(gp_a : gp_b),pr(gp_a : gp_b),'m--','linewidth',1.5);
hold on; patch([t(gp_a : gp_b)';flipud(t(gp_a : gp_b)')],[ci((gp_a : gp_b),1);flipud(ci((gp_a : gp_b),2))],'m','FaceAlpha',0.1);
plot_properties({'Ground Truth','GP Predictions','CI'},'\textbf{Time (t)}','$\mathbf{R}$','',1,0.25);

%%
 
xdv_s_gp(:,1) = [0; 0];
for i = 2:length(t)
    disp(['step ',num2str(i)]);
    av = predict(mdl,[xdv_s_gp(1,i-1), xdv_s_gp(2,i-1)]);

    y1 = xdv_s_gp(1,i-1); y2 = xdv_s_gp(2,i-1);
    misc = disp_vel_filter(y1, y2, f1(i-1), k1, c1, m1,t(i));
        
    xdv_s_gp(:,i) = misc+Cd*av;
end

figure;
subplot(2,1,1); hold on; plot(t,xdv(1,:),'r'); plot(t,xdv_s_gp(1,:),'b--','linewidth',1.5);
plot_properties({'Ground Truth','Projected States'},'\textbf{Time (t)}','$\mathbf{x}$','',1,1);
subplot(2,1,2); hold on; plot(t,xdv(2,:),'r'); plot(t,xdv_s_gp(2,:),'b--','linewidth',1.5);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}}$','',1,1);

%%
 
ug = sin(15*t)*9.81;
f1 = -m1*ug;

xdv_d_gp(:,1) = [0; 0];
for i = 2:length(t)
    disp(['step ',num2str(i)]);
    
    av = predict(mdl,[xdv_d_gp(1,i-1), xdv_d_gp(2,i-1)]);

    y1 = xdv_d_gp(1,i-1); y2 = xdv_d_gp(2,i-1);
    misc = disp_vel_filter(y1, y2, f1(i-1), k1, c1, m1,t(i));
        
    xdv_d_gp(:,i) = misc+Cd*av;
end

m1 = 10; c1 = 2.5; k1 = 100;
f1 = -m1*ug; alpha_dvp = 10;

xdv_d_tr(:,1) = [0; 0];
for i = 2:length(t)
    y1 = xdv_d_tr(1,i-1); y2 = xdv_d_tr(2,i-1);
    y = [y1; y2]; r = [randn; randn];
    xdv_d_tr(:,i) = disp_vel_IT(y, f1(i-1), r, k1, c1, m1);
end

figure;
subplot(2,1,1); hold on; plot(t,xdv_d_tr(1,:),'r'); plot(t,xdv_d_gp(1,:),'b--','linewidth',1.5);
plot_properties({'Ground Truth','Projected States'},'\textbf{Time (t)}','$\mathbf{x}$','',1,1);
subplot(2,1,2); hold on; plot(t,xdv_d_tr(2,:),'r'); plot(t,xdv_d_gp(2,:),'b--','linewidth',1.5);
plot_properties('','\textbf{Time (t)}','$\mathbf{\dot{x}}$','',1,1);

%% Functions

function y = disp_vel_filter(y1, y2, f1, k1, c1, m1,t)
global dt alpha_do

y = [y1; y2];

a = [y2;
    -(alpha_do*y1^3 + k1*y1 - f1 + c1*y2)/m1];

y = y+a*dt;
end

function y = disp_vel_IT(y, f1, r, k1, c1, m1)
global dt alpha_dvp

sigma1 = 0.1; r1 = r(1); r2 = r(2);

y1 = y(1); y2 = y(2);
y = [y1; y2];

a = [y2;
    (- alpha_dvp*y1^3 + k1*y1 + f1 - c1*y2)/m1];

b = [0;
    (sigma1*y1)/m1];

L0a = [(- alpha_dvp*y1^3 + k1*y1 + f1 - c1*y2)/m1;
    (y2*(- 3*alpha_dvp*y1^2 + k1))/m1 - (c1*(- alpha_dvp*y1^3 + k1*y1 + f1 - c1*y2))/m1^2];

L0b = [0;
    (sigma1*y2)/m1];

L1a = [(sigma1*y1)/m1;
    -(c1*sigma1*y1)/m1^2];

dz = (dt^(3/2)*r1)/2 + (3^(1/2)*dt^(3/2)*r2)/6;

dw = dt^(1/2)*r1;

y = y+a*dt+b*dw+L0b*(dw*dt-dz)+L1a*dz+0.5*L0a*dt^2;
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

