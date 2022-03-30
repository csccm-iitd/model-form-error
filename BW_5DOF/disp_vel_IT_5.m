function y = disp_vel_IT_5(y,a_g,r1,r2,m,c,k)
global dt; sigma = [0.01*ones(5,1); 0.001];

y1 = y(1); y2 = y(2); y3 = y(3); y4 = y(4); y5 = y(5); y6 = y(6);
y7 = y(7); y8 = y(8); y9 = y(9); y10 = y(10); y11 = y(11);

m1 = m(1); m2 = m(2); m3 = m(3); m4 = m(4); m5 = m(5);
c1 = c(1); c2 = c(2); c3 = c(3); c4 = c(4); c5 = c(5);
k1 = k(1); k2 = k(2); k3 = k(3); k4 = k(4); k5 = k(5);

sigma1 = sigma(1); sigma2 = sigma(2); sigma3 = sigma(3);
sigma4 = sigma(4); sigma5 = sigma(5); sigma6 = sigma(6);

g = 9.81; Qy = 0.05*(m1+m2+m3+m4+m5)*g; alpha = 1;
beta = 0.5; gamma = 0.5; eta = 1; Dy = 0.013; kr = 1/6; 

y = [y1;y2;y3;y4;y5;y6;y7;y8;y9;y10;y11];

a = [y(6:10);
    (1/m1)*(-m1*g*a_g-(c1*y6+k1*y1+c2*(y6-y7)+k2*(y1-y2)+(1-kr)*Qy*y11));
    (1/m2)*(-m2*g*a_g-(c2*(y7-y6)+k2*(y2-y1)+c3*(y7-y8)+k3*(y2-y3)));
    (1/m3)*(-m3*g*a_g-(c3*(y8-y7)+k3*(y3-y2)+c4*(y8-y9)+k4*(y3-y4)));
    (1/m4)*(-m4*g*a_g-(c4*(y9-y8)+k4*(y4-y3)+c5*(y9-y10)+k5*(y4-y5)));
    (1/m5)*(-m5*g*a_g-(c5*(y10-y9)+k5*(y5-y4)));
    (1/Dy)*(alpha*y6-gamma*y11*abs(y6)*abs(y11)^(eta-1)-beta*y6*abs(y11)^eta)];

b = [zeros(5,1); sigma(1:5)./m; sigma6];

L0a = [-(c1*y6 + k1*y1 + c2*(y6 - y7) + k2*(y1 - y2) + a_g*g*m1 - Qy*y11*(kr - 1))/m1;
    -(c3*(y7 - y8) - c2*(y6 - y7) - k2*(y1 - y2) + k3*(y2 - y3) + a_g*g*m2)/m2;
    -(c4*(y8 - y9) - c3*(y7 - y8) - k3*(y2 - y3) + k4*(y3 - y4) + a_g*g*m3)/m3;
    -(c5*(y9 - y10) - c4*(y8 - y9) - k4*(y3 - y4) + k5*(y4 - y5) + a_g*g*m4)/m4;
    (c5*(y9 - y10) + k5*(y4 - y5) - a_g*g*m5)/m5;
    (k2*y7)/m1 + ((c1 + c2)*(c1*y6 + k1*y1 + c2*(y6 - y7) + k2*(y1 - y2) + a_g*g*m1 - Qy*y11*(kr - 1)))/m1^2 - (y6*(k1 + k2))/m1 - (c2*(c3*(y7 - y8) - c2*(y6 - y7) - k2*(y1 - y2) + k3*(y2 - y3) + a_g*g*m2))/(m1*m2) - (Qy*(kr - 1)*(beta*y6*abs(y11)^eta - alpha*y6 + gamma*y11*abs(y6)*abs(y11)^(eta - 1)))/(Dy*m1);
    (k2*y6)/m2 + (k3*y8)/m2 + ((c2 + c3)*(c3*(y7 - y8) - c2*(y6 - y7) - k2*(y1 - y2) + k3*(y2 - y3) + a_g*g*m2))/m2^2 - (y7*(k2 + k3))/m2 - (c3*(c4*(y8 - y9) - c3*(y7 - y8) - k3*(y2 - y3) + k4*(y3 - y4) + a_g*g*m3))/(m2*m3) - (c2*(c1*y6 + k1*y1 + c2*(y6 - y7) + k2*(y1 - y2) + a_g*g*m1 - Qy*y11*(kr - 1)))/(m1*m2);
    (k3*y7)/m3 + (k4*y9)/m3 + ((c3 + c4)*(c4*(y8 - y9) - c3*(y7 - y8) - k3*(y2 - y3) + k4*(y3 - y4) + a_g*g*m3))/m3^2 - (y8*(k3 + k4))/m3 - (c3*(c3*(y7 - y8) - c2*(y6 - y7) - k2*(y1 - y2) + k3*(y2 - y3) + a_g*g*m2))/(m2*m3) - (c4*(c5*(y9 - y10) - c4*(y8 - y9) - k4*(y3 - y4) + k5*(y4 - y5) + a_g*g*m4))/(m3*m4);
    (k4*y8)/m4 + (k5*y10)/m4 + ((c4 + c5)*(c5*(y9 - y10) - c4*(y8 - y9) - k4*(y3 - y4) + k5*(y4 - y5) + a_g*g*m4))/m4^2 - (y9*(k4 + k5))/m4 - (c4*(c4*(y8 - y9) - c3*(y7 - y8) - k3*(y2 - y3) + k4*(y3 - y4) + a_g*g*m3))/(m3*m4) + (c5*(c5*(y9 - y10) + k5*(y4 - y5) - a_g*g*m5))/(m4*m5);
    (k5*y9)/m5 - (k5*y10)/m5 - (c5*(c5*(y9 - y10) + k5*(y4 - y5) - a_g*g*m5))/m5^2 - (c5*(c5*(y9 - y10) - c4*(y8 - y9) - k4*(y3 - y4) + k5*(y4 - y5) + a_g*g*m4))/(m4*m5);
    ((gamma*abs(y6)*abs(y11)^(eta - 1) + beta*eta*y6*abs(y11)^(eta - 1)*sign(y11) + gamma*y11*abs(y6)*abs(y11)^(eta - 2)*sign(y11)*(eta - 1))*(beta*y6*abs(y11)^eta - alpha*y6 + gamma*y11*abs(y6)*abs(y11)^(eta - 1)))/Dy^2 - (sigma6^2*(2*gamma*abs(y6)*abs(y11)^(eta - 2)*sign(y11)*(eta - 1) + 2*beta*eta*y6*abs(y11)^(eta - 1)*dirac_fn(y11) + beta*eta*y6*abs(y11)^(eta - 2)*sign(y11)^2*(eta - 1) + 2*gamma*y11*abs(y6)*abs(y11)^(eta - 2)*dirac_fn(y11)*(eta - 1) + gamma*y11*abs(y6)*abs(y11)^(eta - 3)*sign(y11)^2*(eta - 1)*(eta - 2)))/(2*Dy) + ((beta*abs(y11)^eta - alpha + gamma*y11*abs(y11)^(eta - 1)*sign(y6))*(c1*y6 + k1*y1 + c2*(y6 - y7) + k2*(y1 - y2) + a_g*g*m1 - Qy*y11*(kr - 1)))/(Dy*m1) - (sigma1*sigma6*(gamma*abs(y11)^(eta - 1)*sign(y6) + beta*eta*abs(y11)^(eta - 1)*sign(y11) + gamma*y11*abs(y11)^(eta - 2)*sign(y6)*sign(y11)*(eta - 1)))/(Dy*m1) - (gamma*sigma1^2*y11*abs(y11)^(eta - 1)*dirac_fn(y6))/(Dy*m1^2)];

L1a = [sigma1/m1;
    sigma2/m2;
    sigma3/m3;
    sigma4/m4;
    sigma5/m5;
    (Qy*sigma6*(kr - 1))/m1 - (sigma1*(c1 + c2))/m1^2 + (c2*sigma2)/(m1*m2);
    (c2*sigma1)/(m1*m2) - (sigma2*(c2 + c3))/m2^2 + (c3*sigma3)/(m2*m3);
    (c3*sigma2)/(m2*m3) - (sigma3*(c3 + c4))/m3^2 + (c4*sigma4)/(m3*m4);
    (c4*sigma3)/(m3*m4) - (sigma4*(c4 + c5))/m4^2 + (c5*sigma5)/(m4*m5);
    (c5*sigma4)/(m4*m5) - (c5*sigma5)/m5^2;
    - (sigma6*(gamma*abs(y6)*abs(y11)^(eta - 1) + beta*eta*y6*abs(y11)^(eta - 1)*sign(y11) + gamma*y11*abs(y6)*abs(y11)^(eta - 2)*sign(y11)*(eta - 1)))/Dy - (sigma1*(beta*abs(y11)^eta - alpha + gamma*y11*abs(y11)^(eta - 1)*sign(y6)))/(Dy*m1)];

dz1 = (dt^(3/2)*r1)/2 + (3^(1/2)*dt^(3/2)*r2)/6; dw = dt^(1/2)*r1;

y = y+a*dt+b*dw+L1a*dz1+0.5*L0a*dt^2;
end
