function y = disp_vel_filter(y,a_g,m,c,k)
global dt; sigma = 0.01*ones(5,1);

y1 = y(1); y2 = y(2); y3 = y(3); y4 = y(4); y5 = y(5); y6 = y(6);
y7 = y(7); y8 = y(8); y9 = y(9); y10 = y(10);

m1 = m(1); m2 = m(2); m3 = m(3); m4 = m(4); m5 = m(5);
c1 = c(1); c2 = c(2); c3 = c(3); c4 = c(4); c5 = c(5);
k1 = k(1); k2 = k(2); k3 = k(3); k4 = k(4); k5 = k(5);

sigma1 = sigma(1);

g = 9.81; co = 10;
y = [y1;y2;y3;y4;y5;y6;y7;y8;y9;y10];

a = [y(6:10);
    (1/m1)*(-m1*g*a_g-(c1*y6+k1*y1+c2*(y6-y7)+k2*(y1-y2)+co*sign(y6)));
    (1/m2)*(-m2*g*a_g-(c2*(y7-y6)+k2*(y2-y1)+c3*(y7-y8)+k3*(y2-y3)));
    (1/m3)*(-m3*g*a_g-(c3*(y8-y7)+k3*(y3-y2)+c4*(y8-y9)+k4*(y3-y4)));
    (1/m4)*(-m4*g*a_g-(c4*(y9-y8)+k4*(y4-y3)+c5*(y9-y10)+k5*(y4-y5)));
    (1/m5)*(-m5*g*a_g-(c5*(y10-y9)+k5*(y5-y4)))];

L0a = [-(c1*y6 + k1*y1 + c2*(y6 - y7) + k2*(y1 - y2) + co*sign(y6) + a_g*g*m1)/m1;
    -(c3*(y7 - y8) - c2*(y6 - y7) - k2*(y1 - y2) + k3*(y2 - y3) + a_g*g*m2)/m2;
    -(c4*(y8 - y9) - c3*(y7 - y8) - k3*(y2 - y3) + k4*(y3 - y4) + a_g*g*m3)/m3;
    -(c5*(y9 - y10) - c4*(y8 - y9) - k4*(y3 - y4) + k5*(y4 - y5) + a_g*g*m4)/m4;
    (c5*(y9 - y10) + k5*(y4 - y5) - a_g*g*m5)/m5;
    ((c1 + c2 + 2*co*dirac_fn(y6))*(c1*y6 + k1*y1 + c2*(y6 - y7) + k2*(y1 - y2) + co*sign(y6) + a_g*g*m1))/m1^2 + (k2*y7)/m1 - (y6*(k1 + k2))/m1 - (c2*(c3*(y7 - y8) - c2*(y6 - y7) - k2*(y1 - y2) + k3*(y2 - y3) + a_g*g*m2))/(m1*m2) - (co*sigma1^2*(-dirac_fn(y6)/y6))/m1^3;
    (k2*y6)/m2 + (k3*y8)/m2 + ((c2 + c3)*(c3*(y7 - y8) - c2*(y6 - y7) - k2*(y1 - y2) + k3*(y2 - y3) + a_g*g*m2))/m2^2 - (y7*(k2 + k3))/m2 - (c3*(c4*(y8 - y9) - c3*(y7 - y8) - k3*(y2 - y3) + k4*(y3 - y4) + a_g*g*m3))/(m2*m3) - (c2*(c1*y6 + k1*y1 + c2*(y6 - y7) + k2*(y1 - y2) + co*sign(y6) + a_g*g*m1))/(m1*m2);
    (k3*y7)/m3 + (k4*y9)/m3 + ((c3 + c4)*(c4*(y8 - y9) - c3*(y7 - y8) - k3*(y2 - y3) + k4*(y3 - y4) + a_g*g*m3))/m3^2 - (y8*(k3 + k4))/m3 - (c3*(c3*(y7 - y8) - c2*(y6 - y7) - k2*(y1 - y2) + k3*(y2 - y3) + a_g*g*m2))/(m2*m3) - (c4*(c5*(y9 - y10) - c4*(y8 - y9) - k4*(y3 - y4) + k5*(y4 - y5) + a_g*g*m4))/(m3*m4);
    (k4*y8)/m4 + (k5*y10)/m4 + ((c4 + c5)*(c5*(y9 - y10) - c4*(y8 - y9) - k4*(y3 - y4) + k5*(y4 - y5) + a_g*g*m4))/m4^2 - (y9*(k4 + k5))/m4 - (c4*(c4*(y8 - y9) - c3*(y7 - y8) - k3*(y2 - y3) + k4*(y3 - y4) + a_g*g*m3))/(m3*m4) + (c5*(c5*(y9 - y10) + k5*(y4 - y5) - a_g*g*m5))/(m4*m5);
    (k5*y9)/m5 - (k5*y10)/m5 - (c5*(c5*(y9 - y10) + k5*(y4 - y5) - a_g*g*m5))/m5^2 - (c5*(c5*(y9 - y10) - c4*(y8 - y9) - k4*(y3 - y4) + k5*(y4 - y5) + a_g*g*m4))/(m4*m5)];

y = y+a*dt+0.5*L0a*dt^2;
end
