%%
clc;
clear all;
close all;
%%

c1=1.65;
c2=2.5;
c3=1.1;
sigma1=0.2;
sigma2=0.1;
sigma3=0.1;
gamma1=2;
gamma2=0.0001;
gamma3=0.000000001;
varepsilon2=10;
varepsilon3=0.0002;

keta=0.9;
J=1.625*10^(-3);
m=0.506;
M0=0.1;
L0=0.1;
R0=0.023;
b=16.25*10^(-3);
L=2*10^(-3);
R=5;
KB=keta;
B=b/keta;
G=9.8;
M=(15*J+5*m*L0^2+15*M0*(L0^2)+6*M0*(R0^2))/(15*keta);
N=m*L0*G/(2*keta)+M0*L0*G/keta; 
B=b/keta;
g2=1/M;
g3=1/L;
theta1=-N/M;
theta2=-B/M;
theta3=-KB/L;
theta4=-R/L;
%% initial values
x1=0.3;
x2=-0.2;
x3=-0.6;
hat_b1=0;
hat_b2=0;
hat_b3=0;
%%
step=0.0001;
n=100000;
t=0:step:step*n;

for i=1:1:n;
      tt=i*step;
      %% desired signal
      yd=0.1+0.5*sin(tt);
      dyd=0.5*cos(tt);
      yds(i)=yd;
      x1s(i)=x1;
      x2s(i)=x2;
      x3s(i)=x3;
% 
%   %% time varying constraint functions
    %% F11=0.4-0.7*sin(tt);
    F11=-0.4+0.2*sin(tt);
  %%  F12=0.8+0.2*sin(tt);
    F12=-0.25*exp(-tt)-0.25*sin(tt)+0.25*cos(tt) + 1.2;
    F21=-0.5*sin(tt)+0.3*cos(tt)-1.8;
    F22=-0.6*exp(-tt/2) - 0.2*sin(tt) + 0.4*cos(tt) + 1;
    F31=-0.5*sin(tt)+0.3*cos(tt)-2;
    F32=-0.6*exp(-tt/2) - 0.2*sin(tt) + 0.4*cos(tt) + 1;
%     %% the boundary of time-varying constraints
    uF11=-0.3;
    uF12=0.6;
    uF21=-0.3;
    uF22=0.6;    
    uF31=-0.3;
    uF32=0.6; 
    %% the derivative of time-varying constraint functions
    dF11=0.5*cos(tt);
    dF12=0.65*exp(-tt)-0.25*cos(tt)-0.25*sin(tt);

    dF21=-cos(tt)+sin(tt);
    dF22=0.3*exp(-tt/2)-0.2*cos(tt)-0.4*sin(tt);

    dF31=-cos(tt)+sin(tt);
    dF32=0.3*exp(-tt/2)-0.2*cos(tt)-0.4*sin(tt);
           
      %% constraint boundaries
      F11s(i)=F11;
      F12s(i)=F12;  
      F21s(i)=F21;
      F22s(i)=F22; 
      F31s(i)=F31;
      F32s(i)=F32; 


      %% transformation
      rho1=(x1-uF11)/(x1-F11)+(x1-uF12)/(F12-x1);
      rho1s(i)=rho1;
      rho2=(x2-uF21)/(x2-F21)+(x2-uF22)/(F22-x2);
      rho2s(i)=rho2;
      rho3=(x3-uF31)/(x3-F31)+(x3-uF32)/(F32-x3);
      rho3s(i)=rho3;

      alpha1f=(-uF11+yd)/(-F11+yd)+(yd-uF12)/(F12-yd);
      eta01=(F11-uF11)/((-F11+yd)^2)+(F12-uF12)/((F12-yd)^2);
      eta02=dF11*(-uF11+yd)/((-F11+yd)^2)-dF12*(yd-uF12)/((F12-yd)^2);
      dalpha1f=eta01*dyd+eta02;
%      alpha1f=yd/((F11+yd)*(F12-yd));
%      dalpha1f=(dyd*(F11*F12+yd^2))/((F11+yd)^2*(F12-yd)^2);

      rho11=(-F11+uF11+F12-uF12)/((x1-F11)*(F12-x1));
      rho21=(-F21+uF21+F22-uF22)/((x2-F21)*(F22-x2));
      rho31=(-F31+uF31+F32-uF32)/((x3-F31)*(F32-x3));

      rho12=(-uF11*F12+F11*uF12)/((x1-F11)*(F12-x1));
      rho22=(-uF21*F22+F21*uF22)/((x2-F21)*(F22-x2));
      rho32=(-uF11*F32+F31*uF32)/((x3-F31)*(F32-x3));

      rho13=(-uF11*F12+F11*uF12)/(F11-uF11+F12-uF12);
      rho23=(-uF21*F22+F21*uF22)/(F21-uF21+F22-uF22);
      rho33=(-uF11*F32+F31*uF32)/(F31-uF31+F32-uF32);

      miu11=(x1^2-F11*F12)/(((x1-F11)^2)*((F12-x1)^2));
      miu21=(x2^2-F21*F22)/(((x2-F21)^2)*((F22-x2)^2));
      miu31=(x3^2-F31*F32)/(((x3-F31)^2)*((F32-x3)^2));

      miu12=((F12*dF11+dF12*F11)*x1-((dF12+dF11)*(x1^2)))/(((x1-F11)^2)*((F12-x1)^2));
      miu22=((F22*dF21+dF22*F21)*x2-((dF22+dF21)*(x2^2)))/(((x2-F21)^2)*((F22-x2)^2));
      miu32=((F32*dF31+dF32*F31)*x3-((dF32+dF31)*(x3^2)))/(((x3-F31)^2)*((F32-x3)^2));


%       %% tracking error       
      z1=rho1-alpha1f;
      z1s(i)=z1;


%       %% true tracking error      
      e=x1-yd;
      es(i)=e;
      phi1=dalpha1f^2+miu12^2+(miu11^2)*((rho22/rho21)^2)+miu11^2;

      
      alpha1=(-c1*z1-z1*phi1)/miu11;
      alpha1s(i)=alpha1;

     

      if i<2
      alpha2f=alpha1;
      end

      alpha2fs(i)=alpha2f;
      z2=rho2-alpha2f;
      z2s(i)=z2;


      dot_alpha2f=(rho21*alpha1-alpha2f)/varepsilon2;
      phi2=miu21^2+(miu21^2)*(rho32/rho31)^2+((1+x2^2)^2)*(miu21^2)+(miu22-dot_alpha2f)^2+((miu11^2)*(z1^2))/(rho21^2);
      alpha2=(-c2*z2-hat_b2*z2*phi2)/miu21;
      alpha2f=alpha2f+dot_alpha2f*step; 

      if i<2
      alpha3f=alpha2;
      end
      
      alpha3fs(i)=alpha3f;
      z3=rho3-alpha3f;
      z3s(i)=z3;

      dot_alpha3f=(rho31*alpha2-alpha3f)/varepsilon3;
      phi3=(miu31^2)*((x2^2)+(x3^2)+1)+(miu32-dot_alpha3f)^2+((miu21^2)*(z2^2)*(z3^2))/(rho31^2);
      alpha3f=alpha3f+dot_alpha3f*step; 

      hat_b2s(i)=hat_b2;
      dot_hatb2=gamma2*(z2^2)*phi2-sigma2*hat_b2;

      hat_b3s(i)=hat_b3;
      dot_hatb3=gamma3*(z3^2)*phi3-sigma3*hat_b3;

     
%     hat_b1=hat_b1+dot_hatb1*step;
      hat_b2=hat_b2+dot_hatb2*step;
      hat_b3=hat_b3+dot_hatb3*step;
      

%% plant
u=(-c3*z3-hat_b3*z3*phi3)/miu31;
us(i)=u;

x1=x1+x2*step;
x2=x2+(g2*x3+theta1*sin(x1)+theta2*x2)*step;
x3=x3+(g3*u+theta3*x2+theta4*x3)*step;
    
end

figure(2);
h=gca;
set(h,'FontSize',14);
plot(t(1:n),x1s,'r',t(1:n),yds,'g--',t(1:n),F11s,'b--',t(1:n),F12s,'m--','linewidth',2);
legend({'$x_1$','$y_d$','$\kappa{_{1l}}$','$\kappa{_{1h}}$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
grid on;

figure(3);
h=gca;
set(h,'FontSize',14);
plot(t(1:n),x2s,'r',t(1:n),F21s,'b--',t(1:n),F22s,'m--','linewidth',2);
legend({'$x_2$','$\kappa{_{2l}}$','$\kappa{_{2h}}$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
grid on;

figure(4);
h=gca;
set(h,'FontSize',14);
plot(t(1:n),x3s,'r',t(1:n),F31s,'b--',t(1:n),F32s,'m--','linewidth',2);
legend({'$x_3$','$\kappa{_{3l}}$','$\kappa{_{3h}}$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
grid on;

figure(5);
subplot(311);
h=gca;
set(h,'FontSize',14);
plot(t(1:n),es,'r--','linewidth',2);
legend({'$z_1$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
grid on;
subplot(312);
h=gca;
set(h,'FontSize',14);
plot(t(1:n),z2s,'r','linewidth',2);
legend({'$z_2$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',10);
grid on;
subplot(313);
h=gca;
set(h,'FontSize',14);
plot(t(1:n),z3s,'r','linewidth',2);
legend({'$z_3$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',10);
grid on;


figure(6);
subplot(211);
h=gca;
set(h,'FontSize',14);
plot(t(1:n),es,'r','linewidth',2);
legend({'$e=x_1-y_d$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
%ylabel('e=x_1-y_d');
grid on;
subplot(212);
h=gca;
set(h,'FontSize',14);
plot(t(1:n),us,'r','linewidth',2);
legend({'$u$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
%ylabel('u','FontName','Times New Roman','FontSize',14);
grid on;

figure(7);
subplot(311);
h=gca;
set(h,'FontSize',14);
plot(t(1:n),alpha1s,'r',t(1:n),F11s,'b--',t(1:n),F12s,'m--','linewidth',2);
legend({'$\alpha_1$','$\kappa{_{2l}}$','$\kappa{_{2h}}$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
grid on;
subplot(312);
h=gca;
set(h,'FontSize',14);
plot(t(1:n),alpha2fs,'r',t(1:n),F11s,'b--',t(1:n),F12s,'m--','linewidth',2);
legend({'$\alpha_2$','$\kappa{_{2l}}$','$\kappa{_{2h}}$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
grid on;
subplot(313);
h=gca;
set(h,'FontSize',14);
plot(t(1:n),alpha3fs,'r',t(1:n),F11s,'b--',t(1:n),F12s,'m--','linewidth',2);
legend({'$\alpha_3$','$\kappa{_{2l}}$','$\kappa{_{2h}}$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
grid on;

figure(8);
h=gca;
set(h,'FontSize',14);
plot3(F11s,F12s,rho1s,'m--','linewidth',2);
% legend({'${\alpha}_{2f}$'},'Interpreter','latex');
% xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
grid on;