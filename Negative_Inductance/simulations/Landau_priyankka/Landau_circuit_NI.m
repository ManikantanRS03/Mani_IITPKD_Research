clear all
close all

alpha = -0.15e4;
beta = 0.75e12;
R=1e4;
RL=0;%1e6;%100e3;
rho=1e-6;%1/RL*0.99;

n=0.5e2;
Ival=3; 

I = [Ival.*ones(1,n),-Ival.*ones(1,n),Ival.*ones(1,n)];
dt=1e-8;
t = linspace(0,dt*length(I),length(I));

%  f=1e7;
% I = Ival*sin(2*pi*f*t);
fc=1e2;
%t=linspace(0,2*(1/fc),fc*2.5);
%I=Ival*sawtooth(2*pi*fc*t,0.5);
dt=1e-9;%1/fc/2.5;
threshold=1e-7;
IL(1)=0;

phi(1)=0;%-5.2e-5;
for i = 2:length(I)
    
    j=1;
    ILguess(1) = 0.5*I(i);%IL(i-1);
    diff =1;
    while(diff>threshold)
        VL(i) = R*(I(i)-ILguess(j)) - ILguess(j)*(RL);
        disp("VL =  ");
        disp(VL(i));
        V(i)= R*(I(i)-ILguess(j));
        phi(i)=phi(i-1) + VL(i)*dt;
        disp("phi =  ");
        disp(phi(i));
        j=j+1;
        ILguess(j) = (phi(i) - phi(i-1))*rho/dt + (2*alpha*phi(i-1) + 4*beta*phi(i-1)^3) ;
        diff=abs(ILguess(j)-ILguess(j-1));
        disp(diff);
    end

%     figure(3)
%     plot(ILguess)
%     hold on
    IL(i)=ILguess(end);
    disp(IL(i));
    L(i)=(phi(i)-phi(i-1))/(IL(i)-IL(i-1));
    U(i)=alpha*phi(i)^2 + beta*phi(i)^4 - IL(i)*phi(i);
    o=i


end

lw=3;

figure(1)
subplot(6,1,1)
plot(t,I,'-','LineWidth',lw);
xlabel('time (s)')
ylabel('I')
xlim([0, t(end)])

subplot(6,1,2)
plot(t,V,'-','LineWidth',lw);
xlabel('time (s)')
ylabel('V')
xlim([0, t(end)])

subplot(6,1,3)
plot(t,VL,'-','LineWidth',lw);
xlabel('time (s)')
ylabel('V_L')
xlim([0, t(end)])

subplot(6,1,4)
plot(t,IL,'- r','LineWidth',lw);
xlabel('time (s)')
ylabel('I_L')
xlim([0, t(end)])

subplot(6,1,5)
plot(t,phi,'- r','LineWidth',lw);
xlabel('time (s)')
ylabel('flux phi')
xlim([0, t(end)])

subplot(6,1,6)
plot(t,L,'-','LineWidth',lw);
xlabel('time (s)')
ylabel('L')
xlim([0, t(end)])

figure(2)
plot(IL,phi)

figure(3)
subplot(3,1,1)
plot(t,I,'-','LineWidth',lw);
xlabel('time (s)')
ylabel('I')
xlim([0, t(end)])

subplot(3,1,2)
plot(t,IL,'- r','LineWidth',lw);
xlabel('time (s)')
ylabel('I_L')
xlim([0, t(end)])

subplot(3,1,3)
plot(t,phi,'- r','LineWidth',lw);
xlabel('time (s)')
ylabel('flux phi')
xlim([0, t(end)])



%% double well

lw=2;
phi_dw = linspace(-10e-5,10e-5,100);
IL_dw=[0.27, -0.27];

for j = 1: length(IL_dw)
    for i = 1:length(phi_dw)
        Energy(j,i) = alpha*phi_dw(i)^2 + beta*phi_dw(i)^4 - IL_dw(j)*phi_dw(i);
    end
    %[Emin,i_min(j)]= min(Energy(j,:));
end



%L = (phi_dw(2)-phi_dw(1))/(IL_dw(2)-IL_dw(1))

% figure(2)
% plot(phi_dw,Energy(2,:),'LineWidth',lw)
% hold on
% plot(phi_dw,Energy(1,:),'LineWidth',lw)
% plot(phi,U,'LineWidth',lw)


