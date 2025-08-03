
global mc_steps Nx Ny Nz  J T m prob 
%H0 = 1; %extreme magnetic field
T=0.2; %tempreture

k=1;
J=0.5;
mu_0=4*pi*1e-7;


%core properties
NuT = 10;
Lc= 0.01;


R=1e-2;


m=0;


%% Start of process



for i =2:length(I)
    clear IL
    IL(1)=I(i)*0.1;%first guess
    j=1;
    lattice_next=lattice(i-1,:,:,:);
    Diff=1;

    while(Diff>1e-9)

        lattice_prev=lattice_next;
        
        %________Input IL to this section and get phi_____
        H(i) = IL(j)*NuT/Lc;
        lattice_next = monte_carlo(H(i)/Hc,lattice_prev);
        M(i) = (sum(sum(sum(lattice_next(:,:,:,:))))/(Nx*Ny*Nz))*Mc;
        B(i) = mu_0*(H(i)+M(i));
        Phi(i)=B(i)*A;
        %_______end of section____________________________
        j=j+1;
        dphi(i)=Phi(i)-Phi(i-1);
        IR(i)= (Phi(i)-Phi(i-1))/dt/R;
        IL(j)= I(i) - IR(i);
        V(i)=IR(i)*R;
        dV(i)=V(i)-V(i-1);
    
        Diff = abs(IL(j)-IL(j-1));
        

    end
    I_L(i)=IL(end);
    lattice(i,:,:,:)=lattice_next;
    L(i)=dphi(i)/(I_L(i)-I_L(i-1));
    g=i
end

% X=linspace(1,Nx,Nx);
% Y=linspace(1,Ny,Ny);
% Z=linspace(1,Nz,Nz);
% 
% x=lattice(320,:,:,:);%40,270,300
% V=squeeze(x);
% fig = figure ; 
% 
% xslice = [1];   
% yslice = [1];
% zslice = [Nz];
% slice(X,Y,Z,V,xslice,yslice,zslice)
% axis vis3d
% camproj('perspective')

lw=3;

%%

alpha = -0.15e4;
beta = 0.75e12;

k = 1;
Time=0.5e-4;
dt=0.2e-6;
t = linspace(0,dt*k,k); %time linspace

for i = 2:length(I)
    U(i) = alpha*Phi(i)^2 + beta*Phi(i)^4 - Phi(i)*I_L(i); 
end
lw=3;

figure(1)
subplot(5,1,1)
plot(t,I,'-','LineWidth',lw);
xlabel('time (s)')
ylabel('I')
xlim([0, t(end)])

subplot(5,1,2)
plot(t,V,'-','LineWidth',lw);
xlabel('time (s)')
ylabel('V')
xlim([0, t(end)])

subplot(5,1,3)
plot(t,V,'-','LineWidth',lw);
xlabel('time (s)')
ylabel('V_L')
xlim([0, t(end)])

subplot(5,1,4)
plot(t,I_L,'- r','LineWidth',lw);
xlabel('time (s)')
ylabel('I_L')
xlim([0, t(end)])

subplot(5,1,5)
plot(t,Phi,'- r','LineWidth',lw);
xlabel('time (s)')
ylabel('flux phi')
xlim([0, t(end)])

figure(2)
plot(I_L,Phi,'LineWidth',lw)

figure(3)
subplot(3,1,1)
plot(t,I,'-','LineWidth',lw);
xlabel('time (s)')
ylabel('I')
xlim([0, t(end)])

subplot(3,1,2)
plot(t,I_L,'- r','LineWidth',lw);
xlabel('time (s)')
ylabel('I_L')
xlim([0, t(end)])

subplot(3,1,3)
plot(t,Phi,'- r','LineWidth',lw);
xlabel('time (s)')
ylabel('flux phi')
xlim([0, t(end)])

figure(4)
plot(Phi,U,'LineWidth',lw)
%% Functions

function lattice = monte_carlo(H,lattice)
global mc_steps Nx J T k Ny Nz m prob    


    for i= 1:mc_steps
        X = randi([1,Nx]);
        Y = randi([1,Ny]);
        Z = randi([1,Nz]);

        % X direction 
        if X~=Nx
            L_right = lattice(1,X+1, Y,Z);
        end
        if X~=1
            L_left = lattice(1,X-1, Y,Z);
        end

        % Y direction
        if Y~=Ny
            L_bottom = lattice(1,X, Y+1,Z);
        end
        if Y~=1
            L_top = lattice(1,X,Y-1,Z);
        end

        % Z direction
        if Z~=Nz
            L_front = lattice(1,X,Y,Z+1) ;
        end
        if Z~=1
            L_back = lattice(1,X,Y,Z-1) ;
        end
  
        if X == 1
            L_left = lattice(1,Nx,Y,Z);
        elseif X==Nx
            L_right = lattice(1,1,Y,Z);
        end

        if Y == 1
            L_top = 0;
        elseif Y==Ny
            L_bottom = 0;
        end

        if Z == 1
            L_back = 0;
        elseif Z==Nz
            L_front = 0;
        end

        a=L_top;
        b=L_bottom;
        c=L_left;
        d=L_right;

        m=m+1;

        dU=2*J*(L_top+L_bottom+L_left+L_right+L_front+L_back)*lattice(1,X, Y,Z) + 2*H*lattice(1,X,Y,Z);

        if(dU<0)
            lattice(1,X,Y,Z)=-lattice(1,X,Y,Z);
        elseif rand < exp(-dU/T)
            lattice(1,X,Y,Z)=-lattice(1,X,Y,Z);
        end
    end
end


function lattice = initialstate(Nx,Ny,Nz)
num = floor((Nx+Ny+Nz)/2);
lattice = -ones(Nx,Ny,Nz);
    for i = 1:num
        X=randi([2,Nx]);
        Y=randi([2,Ny]);
        Z=randi([2,Nz]);
        lattice(X,Y,Z) = -1;
    end
end


