
%units Nx,Ny,Nz,n, -> unitless
%units T -> K
%units k -> J/K (value = 1.38065e-23)
%unit J -> 
%units mu0 -> H/m (value = 4*pi e-7)
%unit Lc -> m
%units A -> m^2
%units R -> 
%units Time,dt -> s


clear all
close all

global mc_steps Nx Ny Nz  J T

Nx=25; 
Ny=25;
Nz=70;

T = 0.2; %tempreture why??
k = 1;
J = 0.5;
mu_0 = 4*pi*1e-7;


%core properties

NuT = 10;
Lc = 0.01;
A = 1e-6;
R = 5e-2; %1e-2 %R parallel

% %scaling
% Hc = 1e3; 
% Mc = 1e5;

%scaling
Hc = 1e3/10; %1e3
Mc = 1e5*2.5; %1e5

%time definition and switching speed need to be defined
Time = 0.5e-4; %0.5e-4
dt = 2e-7; %2e-7
dt1 = 2e-7; %2e-7
n = floor(Time/dt);

sw_speed = 1e7; %1e7
mc_steps = sw_speed*dt ; %number of switches in iteration
%convergence condition
err = 1e-9;  %original 1e-9
iter_limt = 1e5;

Io = 3;

%% Importing Guess value for the previous simulations

file = 'C:\MyDrive\semester 4\oelp\simulations\Ising_results\Ising_guessvalues.csv';
data = readmatrix(file);
ILgs = data(1:end,3);
Phigs = data(1:end,2);
%Latticegs = lattice(1,:,:,:);

% display(I_L);
% display(length(I_L))

%% Start of process


I = [Io*ones(1,n),-Io*ones(1,n),Io*ones(1,n)];
k=3*n;

% Ilin = linspace(0,1,n/5);
% I = [Io*Ilin, Io*ones(1,n),-Io*ones(1,n),Io*ones(1,n)];
% k=3*n+(n/5);

% I = [Io*ones(1,n),-Io*ones(1,n)];
% k=2*n;

% I = zeros(1,k); %for triangular wave
% for v = 1:n
%     I(v) = -Io + (v-1)*2*Io/(n-1);
% end
% for v = 1:n
%     I(v+n) = Io - (v-1)*2*Io/(n-1);
% end
% for v = 1:n
%     I(v+2*n) = -Io + (v-1)*2*Io/(n-1);
% end

t=linspace(0,dt*k,k);

lattice(1,:,:,:) = initialstate(Nx,Ny,Nz);
%lattice(1,:,:,:) = Latticegs;

IL = [0,0];

%pre allocating to increase speed
H = zeros(1,k);
M = zeros(1,k);
B = zeros(1,k);
Phi = zeros(1,k);
IR = zeros(1,k);
V = zeros(1,k);
I_L = zeros(1,k);
L = zeros(1,k);
dIL = zeros(1,k);
j = zeros(1,k);
spin_ratio = zeros(1,k);

j1 = 0;
for i = 2:k
    clear IL
    IL(1) = I(i)*0.1; %first guess
    %IL(1) = ILgs(i); %using values from the previous iterations 
    lattice_next = lattice(i-1,:,:,:);
    Diff = 1;
    j(i) = j1;
    j1 = 0;
    while((Diff>err)&&(j1<iter_limt)) %either stop at convergence or at iteration limit
        j1 = j1+1;
        lattice_prev = lattice_next;
        H(i) = IL(1)*NuT/Lc;
        lattice_next = monte_carlo(H(i)/Hc,lattice_prev);
        M(i) = (sum(sum(sum(lattice_next(:,:,:,:))))/(Nx*Ny*Nz))*Mc;
        B(i) = mu_0*(H(i)+M(i));
        Phi(i) = B(i)*A;
        dphi = Phi(i)-Phi(i-1);
        IR(i) = (Phi(i)-Phi(i-1))/dt1/R; %dt
        IL(2) = I(i) - IR(i);
        V(i) = IR(i)*R;
        Diff = abs(IL(2)-IL(1));
        IL(1) = IL(2);
        %display(IL(1));
    end
    I_L(i) = IL(2);
    lattice(i,:,:,:) = lattice_next;
    dIL(i) = I_L(i)-I_L(i-1);
    L(i) = dphi/(dIL(i));
    spin_ratio(i) = (sum(sum(sum(lattice_next(:,:,:,:))))/(Nx*Ny*Nz));
    disp(i);
%     display(H(i));
%     display(M(i));
%     display(B(i));
%     display(Phi(i));
%     display(IR(i));
%     display(V(i));
end

%% Saving data to CSV

% directory = 'C:\MyDrive\semester 4\oelp\simulations\Ising_results';
% filename = sprintf('Ising_guessvalues.csv');
% columnname = {'time', 'Phi', 'I_L', 'IR', 'I', 'M', 'B', 'V'};
% 
% matrix = [t ;Phi ;I_L ;IR; I ;M ;B ;V];
% matrix = transpose(matrix);
% dataTable = array2table(matrix,'VariableNames', columnname );
% % Save the table to a CSV file
% fullPath = fullfile(directory, filename);
% writetable(dataTable, fullPath, 'Delimiter', ',');

%% Plotting

% % Create a folder to save the frames
% frameFolder = 'spin_frames';
% if ~exist(frameFolder, 'dir')
%     mkdir(frameFolder); %make directory of the file
% end
% 
% X=linspace(1,Nx,Nx);
% Y=linspace(1,Ny,Ny);
% Z=linspace(1,Nz,Nz);
% for itr = 1:150
%     figure(1);
%     clf;
%     x=lattice(itr*5,:,:,:);
%     sqX=squeeze(x);
%     fig = figure ; 
%     figure(1)
%     xslice = 1;   
%     yslice = 1;
%     zslice = Nz;
%     slice(X,Y,Z,sqX,xslice,yslice,zslice)
%     axis vis3d
%     camproj('perspective')
%     hold on;
% end

% %Plotting the lattice
% figure(1);
% clf;
% [x, y, z] = meshgrid(1:Nx, 1:Ny, 1:Nz);
% h = scatter3(x(:), y(:), z(:), 100, lattice(120,:),'filled');  % Use different colors for different spins
% colormap([0.027 0.157 0.529 ; 0.941 0.949 0.184]);
% axis equal;
% grid on;
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('Spin Lattice');
 
% % Update the plot every captureInterval iterations
% for iter = 1:150
%     set(h, 'CData', lattice(iter*5,:));  % Update the color data
%     % Save the frame as a PNG image
%     filename = fullfile(frameFolder, sprintf('frame_%04d.png', iter));
%     saveas(gcf, filename);
% end

lw = 3;
% figure(2)
% subplot(5,1,1)
% plot(t,I,'-','LineWidth',lw);
% xlabel('time (s)')
% ylabel('I')
% xlim([0, t(end)])
% 
% subplot(5,1,2)
% plot(t,V,'-','LineWidth',lw);
% xlabel('time (s)')
% ylabel('V')
% xlim([0, t(end)])
% 
% subplot(5,1,3)
% plot(t,IR,'-','LineWidth',lw);
% xlabel('time (s)')
% ylabel('IR')
% xlim([0, t(end)])
% 
% subplot(5,1,4)
% plot(t,I_L,'- r','LineWidth',lw);
% xlabel('time (s)')
% ylabel('I_L')
% xlim([0, t(end)])
% 
% subplot(5,1,5)
% plot(t,Phi,'- r','LineWidth',lw);
% xlabel('time (s)')
% ylabel('flux phi')
% xlim([0, t(end)])
% 
figure(3)
plot(I_L,Phi,'o--')
title("hysteresis")
xlabel("IL")
ylabel("Phi")

lw = 3;
figure(4)
subplot(2,1,1)
plot(t,I,'-r','LineWidth',lw);
hold on;
% xlabel('Time (s)',FontSize=11,FontName='Arial')
% ylabel('I (A)',FontSize=11,FontName='Arial')
% xlim([0, t(end)])

%subplot(2,1,2)
plot(t,I_L,'-b','LineWidth',lw);
xlabel('Time (s)',FontSize=11,FontName='Arial')
ylabel('I & I_L (A)',FontSize=11,FontName='Arial')
xlim([0, t(end)])
hold off;

subplot(2,1,2)
plot(t,Phi,'-b','LineWidth',lw);
xlabel('Time (s)',FontSize=11,FontName='Arial')
ylabel('\phi (Wb)',FontSize=11,FontName='Arial')
xlim([0, t(end)])

% figure(5)
% plot(Phi,(dIL/dt)./V,'o--')
% 
% figure(6)
% subplot(2,1,1)
% plot(t,I,'-','LineWidth',lw);
% xlabel('time (s)')
% ylabel('I')
% xlim([0, t(end)])
% 
% subplot(2,1,1)
% plot(t,I_L,'- r','LineWidth',lw);
% xlabel('time (s)')
% ylabel('I_L')
% xlim([0, t(end)])
% 
% subplot(2,1,2)
% plot(t,Phi,'- r','LineWidth',lw);
% xlabel('time (s)')
% ylabel('\phi')
% xlim([0, t(end)])

%plots to evaluate convergence

figure(6)
plot(t(1:10),j(1:10));
title("for error < 10^-9")
xlabel("time steps")
ylabel("number of iterations for convergence (1,10)")

figure(7)
plot(t(10:k),j(10:k));
title("for error < 10^-9")
xlabel("time steps")
ylabel("number of iterations for convergence (10,k)")

figure(8)
plot(t,spin_ratio);
title("spin ratio")
xlabel("time")
ylabel("sum of spins/total spins")

%% Functions

function lattice = monte_carlo(H,lattice)
global mc_steps Nx J T k Ny Nz   

    for i= 1:mc_steps
        X = randi([1,Nx]);
        Y = randi([1,Ny]);
        Z = randi([1,Nz]);

        % X direction 
        if X ~= Nx
            L_right = lattice(1,X+1, Y,Z);
        end
        if X ~= 1
            L_left = lattice(1,X-1, Y,Z);
        end

        % Y direction
        if Y ~= Ny
            L_bottom = lattice(1,X, Y+1,Z);
        end
        if Y ~= 1
            L_top = lattice(1,X,Y-1,Z);
        end

        % Z direction
        if Z ~= Nz
            L_front = lattice(1,X,Y,Z+1) ;
        end
        if Z ~= 1
            L_back = lattice(1,X,Y,Z-1) ;
        end
  
        if X == 1
            L_left = lattice(1,Nx,Y,Z);
        elseif X == Nx
            L_right = lattice(1,1,Y,Z);
        end

        if Y == 1
            L_top = 0;
        elseif Y == Ny
            L_bottom = 0;
        end

        if Z == 1
            L_back = 0;
        elseif Z == Nz
            L_front = 0;
        end

        dU=2*J*(L_top+L_bottom+L_left+L_right+L_front+L_back)*lattice(1,X, Y,Z) + 2*H*lattice(1,X,Y,Z);
        if(dU<0)
            lattice(1,X,Y,Z) = -lattice(1,X,Y,Z);
        elseif rand < exp(-dU/T)
            lattice(1,X,Y,Z) = -lattice(1,X,Y,Z);
        end
    end
end


function lattice = initialstate(Nx,Ny,Nz)
num = floor((Nx+Ny+Nz)/2);
lattice = -ones(Nx,Ny,Nz);
    for i = 1:num
        X = randi([2,Nx]);
        Y = randi([2,Ny]);
        Z = randi([2,Nz]);
        lattice(X,Y,Z) = -1;
    end
end

