%% constants
clear all;
close all;

global mc_steps Nx Ny Nz  J kT m prob 

% Initialization
J = 0.5;  %why
kT = 0.2; %why
mu_0 = 4*pi*1e-7; 

%core properties

NofT = 10;  %Number of turns
Lc = 0.01;  %length of core
A = 1e-6;   %Area of the core
R = 1e-2;   %Resistance in parallel

%scaling
Hc = 1e3; 
Mc = 1e5;

%time definition and switching speed need to be defined
T = 0.5e-4;
dt = 0.2e-6;
n = floor(T/dt);

sw_speed = 1e7; %number of switches in one second
mc_steps = sw_speed*dt ;

% Create a folder to save the frames
frameFolder = 'spin_frames';
if ~exist(frameFolder, 'dir')
    mkdir(frameFolder); %make directory of the file
end
%% Start of process

Io = 3;
I = [Io*ones(1,n),-Io*ones(1,n),Io*ones(1,n)];
k = 3*n; %since there are 3 cycle of step
t = linspace(0,dt*k,k); %time linspace

Nx = 25;
Ny = 25;
Nz = 70;
probSpinUp = 0.5;
spin(1,:,:,:) = sign(probSpinUp - rand(Nx,Ny,Nz));    %create a spin array

%%
I_L = [0,0];
dphi = 0;

%pre allocating to increase speed
H = zeros(1,k);
M = zeros(1,k);
B = zeros(1,k);
phi = zeros(1,k);
IR = zeros(1,k);
V = zeros(1,k);
IL = zeros(1,k);
L = zeros(1,k);

for i = 2:10
    clear I_L;
    I_L(1) = I(i)*0.1;
    diff = 1;
    spin_new = spin(i-1,:,:,:);
    while(diff>10^-9)
        spin_prev = spin_new;
        H(i) = (I_L(1)*NofT/Lc)/Hc;
        spin_new = monte_carlo(H(i),spin_prev);  %compute the new spin using given conditions
        M(i) = (sum(sum(sum(spin_new(:,:,:))))/(Nx*Ny*Nz))*Mc;
        B(i) = mu_0*(M(i) + H(i));
        phi(i) = B(i)*A;
        IR(i) = ((phi(i)-phi(i-1))/dt)/R;
        I_L(2) = I(i) - IR(i);
        diff = abs(I_L(2)-I_L(1));
        I_L(1) = I_L(2);
    end
    spin(i,:,:,:) = spin_new;
    V(i) = IR(i)*R;
    IL(i) = I_L(2);
    L(i) = dphi/(IL(i)-IL(i-1));
    disp(i);
    display(H(i));
    display(M(i));
    display(B(i));
    display(phi(i));
    display(IR(i));
    display(V(i));
end
lw = 3;

figure(2)
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
plot(t,IL,'- r','LineWidth',lw);
xlabel('time (s)')
ylabel('IL')
xlim([0, t(end)])

subplot(5,1,5)
plot(t,phi,'- r','LineWidth',lw);
xlabel('time (s)')
ylabel('flux phi')
xlim([0, t(end)])

figure(3)
plot(IL,phi,'LineWidth',lw)

%% Function
function spin = monte_carlo(H,spin)
global mc_steps Nx J kT k Ny Nz m prob  
    % Metropolis algorithm  
    for iter = 1:mc_steps
    
        linearIndex = randi(numel(spin));   %select a spin randomly
        [row, col, z] = ind2sub(size(spin), linearIndex);   %convert the number to equalent position
    
        % Find its nearest neighbors
        above = mod(row - 1 - 1, size(spin, 1)) + 1;
        below = mod(row + 1 - 1, size(spin, 1)) + 1;
        left = mod(col - 1 - 1, size(spin, 2)) + 1;
        right = mod(col + 1 - 1, size(spin, 2)) + 1;
        front = mod(z - 1 - 1, size(spin, 3)) + 1;
        back = mod(z + 1 - 1, size(spin, 3)) + 1;
    
        neighbors = [spin(above, col, z); spin(below, col, z); spin(row, left, z);
                     spin(row, right, z); spin(row, col, front);spin(row, col, back)];
    
        % Calculate energy change if this spin is flipped
        dE = 2 * J * spin(row, col, z) .* sum(neighbors) + 2*H*spin(row, col, z);
    
        % Boltzmann probability of flipping
        prob = exp(-dE / kT);
    
        % Spin flip condition
        if dE <= 0 || rand() <= prob
            spin(row, col, z) = -spin(row, col, z);
        end
    end
end

%     captureInterval = numel(spin);     % Capture a frame every numel(spin) iterations
%     % Update the plot every captureInterval iterations
%     if mod(iter, captureInterval) == 0
%         set(h, 'CData', spin(:));  % Update the color data
% 
%         % Save the frame as a PNG image
%         filename = fullfile(frameFolder, sprintf('frame_%04d.png', iter));
%         saveas(gcf, filename);
%     end
