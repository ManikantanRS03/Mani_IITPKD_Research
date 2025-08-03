% Initialization
J = 1;
numSpinsPerDim = 2^5;
probSpinUp = 0.5;
spin = sign(probSpinUp - rand(numSpinsPerDim, numSpinsPerDim, numSpinsPerDim));
kT = 1;

captureInterval = numel(spin);     % Capture a frame every numel(spin) iterations

% Create a folder to save the frames
frameFolder = 'spin_frames';
if ~exist(frameFolder, 'dir')
    mkdir(frameFolder);
end

% Metropolis algorithm
numIters = 2^7 * numel(spin);
for iter = 1:numIters
    % Pick a random spin
    linearIndex = randi(numel(spin));
    [row, col, z] = ind2sub(size(spin), linearIndex);

    % Find its nearest neighbors
    above = mod(row - 1 - 1, size(spin, 1)) + 1;
    below = mod(row + 1 - 1, size(spin, 1)) + 1;
    left = mod(col - 1 - 1, size(spin, 2)) + 1;
    right = mod(col + 1 - 1, size(spin, 2)) + 1;
    front = mod(z - 1 - 1, size(spin, 3)) + 1;
    back = mod(z + 1 - 1, size(spin, 3)) + 1;

    neighbors = [spin(above, col, z);
                 spin(below, col, z);
                 spin(row, left, z);
                 spin(row, right, z);
                 spin(row, col, front);
                 spin(row, col, back)];

    % Calculate energy change if this spin is flipped
    dE = 2 * J * spin(row, col, z) * sum(neighbors);

    % Boltzmann probability of flipping
    prob = exp(-dE / kT);

    % Spin flip condition
    if dE <= 0 || rand() <= prob
        spin(row, col, z) = -spin(row, col, z);
    end

    % Capture a frame every captureInterval iterations
    if mod(iter, captureInterval) == 0
        display(iter / captureInterval);

        % Plot the 3D lattice
        figure('Renderer', 'zbuffer');
        clf;
        xslice = [];  % Indices of X-axis slices (empty to include all)
        yslice = [];  % Indices of Y-axis slices (empty to include all)
        zslice = 1:size(spin, 3);  % Indices of Z-axis slices (all)

        slice(spin, xslice, yslice, zslice);
        shading interp;
        colormap([1 1 1; 0 0 1]);

        % Save the frame as a PNG image
        filename = fullfile(frameFolder, sprintf('frame_%04d.png', iter));
        saveas(gcf, filename);

        % Close the figure to avoid overlapping plots
        close(gcf);
    end
end