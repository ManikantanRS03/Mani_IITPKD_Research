% Example data
x = linspace(-10, 10, 100);
y1 = x.^2;                % Curve 1
y2 = -x.^2 + 20;          % Curve 2

% Calculate the slopes of the curves
dy1_dx = diff(y1) ./ diff(x);
dy2_dx = diff(y2) ./ diff(x);

% Plot the curves
plot(x, y1, 'b', 'LineWidth', 2);
hold on;
plot(x, y2, 'r', 'LineWidth', 2);
xlabel('x');
ylabel('y');
title('Plot with Highlighted Region');
legend('Curve 1', 'Curve 2');

% Find the region where slope of Curve 1 is negative and slope of Curve 2 is positive
highlightIndices = find(dy1_dx < 0 & dy2_dx > 0);
highlightX = x(highlightIndices);
highlightY1 = y1(highlightIndices);
highlightY2 = y2(highlightIndices);

% Plot the highlighted region
plot(highlightX, highlightY1, 'g', 'LineWidth', 2);
plot(highlightX, highlightY2, 'g', 'LineWidth', 2);

% Add a legend for the highlighted region
legend('Curve 1', 'Curve 2', 'Highlighted Region');
