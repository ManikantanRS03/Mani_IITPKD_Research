% Parameters
Ts = 1/50;           % Sampling time
Kp = 0.2;            % Proportional gain (you can tune this)
Ki = 0.08;            % Integral gain (you can tune this)
N = 100;              % Moving average window length

% Discrete plant (gain only)
G = tf([0 0.7], 1, Ts);  % Plant is just a constant gain

% PI Controller with Tustin (bilinear) approximation for integral
s = tf('s');
PI_cont = Kp + Ki/s;                   % Continuous-time PI
PI_disc = c2d(PI_cont, Ts, 'tustin');  % Discretize using Tustin method

% Moving Average Filter in Feedback Path
b = ones(1, N)/N;         % Numerator of MA filter
a = 1;                    % Denominator (FIR filter)
H = tf(b, a, Ts);         % Feedback filter transfer function

P_TRNG = G*H;

% Create open-loop system with feedback filter
% Structure: u --> [PI] --> [G] --> y --> [H] --> back to error
Loop = series(PI_disc, G);
OpenLoop_with_H = feedback(Loop, H);  % Close loop with feedback filter
% sisotool(Loop)  % You tune PI here, H is fixed in external feedback

%%
sisotool(P_TRNG,PI_disc,1) %using sisotool to set the control parameters.