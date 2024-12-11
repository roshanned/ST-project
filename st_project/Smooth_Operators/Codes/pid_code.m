% Constants
l1 = 0.2;
l2 = 0.1;
m1 = 10;
m2 = 5;
g = 9.81;

% Initial State Variables
e1 = 0;
e2 = 0;
positions = [0.1, 0.1]; % q1 and q2
velocities = [0, 0]; % q1_dot and q2_dot
initial_q = [e1, e2, positions, velocities];

% PID Controller gains
kp = 100;
ki = 150;
kd = 20;

% Final angles
qfinal = [0; 0];

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
[time, s] = ode45(@(t, q) pid(t, q, qfinal, kp,kd,ki), [0,30], initial_q, options);

figure;
% subplot(2, 1, 1);
plot(time, s(:, 3));
xlabel("Time");
ylabel("q_1");
hold on;
% subplot(2, 1, 2);
plot(time, s(:, 4), 'r');
xlabel("Time");
ylabel("q_2");

sgtitle('Joint Angles vs Time - PID Controller');

% PID Controller Function
function q_derivative = pid(~, q, qfinal, kp,kd,ki)
    % Constants
    l1 = 0.2;
    l2 = 0.1;
    m1 = 10;
    m2 = 5;
    g = 9.81;
    
    % Error calculation
    e1 = qfinal(1) - q(3);
    e2 = qfinal(2) - q(4);

    % Mass Matrix
    M11 = (m1 + m2) * l1^2 + m2 * l2 * (l2 + 2 * l1 * cos(q(4)));
    M22 = m2 * l2^2;
    M12 = m2 * l2 * (l2 + l1 * cos(q(4)));
    M21 = M12;
    M = [M11, M12; M21, M22];

    % Torque
    f1 = kp * e1 - kd * q(5) + ki * q(1);
    f2 = kp * e2 - kd * q(6) + ki * q(2);
    f=[f1;f2];

    % Coriolis Matrix
    c11 = -m2 * l1 * l2 * sin(q(4)) * q(5);
    c12 = -m2 * l1 * l2 * sin(q(4)) * (q(5) + q(6));
    c21 = 0;
    c22 = m2 * l1 * l2 * sin(q(4)) * q(6);
    C = [c11, c12; c21, c22];

    % Gravitational Matrix
    G1 = m1 * l1 * g * cos(q(3)) + m2 * g * (l2 * cos(q(3) + q(4)) + l1 * cos(q(3)));
    G2 = m2 * l2 * g * cos(q(3) + q(4));
    G = [G1; G2];

    % State Derivatives
    dq = [q(5); q(6)];
    M_inv = inv(M);
    W = (M_inv) * ((-C * dq - G)) + f;

    % Next state
    q_derivative = zeros(6, 1);
    q_derivative(1) = e1;
    q_derivative(2) = e2;
    q_derivative(3) = q(5);
    q_derivative(4) = q(6);
    q_derivative(5:6) = W;
end