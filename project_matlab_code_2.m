clear all;
close all;
clc;

syms S
C_m =10; C_d = 3; J =5 ; C_t=8; d = 2 ; tau = 0.5;
A = [0 1 0;
     -C_m/J -C_d/J C_t*d/J;
     0 0 -1/tau]
B = [0;0;1/tau]

C = [1 0 0];
D = [0];
I = eye(3);
SI = (S*I);
disp("Find the characteristic equation for the Accf")
actual_characteristic_eqn = det(SI-A)
% Extract coefficients of the polynomial
actual_coeffs_poly = coeffs(actual_characteristic_eqn, S);

% Convert symbolic fractions to decimal
actual_coeffs_decimal = double(actual_coeffs_poly)

p_Mat = [0 0 6.4;
         0 6.4 -16.64;
         2 -4 8];
pCCFinv_Mat = [3.2 2.6 1;
                2.6 1 0;
                1 0 0];
pCCF = inv(pCCFinv_Mat);
det(pCCF);

tccf = p_Mat * pCCFinv_Mat;

tccf_inv = inv(tccf);
Accf = tccf_inv * A  * tccf;
Bccf = tccf_inv * B;
Cccf = C * tccf;

pretty(det(SI-Accf))
[V,D] = eig (A);
eigVals = diag(D);
eqn = (S-eigVals(1)) * (S-eigVals(2)) * (S-eigVals(3));






overshoot = 1;
settlingtime = 0.5;
zeta_desired = sqrt((log(overshoot/100))^2/ (pi^2 + (log(overshoot/100))^2))
omega_desired = 4 /(zeta_desired* settlingtime)


% Open−loop Response 
s = tf('s');

%%% Second order approximation:
disp("The third order transfer function is: ")
H = (C_t*d)/((tau*s+1)*(J*s^2+C_m+C_d*s))
% Define the coefficients of the denominator polynomial
denom = [J, C_d, C_m, tau];  % Coefficients of J*s^2 + C_d*s + C_m + tau*s + 1

% Find the roots of the denominator (which are the poles)
disp("The roots of the characteristic equation is:")
poles = roots(denom)
% poles(H)We c




%2ndorder
char_eqn_desired = S^2 + (2 * zeta_desired * omega_desired * S ) + (omega_desired^2)
roots_char_eqn_desired = roots(char_eqn_desired)

%3rd order 
full_char_eqn_desired = char_eqn_desired * (S + 80);

% Expand the equation to get a polynomial
expanded_eqn = expand(full_char_eqn_desired)


% Extract coefficients of the polynomial
coeffs_poly = coeffs(expanded_eqn, S);

% Convert symbolic fractions to decimal
alpha = double(coeffs_poly)
% coeffs_decimal(1)

% Find the Kccf matrix
Kccf = [alpha(1)-actual_coeffs_decimal(1) alpha(2)-actual_coeffs_decimal(2) alpha(3)-actual_coeffs_decimal(3)]

K = Kccf*inv(p_Mat*pCCFinv_Mat)
char_eqn__solution = solve(char_eqn_desired,S);


%H_OL = 1/(m*s^2+c*s+k); 
H_OL = (C_t*d)/(2.5*s^3 + 6.5 *s^2 + 8*s+ C_m); 
figure(1)
step(H_OL,20)
[y_OL,t1]=step(H_OL,20);


% Calculating the G value
DC_gain = 1;
G = (-1*DC_gain)/(C*inv(A-(B*K))*B)


% Closed−loop  
s = tf('s');
% H_CL = G/(2.5*s^3 + (6.5+K(2)) *s^2 + (8+(1))*s+ (10+K(3))); 
% H_CL = (C-D*K_1)*inv((s.*eye(3))-(A-(B*K_1)))*B.*G + (D*G)
H_CL = ((C)*inv((s.*eye(3))-(A-(B*K)))*B)*G;
[y_CL,t2]=step(H_CL,20);
figure(2)
step(H_CL,20);



Q = eye(3) * 5;
R = 1;
K_Lqr = lqr(A,B,Q,R)
DC_gain = 1;
G_Lqr = (-1*DC_gain)/(C*inv(A-(B*K_Lqr))*B)
H_CL_LQR = ((C)*inv((s.*eye(3))-(A-(B*K_Lqr)))*B)*G_Lqr  
figure(3)
step(H_CL_LQR,20)
[y_CLLQR,t3]=step(H_CL_LQR,20);


figure(4)
plot(t1,y_OL,'b','LineWidth',2)
hold
plot(t2,y_CL,'r','LineWidth',2) 
plot(t3,y_CLLQR,'g','LineWidth',2)
xlabel('$t$ (s)', 'Interpreter','latex')  
ylabel('$y$', 'Interpreter','latex')  
legend('Open−Loop','PID Controller','LQR')
set(gca,'linewidth',2,'fontsize',20,'fontname','Times');
set(gcf,'color','white')


% Animation
% % Airplane shape (triangle)
% airplane_x = [0, -0.5, -0.5, 0]; % Airplane body x-coordinates
% airplane_y = [0, 0.25, -0.25, 0]; % Airplane body y-coordinates
% Define airplane shape
% airplane_x = [0, -0.5, -0.7, -1, -1.2, -1, -0.7, -0.5, 0]; % X-coordinates
airplane_x = [0, -0.1, -0.7, -0.9, -0.9,  0]; % X-coordinates
airplane_x = airplane_x+0.5;
airplane_y = [0,  0.2,  0.2,  0.5,    0,  0]; % Y-coordinates


% Set up figure
figure;
axis equal;
hold on;
grid on;
xlim([-2, 2]); % X-axis range for visualization
ylim([-2, 2]); % Y-axis range for visualization
xlabel('X-axis');
ylabel('Y-axis');
title('Airplane Response to Step Input');

% Plot the initial airplane
h = fill(airplane_x, airplane_y, 'b'); % Initial airplane plot

disp("Starting animation")
% Animate the airplane
for i = 1:200%(length(y_OL))
    % Calculate the rotation matrix for pitch angle
    theta = y_OL(i) % Pitch angle from step response (in radians)
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)]; % 2D rotation matrix
    disp(i);
    % Rotate the airplane
    rotated_coords = R * [airplane_x; airplane_y];
    set(h, 'XData', rotated_coords(1, :), 'YData', rotated_coords(2, :)); % Update airplane position

    % Pause to create animation effect
    pause(0.03);
    if i == 300 
        i = length(y_OL);
    end

end


%% Closed loop system animation
% Set up figure
figure;
axis equal;
hold on;
grid on;
xlim([-2, 2]); % X-axis range for visualization
ylim([-2, 2]); % Y-axis range for visualization
xlabel('X-axis');
ylabel('Y-axis');
title('Airplane Response to Step Input');

% Plot the initial airplane
h = fill(airplane_x, airplane_y, 'b'); % Initial airplane plot

disp("Starting animation")
% Animate the airplane
for i = 1:200%(length(y_CL))
    % Calculate the rotation matrix for pitch angle
    disp(i);
    theta = y_CL(i) % Pitch angle from step response (in radians)
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)]; % 2D rotation matrix

    % Rotate the airplane
    rotated_coords = R * [airplane_x; airplane_y];
    set(h, 'XData', rotated_coords(1, :), 'YData', rotated_coords(2, :)); % Update airplane position

    % Pause to create animation effect
    pause(0.03);
    if i == 300 
        i = length(y_CL);
    end
end

%% LQR closed loop system
%% Closed loop system animation
% Set up figure
figure;
axis equal;
hold on;
grid on;
xlim([-2, 2]); % X-axis range for visualization
ylim([-2, 2]); % Y-axis range for visualization
xlabel('X-axis');
ylabel('Y-axis');
title('Airplane Response to Step Input');

% Plot the initial airplane
h = fill(airplane_x, airplane_y, 'b'); % Initial airplane plot

disp("Starting animation")
% Animate the airplane
for i = 1:200%length(y_CL)
    % Calculate the rotation matrix for pitch angle
    theta = y_CLLQR(i) % Pitch angle from step response (in radians)
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)]; % 2D rotation matrix

    % Rotate the airplane
    rotated_coords = R * [airplane_x; airplane_y];
    set(h, 'XData', rotated_coords(1, :), 'YData', rotated_coords(2, :)); % Update airplane position

    % Pause to create animation effect
    pause(0.03);
end

