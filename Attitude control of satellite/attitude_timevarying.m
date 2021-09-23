% M-script for numerical integration of the attitude dynamics of a rigid 
% body represented by unit quaternions. The MSS m-files must be on your
% Matlab path in order to run the script.
%
% System:                      .
%                              q = T(q)w
%                              .
%                            I w - S(Iw)w = tau
% Control law:
%                            tau = constant
% 
% Definitions:             
%                            I = inertia matrix (3x3)
%                            S(w) = skew-symmetric matrix (3x3)
%                            T(q) = transformation matrix (4x3)
%                            tau = control input (3x1)
%                            w = angular velocity vector (3x1)
%                            q = unit quaternion vector (4x1)
%
% Author:                   2018-08-15 Thor I. Fossen and Håkon H. Helgesen

%% USER INPUTS
h = 0.1;                     % sample time (s)
N  = 4000;                    % number of samples. Should be adjusted

% model parameters
m = 180;
r = 2;
I = m*r^2*eye(3);            % inertia matrix
I_inv = inv(I);

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

phi = -5*deg2rad;            % initial Euler angles
theta = 10*deg2rad;
psi = -20*deg2rad;

q = euler2q(phi,theta,psi);   % transform initial Euler angles to q

w = [0 0 0]';                 % initial angular rates

table = zeros(N+1,14);        % memory allocation

% control params
K_d = 400 * eye(3);
K_p = 20;

% initial value
q_error = [1,0,0,0]';
w_error = [0,0,0]';

%% FOR-END LOOP
for i = 1:N+1
   t = (i-1)*h;                  % time
   
   % Reference signals
   phi_ref = 0 * deg2rad;
   theta_ref = 15*cos(0.1*t) * deg2rad;
   psi_ref = 10*sin(0.05*t) * deg2rad;
   
   phidot_ref = 0 * deg2rad;
   thetadot_ref = -1.5*sin(0.1*t) * deg2rad;
   psidot_ref = 0.5*cos(0.05*t) * deg2rad;
   
   eulerdot_ref = [phidot_ref, thetadot_ref, psidot_ref]';
   
   q_ref = euler2q( phi_ref, theta_ref, psi_ref);
   q_ref = q_ref/norm(q_ref);
   
   T_inv = [...
            1 0 -sin(theta_ref);
            0 cos(phi_ref) cos(theta_ref)*sin(phi_ref)
            0 -sin(phi_ref) cos(theta_ref)*cos(phi_ref)]; % Use inv(Tzyx) next time
        
   w_ref = T_inv * eulerdot_ref; % See (2.83) 6-DOF, Fossen 
   
   % simulation
   eps = q_error(2:4);
   tau = -K_d*w_error - K_p*eps;       % control law

   [phi,theta,psi] = q2euler(q); % transform q to Euler angles
   [J,J1,J2] = quatern(q);       % kinematic transformation matrices
   
   q_dot = J2*w;                        % quaternion kinematics
   w_dot = I_inv*(Smtrx(I*w)*w + tau);  % rigid-body kinetics
   
   table(i,:) = [t q' phi theta psi w' tau'];  % store data in table
   
   q = q + h*q_dot;	             % Euler integration
   w = w + h*w_dot;
   
   q = q/norm(q);               % unit quaternion normalization
   
   % Error
   q_error = quaternion_error(q_ref,q);
   q_error = q_error/norm(q_error);
   w_error = w - w_ref;
   
end 

%% PLOT FIGURES
t       = table(:,1);  
q       = table(:,2:5); 
phi     = rad2deg*table(:,6);
theta   = rad2deg*table(:,7);
psi     = rad2deg*table(:,8);
w       = rad2deg*table(:,9:11);  
tau     = table(:,12:14);


figure (1); clf;
hold on;
plot(t, phi, 'b');
plot(t, theta, 'r');
plot(t, psi, 'g');
hold off;
grid on;
legend('\phi', '\theta', '\psi');
title('Euler angles');
xlabel('time [s]'); 
ylabel('angle [deg]');

figure (2); clf;
hold on;
plot(t, w(:,1), 'b');
plot(t, w(:,2), 'r');
plot(t, w(:,3), 'g');
hold off;
grid on;
legend('p', 'q', 'r');
title('Angular velocities');
xlabel('time [s]'); 
ylabel('angular rate [deg/s]');

figure (3); clf;
hold on;
plot(t, tau(:,1), 'b');
plot(t, tau(:,2), 'r');
plot(t, tau(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Control input');
xlabel('time [s]'); 
ylabel('input [Nm]');

function q_diff = quaternion_error(x,y)
    x = [x(1); -x(2:4)]; % conjugate multiplied with y gives error

    nd = x(1);
    n = y(1);
    ed = x(2:4);
    e = y(2:4);
    nu = nd*n - ed'* e;
    eps = n*ed + nd*e + Smtrx(ed)*e;
    q_diff = [nu; eps];
end

