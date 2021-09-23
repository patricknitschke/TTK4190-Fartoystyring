%% Problem 1.1
%Constants
m = 180;
R_33 = 2.0;
I_g = m*R_33^2 * eye(3);
tau = sym('t', [3,1]);

%States
q = sym('q', [4,1]);
e = q(2:4);
w = sym('w', [3,1]);
x = [e;w];

%Eq1
f1 = Tquat(q)*w;
q_dot = f1(2:4);

%Eq2
Iw = (I_g * w);
S_Iw = [0 -Iw(3) Iw(2) ; 
        Iw(3) 0 -Iw(1) ; 
        -Iw(2) Iw(1) 0 ];
  
w_dot = (I_g) \ (tau + S_Iw * w); % I_inv * (...)

%Linearizations
A = [jacobian(q_dot, x);
    jacobian(w_dot, x)]

B = [jacobian(q_dot, tau); 
    jacobian(w_dot, tau)]

%%Problem 1.2
k_d = 40;
k_p = 2;

K_d = k_d * eye(3);

K = [k_p*eye(3) K_d];
tau = -K*x;

%System stability with state feedback given by eig(A-BK)
%around equilibrium
A_BK = [eye(3)*0, eye(3)*0.5;
       -eye(3)*k_p/720, -eye(3)*k_d/720]
eig(A_BK)



