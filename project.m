clc; clear; close all;

% System constants
mb = 0.2445;
mw = 0.0343;
R = 0.035;
l = 0.035;
Jb = 0.5 * mb * l^2;
Jw = 0.5 * mw * R^2;
g = 9.81;

% function calls for system design and evaluation
[sys_ol, A_ol, B_ol, C_ol, D_ol] = system_definition(mb, mw, R, l, Jb, Jw, g);
system_analysis(sys_ol);

[Co, Ob] = control_observe(A_ol, B_ol, C_ol);

[sys_cl, A_cl, K_cl] = state_feedback_stabilizer_design(A_ol, B_ol, C_ol, D_ol, sys_ol);
system_analysis(sys_cl);

[A_obs, B_obs, C_obs, D_obs, sys_obs] = observer_design(A_ol, B_ol, C_ol, D_ol);

%% FUNCTIONS
function [sys_ol, A, B, C, D] = system_definition(mb, mw, R, l, Jb, Jw, g)

    % Linearized state space system definition

    c1 = (mb^2 * l^2 * g)/((mb + 2 * (mw + (Jw/R^2) - (mb^2 * l^2)/(mb * l^2 + Jb))*(mb^2 * l^2 * g)/(mb * l^2 + Jb)));
    c2 = (2/(mb + 2 * (mw + Jw/R^2) - (mb^2 * l^2)/(mb * l^2 + Jb))) * (1/R + (mb * l)/(mb * l^2 + Jb));
    c3 = (mb * g * l) * (1/Jb + 1/(mb * l^2) - ((mb + 2 * (mw + Jw/R^2)/(mb^2 * l^2))));
    c4 = (2/R) * (1/Jb + 1/(mb * l^2) - ((mb + 2 * (mw + Jw/R^2)/(mb^2 * l^2)))) * (R + ((mb * l)/(mb + 2 * (mw + Jw/R^2))));
    
    A = [0 1 0 0; 0 0 -c1 0; 0 0 0 1; 0 0 c3 0];
    B = [0; c2; 0; -c4];
    C = eye(4);
    D = [0; 0; 0; 0];
    
    sys_ol = ss(A,B,C,D);
end

function system_analysis(sys)
    % Step response and Pole-Zero map of linearized system
    figure()
    step(sys);

    figure()
    pzmap(sys);
end

function [Co, Ob] = control_observe(A, B, C)
    Co = ctrb(A,B);
    Ob = obsv(A,C);
    disp("Open loop system is:")
    if rank(Co) == length(A)
        disp("Controllable")
    else
        disp("Not Controllabe")
    end
    if rank(Ob) == length(A)
        disp("Observable")
    else
        disp("Not Observable")
    end
end

function [sys_cl, A_cl, K] = state_feedback_stabilizer_design(A, B, C, D, sys_ol)
    P = [-3 -6 -30 -50];
    K = place(A,B,P);
    Q = diag([1 1 3 1]);
    % [K, S, poles_cl] = lqr(sys_ol, Q, 2.0)
    A_cl = A-B*K;
    
    sys_cl = ss(A_cl,B,C,D);

end

function [A_obs, B_obs, C_obs, D_obs, sys_obs] = observer_design(A, B, C, D)
    poles_obs = [-15 -30 -150 -250];
    L = place(A',C',poles_obs).';

    A_obs = A-L*C;
    B_obs = [B L];
    C_obs = eye(4); 
    D_obs = [D D D D D];
    sys_obs = ss(A_obs,B_obs,C_obs,D_obs);
end
