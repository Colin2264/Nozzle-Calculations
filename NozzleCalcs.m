clear
clc

% Given values
P_c = 5e6;           % Chamber pressure (Pa)
T_c = 3500;          % Chamber temperature (K)
R = 377;             % Specific gas constant (J/kg*K)
gamma = 1.3;         % Ratio of specific heats
mass_flow_rate = 8;  % Mass flow rate (kg/s)
L_star = 1.35;       % Characteristic length (m)
Expansion_ratio = 6.63;

% Calculate Temperature in the throat
T_t = T_c * (2 / (gamma + 1));

% Calculate speed of sound at the throat
a_t = sqrt(gamma * R * T_t);

% Calculate Pressure at the throat
P_t = P_c * (2 / (gamma + 1)) ^ (gamma / (gamma - 1));

% Calculate density at the throat
rho_t = P_t / (R * T_t);

% Calculate Area of the throat
A_t = mass_flow_rate / (rho_t * a_t);

% Calculate Area at the exit
A_e = A_t * Expansion_ratio;

% Calculate combustion chamber density
rho_c = P_c / (R * T_c);

% R_e and R_t
R_e=sqrt(A_e/pi);
R_t=sqrt(A_t/pi);

% Length of exit line
L_e=(2*R_e-2*R_t)/(2*tand(15));

% Radius of Curvatures
R1=0.5*R_t;
R2=1.5*R_t;
R3=1.5*R_t;

% Initial guess for A_c and v_c
A_c_guess = A_t;          % Start with throat area as an initial guess
v_c_guess = mass_flow_rate / (rho_c * A_c_guess);  % Use the guess for A_c to estimate v_c

% Initial guess for L_c
L_c_guess = 1.0;  % Start with an arbitrary guess for L_c

% Set a tolerance for the iterative process
tolerance = 1e-6;
max_iterations = 100;
iterations = 0;

% Iterative process to adjust A_c, v_c, and L_c
while iterations < max_iterations
    iterations = iterations + 1;
    
    % Calculate chamber volume V_c based on A_c and v_c
    V_c = mass_flow_rate / (rho_c * v_c_guess);
    
    % Update A_c using mass flow rate and throat conditions (continuity equation)
    A_c_new = mass_flow_rate / (rho_c * v_c_guess);
    
    % Calculate L_c from the updated A_c
    L_c_new = (L_star * A_c_new) / (V_c / A_c_new);
    
    % Check if the new value for L_c is within tolerance
    if abs(L_c_new - L_c_guess) < tolerance
        break;
    end
    
    % Update A_c, v_c, and L_c for the next iteration
    A_c_guess = A_c_new;
    v_c_guess = mass_flow_rate / (rho_c * A_c_guess);
    L_c_guess = L_c_new;
end

R_c=sqrt(A_c_new/pi);

%Adjusted CC values based on typical engine ratios
R_c=1.2*R_t;
L_c=3*L_e;
A_c=pi*(R_c^2);

% Display results
fprintf('Chamber Area (A_c): %.4f m^2\n', A_c);
fprintf('Chamber Velocity (v_c): %.4f m/s\n', v_c_guess);
fprintf('Chamber Length (L_c): %.4f m\n', L_c);
fprintf('At: %d\n', A_t);
fprintf('Rt: %d\n', R_t);
fprintf('Ae: %d\n', A_e);
fprintf('Re: %d\n', R_e);
fprintf('Rc: %d\n', R_c);
fprintf('R1: %d\n', R1);
fprintf('R2 and R3: %d\n', R2);
fprintf('Exit Line Length: %d\n', L_e);


