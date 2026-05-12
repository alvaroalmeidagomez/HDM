clearvars;
close all;
clc;

fprintf('\n=================================================================\n');
fprintf(' Project supported by the Centro de Modelamiento Matematico (Chile)\n');
fprintf('=================================================================\n\n');

%% ===================== USER INPUT & VALIDATION =====================
m_in = input('Enter the number of sample points on the sphere (e.g., 1000): ');

if isempty(m_in) || m_in < 3 || m_in ~= floor(m_in)
    warning('Invalid input. Defaulting to m = 1000.');
    m = 1000;
else
    m = m_in;
end

%% ===================== PARAMETER DEFINITION =====================
t = 1 / sqrt(log(log(m)));

%% ===================== GENERATE POINTS ON SPHERE =====================
rng(1); % Reproducibility

X = randn(m, 3);
X = X ./ vecnorm(X, 2, 2); % Normalize rows to lie on unit sphere

%% ===================== TANGENT STRUCTURES =====================
% NOTE: Ensure external functions are in MATLAB path

Tang_M      = tangvect(X);
[A_mat, D]  = compA(X, Tang_M, t);

%% ===================== OPERATORS =====================
inv_t2      = 1 / (t^2);

ED_matrix   = inv_t2 * ED1(Tang_M, A_mat, D);
Grad_matrix = inv_t2 * grad_fun(A_mat, D);

%% ===================== INITIAL VECTOR FIELD =====================
vF_test = VF(X, Tang_M);

visualize_tangent_field(X, vF_test, Tang_M, 'Initial Field');

%% ===================== HODGE LAPLACIAN =====================
Hodge_Matrix = (ED_matrix.' * ED_matrix) + (Grad_matrix * Grad_matrix.');

%% ===================== NUMERICAL SOLUTIONS OF THE HEAT EQUATION =====================

% -------------------------------------------------------------------------
% PARAMETERS
% -------------------------------------------------------------------------
NUM_STEPS       = 4;
TIME_MULTIPLIER = 2;

% -------------------------------------------------------------------------
% COMPUTE HEAT SOLUTIONS
% -------------------------------------------------------------------------
heat_sol = cell(1, NUM_STEPS);

for i = 1:NUM_STEPS
    iteration = TIME_MULTIPLIER * i;
    heat_sol{i} = hodgeHeatEuler(vF_test, Hodge_Matrix, iteration);
end

% -------------------------------------------------------------------------
% VISUALIZATION
% -------------------------------------------------------------------------
visualize_heat_magnitude(heat_sol, Tang_M, X, TIME_MULTIPLIER);


%% ===================== ERROR ANALYSIS ACROSS ITERATIONS =====================
% =========================================================================
% HEAT EQUATION ? MAGNITUDE ERROR ANALYSIS
% Tracks the L2 magnitude decay of the Hodge heat flow over Euler steps.
% =========================================================================

%% --- Configuration ------------------------------------------------------

NUM_STEPS = 100;   % Number of Euler iterations

%% --- Initial Vector Field -----------------------------------------------

% Evaluate vector field on the mesh
vF_test = VF(X, Tang_M);

%% --- Preallocate Error Array ---------------------------------------------

heatMagnitudeError = zeros(1, NUM_STEPS);

%% --- Main Loop: Magnitude Error Over Iterations -------------------------

for k = 1:NUM_STEPS

    % Advance heat equation one Euler step
    currentSolution = hodgeHeatEuler(vF_test, Hodge_Matrix, k);

    % L2 magnitude error (normalized by mesh size)
    heatMagnitudeError(k) = norm(currentSolution) / m_in;

end

%% --- Plot: Heat Magnitude Error -----------------------------------------

figure('Color', 'w', 'Position', get(0, 'ScreenSize'), ...
       'Name', 'Heat Magnitude Error');

semilogy(1:NUM_STEPS, heatMagnitudeError, '-o', ...
    'LineWidth', 1.5, 'MarkerSize', 4);

grid on;  box on;
xlabel('Iteration', 'FontSize', 11);
ylabel('L^2 Error', 'FontSize', 11);
title('Heat Magnitude Error', 'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'FontSize', 10);

%% --- Hodge Distance -----------------------------------------------------

hodge_distance(X, Hodge_Matrix, m_in);