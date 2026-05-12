function v_out = hodgeHeatEuler(v0, A, num_iters)
    % HODGEHEATEULER Solves the heat equation using the Hodge Laplacian
    % via the explicit Euler method.
    %
    % Inputs:
    %   v0        : Initial state vector (column vector).
    %   A         : Hodge Laplacian matrix (should be a sparse matrix).
    %   num_iters : Number of explicit Euler iterations.
    %
    % Output:
    %   v_out     : Final state vector after num_iters.

    % Ensure the input is a column vector
    if ~iscolumn(v0)
        v0 = v0(:);
    end

    % 1. OPTIMIZED EIGENVALUE COMPUTATION
    % The original 'eig(A)' computes ALL eigenvalues. For an N x N matrix, 
    % this is an O(N^3) operation and will run out of memory for large meshes.
    % We only need the largest eigenvalue, so we use 'eigs' which is much faster.
    opts.disp = 0; % Suppress command window output from eigs
    try
        % 'lm' asks for the single Largest Magnitude eigenvalue
        max_eig = eigs(A, 1, 'lm', opts); 
    catch
        % Fallback just in case the matrix is extremely small (e.g., < 6x6)
        % where 'eigs' might fail to converge.
        max_eig = max(eig(full(A)));
    end

    % Calculate stable time step.
    % The absolute mathematical limit is 2/max_eig. We use 1.99 to give 
    % a tiny safety margin against floating-point precision drift.
    dt = 1.9 / max_eig;

    % Initialize output vector
    v_out = v0;

    % 2. OPTIMIZED LOOP PERFORMANCE
    for i = 1:num_iters
        % CRITICAL FIX: Added parentheses around (A * v_out).
        % 
        % Why this matters:
        % Without parentheses, MATLAB evaluates 'dt * A * v_out' left-to-right.
        % It computes (dt * A) first, which forces MATLAB to allocate memory 
        % for a brand new, scaled matrix during EVERY SINGLE LOOP ITERATION.
        % 
        % By using 'dt * (A * v_out)', MATLAB performs matrix-vector multiplication 
        % first (which outputs a small vector), and then scales that vector by dt. 
        % This simple change makes the loop significantly faster.
        v_out = v_out - dt * (A * v_out);
    end
end