function visualize_heat_magnitude(heat_sol, Tang_M, X, TIME_MULTIPLIER)

    NUM_STEPS = numel(heat_sol);

    % ---------------------------------------------------------------------
    % Compute global magnitude range
    % ---------------------------------------------------------------------
    mag_all    = cell(1, NUM_STEPS);

    global_min = inf;
    global_max = -inf;

    for i = 1:NUM_STEPS

        vF_current = heat_sol{i};

        % Magnitude of vector field
        mag = vecnorm(reshape(vF_current,2,[]), 2, 1);

        mag_all{i} = mag;

        global_min = min(global_min, min(mag));
        global_max = max(global_max, max(mag));
    end

    % ---------------------------------------------------------------------
    % Create figure
    % ---------------------------------------------------------------------
    figure('Name', 'Magnitude of Heat Equation Solutions', ...
           'NumberTitle', 'off', ...
           'Color', 'w', ...
           'Units', 'normalized', ...
           'OuterPosition', [0.02 0.02 0.95 0.95]);

    % ---------------------------------------------------------------------
    % Plot each solution
    % ---------------------------------------------------------------------
    for i = 1:NUM_STEPS

        iteration = TIME_MULTIPLIER * i;

        vF_current = heat_sol{i};

        % Tangent vectors
        Vecttang = Coor2Vect(vF_current, Tang_M);

        % Projection to square coordinates
        [UV, ~] = sphere_to_square_with_tangent(X, Vecttang);

        subplot(2,2,i);

        scatter(UV(:,1), ...
                UV(:,2), ...
                35, ...
                mag_all{i}, ...
                'filled');

        axis equal tight;

        xlabel('U', 'FontSize', 11);
        ylabel('V', 'FontSize', 11);

        title(sprintf('Iteration %d', iteration), ...
              'FontSize', 13, ...
              'FontWeight', 'bold');

        caxis([global_min global_max]);

        colorbar;

        grid on;
    end

    % ---------------------------------------------------------------------
    % Global title
    % ---------------------------------------------------------------------
    annotation('textbox', ...
               [0 0.95 1 0.05], ...
               'String', ...
               'Magnitude of the Heat Equation Solutions', ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'FontSize', 18, ...
               'FontWeight', 'bold');
end


