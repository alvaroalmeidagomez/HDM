function visualize_tangent_field(X, vF_test, Tang_M, stage_name)
    % ---- Magnitude and 3D tangent vectors ----
    mag       = vecnorm(reshape(vF_test, 2, []), 2, 1);
    Vecttang  = Coor2Vect(vF_test, Tang_M);
    % ---- Normalize 3D vectors ----
    Vecttang_norm = normalize_vectors_safe(Vecttang);

    figure('Name', sprintf('%s - All Views', stage_name), 'Color', 'w');
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    % ---- Subplot 1: 3D vector field ----
    subplot(2, 2, 1);
    quiver3(X(:,1), X(:,2), X(:,3), ...
            Vecttang_norm(:,1), Vecttang_norm(:,2), Vecttang_norm(:,3), ...
            0.5, 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    axis equal;
    grid on;
    view(3);
    title(sprintf('%s Normalized Tangent Vector Field', stage_name));
    xlabel('X');
    ylabel('Y');
    zlabel('Z');

    % ---- Subplot 2: 3D magnitude ----
    subplot(2, 2, 2);
    scatter3(X(:,1), X(:,2), X(:,3), 30, mag, 'filled');
    colorbar;
    title(sprintf('%s Point Magnitude Distribution', stage_name));
    axis equal;
    grid on;
    view(3);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');

    % ---- Projection to square ----
    [UV, UV_vec] = sphere_to_square_with_tangent(X, Vecttang);
    UV_vec_norm = normalize_vectors_safe(UV_vec);

    % ---- Subplot 3: 2D magnitude ----
    subplot(2, 2, 3);
    scatter(UV(:,1), UV(:,2), 30, mag, 'filled');
    colorbar;
    title(sprintf('%s Projected Magnitude', stage_name));
    axis equal;
    grid on;
    xlabel('u');
    ylabel('v');

    % ---- Subplot 4: 2D vector field ----
    subplot(2, 2, 4);
    quiver(UV(:,1), UV(:,2), UV_vec_norm(:,1), UV_vec_norm(:,2), 0.7);
    axis equal;
    grid on;
    title(sprintf('%s Normalized Projected Vector Field', stage_name));
    xlabel('u');
    ylabel('v');
end