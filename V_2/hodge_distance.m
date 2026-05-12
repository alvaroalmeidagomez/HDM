function hodge_distance(X, HodgeMatrix, m_in)




    % Preallocate output
    hodge_distance = zeros(m_in,1);

    % Size of each block
    blockSize = 2;

    % Frobenius norm of the first diagonal block
    block11 = HodgeMatrix(1:blockSize, 1:blockSize);
    var1_sq = norm(block11, 'fro')^2;

    % Compute Hodge distances
    for i = 1:m_in

        idx = (i-1)*blockSize + (1:blockSize);

        Hii = HodgeMatrix(idx, idx);
        H1i = HodgeMatrix(1:blockSize, idx);

        var2_sq = norm(Hii, 'fro')^2;
        var3_sq = norm(H1i, 'fro')^2;

        hodge_distance(i) = var1_sq + var2_sq - var3_sq;
    end

    % ---- Map sphere to square ----
    x = X(:,1);
    y = X(:,2);
    z = X(:,3);

    theta = atan2(y,x);
    phi   = acos(z);

    u = theta/(2*pi) + 0.5;
    v = phi/pi;

    UV = [u v];

    % ================= FIGURE =================
    figure('Name','Hodge distance','NumberTitle','off');

    % Expand to full screen (MATLAB 2017b compatible)
    set(gcf,'Units','normalized','OuterPosition',[0 0 1 1]);

    % -------- Subplot 1: 2D map --------
    subplot(1,2,1)

    plot(UV(1,1), UV(1,2), 'rx', 'MarkerSize', 20, 'LineWidth', 4);
    hold on
    scatter(UV(:,1), UV(:,2), 30, hodge_distance, 'filled');

    colorbar;
    title('Hodge distance in 2D coordinate system');
    xlabel('u'); ylabel('v');
    axis tight
    grid on

    % -------- Subplot 2: 3D sphere --------
    subplot(1,2,2)
    
    scatter3(X(1,1), X(1,2), X(1,3), 100, 'rx', 'LineWidth', 10);
    hold on
    scatter3(X(:,1), X(:,2), X(:,3), 30, hodge_distance, 'filled');

    colorbar;
    title('Hodge distance in 3D unit sphere');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on
    axis equal
    view(3)

end
