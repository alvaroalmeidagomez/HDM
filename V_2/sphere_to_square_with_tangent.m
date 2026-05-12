function [UV, UV_vec] = sphere_to_square_with_tangent(X, V)

    m = size(X,1);

    x = X(:,1);
    y = X(:,2);
    z = X(:,3);

    % ---- Spherical coordinates ----
    theta = atan2(y,x);
    phi   = acos(z);

    % ---- Map to unit square ----
    u = theta/(2*pi) + 0.5;
    v = phi/pi;

    UV = [u v];

    % ---- Tangent basis ----
    sinphi = sin(phi);
    cosphi = cos(phi);
    costh  = cos(theta);
    sinth  = sin(theta);

    dtheta = [-sinth.*sinphi, costh.*sinphi, zeros(m,1)];
    dphi   = [ costh.*cosphi, sinth.*cosphi, -sinphi];

    UV_vec = zeros(m,2);

    for i = 1:m

        B = [dtheta(i,:)' dphi(i,:)'];

        coeff = B \ V(i,:)';

        a = coeff(1);
        b = coeff(2);

        u_dot = a/(2*pi);
        v_dot = b/pi;

        UV_vec(i,:) = [u_dot v_dot];

    end

end