function plot_deformation_global(X, connectivity, U, Ndof, scale)
figure; hold on; axis equal; grid on;
title('Deformed Structure (Global Coordinates)');
xlabel('X (m)'); ylabel('Y (m)');

num_elem = size(connectivity, 1);
nnode = size(X, 1);

% === 每个节点的全局 Δx Δy 位移 ===
U_global_nodal = zeros(nnode, 2);
for n = 1:nnode
    ux = U(Ndof(n,1));  % x-direction DOF
    uy = U(Ndof(n,2));  % y-direction DOF
    U_global_nodal(n,:) = [ux, uy];
end

for e = 1:num_elem
    n1 = connectivity(e,1);
    n2 = connectivity(e,2);
    x1 = X(n1,1); y1 = X(n1,2);
    x2 = X(n2,1); y2 = X(n2,2);

    % 原始结构
    plot([x1, x2], [y1, y2], 'k-', 'LineWidth', 1.5);

    % 位移后结构
    dx1 = U_global_nodal(n1,1) * scale;
    dy1 = U_global_nodal(n1,2) * scale;
    dx2 = U_global_nodal(n2,1) * scale;
    dy2 = U_global_nodal(n2,2) * scale;

    x_def = linspace(x1 + dx1, x2 + dx2, 11);
    y_def = linspace(y1 + dy1, y2 + dy2, 11);
    plot(x_def, y_def, 'r--', 'LineWidth', 1.5);
end

% ==== 设置边界 ====
xlimits = xlim;
ylimits = ylim;
xpad = 0.02 * (xlimits(2) - xlimits(1));
ypad = 0.02 * (ylimits(2) - ylimits(1));
xlim([xlimits(1) - xpad, xlimits(2) + xpad]);
ylim([ylimits(1) - ypad, ylimits(2) + ypad]);


legend('Original Structure', 'Deformed Shape');
end
