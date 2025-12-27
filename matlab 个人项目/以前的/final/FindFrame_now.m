function FindFrame_now

%Structural Mechanics
%Matrix analysis of frame
%Program created by Dr. Antonio Gil
%Program Developed by Dr. Oubay Hassan and Shen Ma  2024.6 - 2025.2


clc;close all;

datafilename = input('Enter the input file name:');

% =========================================================================
% Data input phase. Read information from input data file
% =========================================================================
[ndof,X,connectivity,Ndof,Nelem,totaldof,E,A,I,H,F,U,alldof,Elastic,Bnode,w,forceunits,lengthunits,UDLFinternal]=Read_input(datafilename);


%assemble global stiffness matrix

[K,T1]=compute_stiffness_matrix(ndof,X,Ndof,connectivity,E,A,I,totaldof,Nelem);

% =========================================================================
%  Modify matrix phase. Modify matrix according boundary conditions
% =========================================================================

%index fix and free degrees of freedom

fixdof  = find(alldof == 1);
freedof = find(alldof == 0);
%modify the K matrix according to fixed support
K_mod = K;
F_mod = F;
T1_mod = T1;
%zero the coefficients of the  corresponding line and column
K_mod(fixdof, :) = 0;
K_mod(:, fixdof) = 0;

%place 1 in the diagonal location of the same equation.
for iq = 1:length(fixdof)
    K_mod(fixdof(iq), fixdof(iq)) = 1;
end
%modify the F matrix according to prescribed displacement
% Now U1 is all prescribed displacement
U1 = U;

dpredof = find(U == 0);
predof  = find(U ~= 0);

if predof ~= 0
    for iu = 1:length(predof)
        for iux = 1:length(freedof)
            % F_mod(dpredof(iux)) is degree of freedom no prescribed displacement
            F_mod(freedof(iux)) = F_mod(freedof(iux)) - U(predof(iu))*K(freedof(iux),predof(iu));
            F_mod(predof)  = U(predof);
        end
    end
else
end

%Modify the K matrix according to elastic supports
El = Elastic(fixdof);
% Place elastic support value in the diagonal location of the same equation.
for iel = 1:length(fixdof)
    K_mod(fixdof(iel), fixdof(iel)) = K_mod(fixdof(iel), fixdof(iel))+El(iel);
end
% =========================================================================
%  Computing unknown value phase. Computing unknown displacements, final
%  displacements, unknown reactions, internal forces
% =========================================================================

%computing unknown displacements

U = K_mod \ F_mod;
U(fixdof) = 0;
% unknown displacements + prescribed displacement
U2 = U+U1;
F_mod = K_mod * U2;
R = K(fixdof,:) * U2 - F(fixdof);
U = U2;
UFINAL = U;
Ud = zeros(size(U));

F = F_mod;
F(fixdof) = R;

%computing internal forces
[Finternal]=compute_internal_forces(Nelem,X,connectivity,E,A,I,U2,UDLFinternal);


% ====== Menu System ======
fprintf('\nFrame analysis complete.\n');

while true
    disp('---------------------------------------------');
    disp('What do you want to do?');
    disp('1 - Plot Bending Moment Diagram');
    disp('2 - Plot Shear   Force  Diagram');
    disp('3 - Plot  Deformation   Diagram');
    disp('4 - Display Displacements results');
    disp('5 - Display Node Forces');
    disp('6 - Output result to file');
    disp('7 - Element Information');
    disp('8 - Exit');
    choice = input('Select option: ');

    switch choice

        %-------- Case 1: Plot Bending Moment Diagram --------


        case 1
            disp('Select Bending Moment diagram scale:');
            disp('1 - 0.2 (large display)');
            disp('2 - 0.0001 (small display)');
            disp('3 - 0.005 (around 500 KN)');
            scale_choice = input('Select option [1/3]: ');
            if scale_choice == 1
                scale = 0.2;
            elseif scale_choice == 2
                scale = 0.0001;
            else
                scale = 0.005;
            end

            plot_BM_only(X, connectivity, Finternal, w, scale, true)


            %-------- Case 2: Plot Shear Force Diagram --------
        case 2
            disp('Select Shear Force diagram scale:');
            disp('1 - 0.2 (large display)');
            disp('2 - 0.0001 (small display)');
            disp('3 - 0.005 (around 500 KN)');
            scale_choice = input('Select option [1/3]: ');
            if scale_choice == 1
                scale = 0.2;
            elseif scale_choice == 2
                scale = 0.0001;
            else
                scale = 0.005;
            end
            plot_SF_only(X, connectivity, Finternal, w, scale, true);

            %-------- Case 3: Plot Global Deformation Diagram --------
        case 3
            disp('1 - 20 (large display)');
            disp('2 - 1000 (small display)');
            disp('3 - 5000 (around 500 KN)');
            scale_choice = input('Select option [1/3]: ');
            if scale_choice == 1
                scale = 20;

            elseif scale_choice == 2
                scale = 1000;
            else
                scale = 5000;
            end
            [~, ~, TotalGlobal] = plot_deformation_local(X, connectivity, U, Ud, Nelem, w, E, I, A, Finternal, fixdof, UFINAL);
            plot_deformation_global(X, connectivity, UFINAL, Ndof, TotalGlobal, scale)

            %-------- Case 4: Output Node Displacements --------
        case 4
            disp('Display displacements in');
            disp('1 - meters (m)');
            disp('2 - millimeters (mm)');
            unit_choice = input('Select option [1/2]: ');

            if unit_choice == 2
                scale_disp = 1e3;
                unit_str = 'mm';
            else
                scale_disp = 1;
                unit_str = 'm';
            end

            fprintf('Node Displacements [U (%s)       V (%s)     Theta (rad)]:\n\n', unit_str, unit_str);
            for e1 = 1:size(connectivity,1)
                n13 = connectivity(e1,1);
                n23 = connectivity(e1,2);
                U3 = U(Nelem(:,:,e1));
                fprintf('Element %d\n', e1);
                fprintf('  Node %d       %12.5g %12.5g %12.5g\n', n13, U3(1), U3(2), U3(3));
                fprintf('  Node %d       %12.5g %12.5g %12.5g\n\n', n23, U3(4), U3(5), U3(6));
            end

            %-------- Case 5: Output Reaction Forces --------
        case 5
            disp('Reaction Forces       [ Fx (KN)   Fy (KN)    M (KN·m)]');
            disp('----------------------------------------------------------');
            Fcase4 = F;
            Fcase4(abs(Fcase4) < 1e-5) = 0;
            for i1 = 1:length(Bnode)
                node_id = Bnode(i1);
                idof1 = Ndof(node_id, :);
                coord = X(node_id, :);
                fprintf('Node%-3d(%4.1f, %4.1f)  %10.3f %10.3f %10.3f\n', ...
                    node_id, coord(1), coord(2), ...
                    Fcase4(idof1(1)), Fcase4(idof1(2)), Fcase4(idof1(3)));
            end

            %-------- Case 6: Output results --------
        case 6
            [~, name, ~] = fileparts(datafilename);
            outputname = ['Output_' name '.txt'];
            fid = fopen(outputname, 'w');

            %% -------- Element-wise Nodal Displacements --------
            fprintf(fid, '========== Nodal Displacements ==========\n');
            fprintf(fid, 'Units: mm (U, V), rad (θ)\n\n');

            for e1 = 1:size(connectivity,1)
                n13 = connectivity(e1,1);
                n23 = connectivity(e1,2);
                U3 = U(Nelem(:,:,e1));
                U3(abs(U3) < 1e-6) = 0;
                fprintf(fid, 'Element %d\n', e1);
                fprintf(fid, '  Node %d   %12.5g %12.5g %12.5g\n', n13, U3(1)*1e3, U3(2)*1e3, U3(3));
                fprintf(fid, '  Node %d   %12.5g %12.5g %12.5g\n\n', n23, U3(4)*1e3, U3(5)*1e3, U3(6));
            end

            %% --------  Reaction Forces --------
            fprintf(fid, '\n========== Reaction Forces ==========\n');
            fprintf(fid, 'Units: KN, KN·m\n\n');
            fprintf(fid, 'Node   (X, Y) [mm]           Fx (KN)    Fy (KN)     M (KN·m)\n');
            fprintf(fid, '--------------------------------------------------------------\n');

            Fcase4 = F;
            Fcase4(abs(Fcase4) < 1e-5) = 0;

            for i1 = 1:length(Bnode)
                node_id = Bnode(i1);
                idof1 = Ndof(node_id, :);
                coord = X(node_id, :) * 1e3;

                fprintf(fid, ' %-4d (%6.1f, %6.1f)    %10.3f %10.3f %10.3f\n', ...
                    node_id, coord(1), coord(2), ...
                    Fcase4(idof1(1)), Fcase4(idof1(2)), Fcase4(idof1(3)));
            end

            %% --------  Element Internal Info --------
            fprintf(fid, '\n========== Element Internal Forces & Displacement ==========\n');
            fprintf(fid, 'Units: Position (m), BM (KN·m), SF (KN), Disp (mm)\n\n');

            BM_all = plot_BM_only(X, connectivity, Finternal, w, 0, false);
            SF_all = plot_SF_only(X, connectivity, Finternal, w, 0, false);
            [TotalYd, TotalXd, TotalGlobal] = plot_deformation_local(X, connectivity, U, Ud, Nelem, w, E, I, A, Finternal, fixdof);
            axial_all = Finternal(:,2,:);

            GlobalX = squeeze(TotalGlobal(:,1,:));
            GlobalY = squeeze(TotalGlobal(:,2,:));

            for e6 = 1:size(connectivity, 1)
                n51 = connectivity(e6,1); n52 = connectivity(e6,2);
                x1c = X(n51,1); y51 = X(n51,2);
                x2c = X(n52,1); y52 = X(n52,2);
                L6 = norm([x2c - x1c, y52 - y51]);
                x_local6 = linspace(0, L6, 11);

                fprintf(fid, 'Element %d (Length: %.3f m)\n', e6, L6);
                fprintf(fid, '%-15s %-15s %-15s %-15s %-15s %-15s\n', ...
                    'Position (m)', 'AF (KN)', 'BM (KN·m)', 'SF (KN)', 'X Disp (mm)', 'Y Disp (mm)');

                for k5 = 1:length(x_local6)
                    xd = GlobalX(k5,e6) * 1e3;
                    yd = GlobalY(k5,e6) * 1e3;
                    if abs(xd) < 1e-6, xd = 0; end
                    if abs(yd) < 1e-6, yd = 0; end

                    fprintf(fid, '%-15.3f %-15.3f %-15.5f %-15.5f %-15.5f %-15.5f\n', ...
                        x_local6(k5), axial_all(e6), BM_all(k5,e6), SF_all(k5,e6), xd, yd);
                end
                fprintf(fid, '\n');
            end


            %% -------- Close file --------
            fclose(fid);
            fprintf('Results written to %s\n', outputname);


            %-------- Case 8: Exit --------

        case 8
            disp('Exit.');
            break;


            %-------- Case 7: Element Information --------

        case 7
            elem_id = input('Enter element number to inspect: ');
            if elem_id < 1 || elem_id > size(connectivity,1)
                disp('Invalid element number.');
                continue;
            end

            % Node and geometry info
            n1_info = connectivity(elem_id,1);
            n2_info = connectivity(elem_id,2);
            x1_info = X(n1_info,1); y1_info = X(n1_info,2);
            x2_info = X(n2_info,1); y2_info = X(n2_info,2);
            L_info = norm([x2_info - x1_info, y2_info - y1_info]);

            fprintf('\nElement %d connects:\n', elem_id);
            fprintf('  Node %d at (%.3f, %.3f)\n', n1_info, x1_info, y1_info);
            fprintf('  Node %d at (%.3f, %.3f)\n', n2_info, x2_info, y2_info);
            fprintf('  Length = %.3f m\n', L_info);

            x_local_info = linspace(0, L_info, 11);

            BM_all = plot_BM_only(X, connectivity, Finternal, w, 0, false);
            SF_all = plot_SF_only(X, connectivity, Finternal, w, 0, false);
            [TotalYd, TotalXd, TotalGlobal] = plot_deformation_local(X, connectivity, U, Ud, Nelem, w, E, I, A, Finternal, fixdof);


            GlobalX = squeeze(TotalGlobal(:,1,:));

            GlobalY = squeeze(TotalGlobal(:,2,:));


            % Display table
            fprintf('\n%-15s %-15s %-15s %-15s %-15s\n', ...
                'Position (m)', 'BM (KN·m)', 'SF (KN)', 'X Disp  (mm)','Y Disp  (mm)');
            for i11 = elem_id
                for k_info = 1:length(x_local_info)
                    fprintf('%-15.3f %-15.5f %-15.5f %-15.5f %-15.5f\n', ...
                        x_local_info(k_info), BM_all(k_info,i11), SF_all(k_info,i11), GlobalX(k_info,i11)*1e3, GlobalY(k_info,i11)*1e3);
                end
            end

        otherwise
            disp('Invalid option.');
    end
end


%plot structure displacements
% =========================================================================================================================================
%  Sub Function
% =========================================================================================================================================


    function [K,T1]=compute_stiffness_matrix(ndof,X,Ndof,connectivity,E,A,I,totaldof,Nelem)
        T1 = zeros(totaldof,totaldof);
        K = zeros(totaldof,totaldof);
        %loop over elements
        for ielem=1:size(connectivity,1)
            %extract connectivity information
            i=connectivity(ielem,1);
            j=connectivity(ielem,2);
            idof = Nelem(:,1:3,ielem);
            jdof = Nelem(:,4:6,ielem);
            %compute geometric information
            Xij=X(j,:)-X(i,:);
            L=norm(Xij,2); c=Xij(1)/L; s=Xij(2)/L;
            %frame transformation matrix
            T=[c s 0; -s c 0; 0 0 1];
            %compute element local stiffness blocks
            kii=[E(ielem)*A(ielem)/L               0                            0;...
                0             12*E(ielem)*I(ielem)/L^3     6*E(ielem)*I(ielem)/L^2;...
                0              6*E(ielem)*I(ielem)/L^2     4*E(ielem)*I(ielem)/L];

            kij=[-E(ielem)*A(ielem)/L,               0,                            0;...
                0,             -12*E(ielem)*I(ielem)/L^3,     6*E(ielem)*I(ielem)/L^2;...
                0,             - 6*E(ielem)*I(ielem)/L^2,     2*E(ielem)*I(ielem)/L];

            kji=[-E(ielem)*A(ielem)/L,               0,                            0;...
                0,             -12*E(ielem)*I(ielem)/L^3,    -6*E(ielem)*I(ielem)/L^2;...
                0,               6*E(ielem)*I(ielem)/L^2,     2*E(ielem)*I(ielem)/L];

            kjj=[E(ielem)*A(ielem)/L,               0,                             0;...
                0,             12*E(ielem)*I(ielem)/L^3,     -6*E(ielem)*I(ielem)/L^2;...
                0,             -6*E(ielem)*I(ielem)/L^2,      4*E(ielem)*I(ielem)/L];

            %compute element global stiffness blocks
            Kii=T'*kii*T;Kij=T'*kij*T;Kji=T'*kji*T;Kjj=T'*kjj*T;
            %assemble global stiffness matrix
            K(idof,idof)=K(idof,idof)+Kii; K(idof,jdof)=K(idof,jdof)+Kij;
            K(jdof,idof)=K(jdof,idof)+Kji; K(jdof,jdof)=K(jdof,jdof)+Kjj;

            T1(idof,idof) = T1(idof,idof)+T; T1(idof,jdof) = T1(idof,jdof)+T;
            T1(jdof,idof) = T1(jdof,idof)+T; T1(jdof,jdof) = T1(jdof,jdof)+T;

        end
    end


    function [Finternal]=compute_internal_forces(Nelem,X,connectivity,E,A,I,U2,UDLFinternal)
        %initial allocation
        Finternal=zeros(size(connectivity,1),2,3);

        %loop over elements
        for ielem=1:size(connectivity,1)
            %extract connectivity information
            i=connectivity(ielem,1);
            j=connectivity(ielem,2);
            idof = Nelem(:,1:3,ielem);
            jdof = Nelem(:,4:6,ielem);
            %compute geometric information
            Xij=X(j,:)-X(i,:);
            L=norm(Xij,2); c=Xij(1)/L; s=Xij(2)/L;
            %frame transformation matrix
            T=[c s 0; -s c 0; 0 0 1];
            %local displacements
            ui=U2(idof,1);uj=U2(jdof,1);
            %compute element local stiffness blocks
            kii=[E(ielem)*A(ielem)/L               0                            0;...
                0             12*E(ielem)*I(ielem)/L^3     6*E(ielem)*I(ielem)/L^2;...
                0              6*E(ielem)*I(ielem)/L^2     4*E(ielem)*I(ielem)/L];

            kij=[-E(ielem)*A(ielem)/L,               0,                            0;...
                0,             -12*E(ielem)*I(ielem)/L^3,     6*E(ielem)*I(ielem)/L^2;...
                0,             - 6*E(ielem)*I(ielem)/L^2,     2*E(ielem)*I(ielem)/L];

            kji=[-E(ielem)*A(ielem)/L,               0,                            0;...
                0,             -12*E(ielem)*I(ielem)/L^3,    -6*E(ielem)*I(ielem)/L^2;...
                0,               6*E(ielem)*I(ielem)/L^2,     2*E(ielem)*I(ielem)/L];

            kjj=[E(ielem)*A(ielem)/L,               0,                             0;...
                0,             12*E(ielem)*I(ielem)/L^3,     -6*E(ielem)*I(ielem)/L^2;...
                0,             -6*E(ielem)*I(ielem)/L^2,      4*E(ielem)*I(ielem)/L];


            kii=T'*kii*T;kij=T'*kij*T;kji=T'*kji*T;kjj=T'*kjj*T;
            %compute local internal forces
            fi=kii*ui+kij*uj;fj=kji*ui+kjj*uj;
            fi=T*fi;fj=T*fj;
            %extract axial, shear and bending moment for Node 1
            Finternal(ielem,1,1)=fi(1,1);
            Finternal(ielem,1,2)= fi(2,1);
            Finternal(ielem,1,3)=fi(3,1);
            %extract axial, shear and bending moment for Node 2
            Finternal(ielem,2,1)= fj(1,1);
            Finternal(ielem,2,2)= fj(2,1);
            Finternal(ielem,2,3)= fj(3,1);
        end

        Finternal = Finternal + UDLFinternal;

        Finternal(:,2,1) = - Finternal(:,2,1);
        Finternal(:,2,3) = - Finternal(:,2,3);

    end

% Function to plot Bending Moment and Shear Force with scaling

    function plot_BM_and_SF(X, connectivity, Finternal, w, scale)

        figure(1); clf; hold on;
        title('Bending Moment Diagram (Global)');
        xlabel('X (m)'); ylabel('Y (m)');
        axis equal; grid on;

        figure(2); clf; hold on;
        title('Shear Force Diagram (Global)');
        xlabel('X (m)'); ylabel('Y (m)');
        axis equal; grid on;

        for e = 1:size(connectivity, 1)
            % === node ===
            i = connectivity(e,1);
            j = connectivity(e,2);
            xi = X(i,1); yi = X(i,2);
            xj = X(j,1); yj = X(j,2);

            L = norm([xj - xi, yj - yi]);
            dx = (xj - xi)/L; dy = (yj - yi)/L;
            R1 = [dx, -dy; dy, dx];  % local → global

            % === Internal force ===
            MA = Finternal(e,1,3);  VA = Finternal(e,1,2);
            MB = Finternal(e,2,3);  VB = Finternal(e,2,2);

            x_local = linspace(0, L, 50);

            % === BM ===
            M_local = -MA + VA * x_local + 0.5 * w(e) * x_local.^2;
            yM_local = -M_local * scale;

            coords_M = R1 * [x_local; yM_local] + [xi; yi];

            figure(1);
            plot([xi, xj], [yi, yj], 'k-', 'LineWidth', 1.5);
            plot(coords_M(1,:), coords_M(2,:), 'r-', 'LineWidth', 1.5);


            epsilon = 1e-4 * (-1)^e;

            for k = 1:5:length(x_local)
                xk = x_local(k);
                yk = yM_local(k);
                p0 = R1 * [xk; 0] + [xi; yi];
                p1 = R1 * [xk; yk] + [xi; yi];
                plot([p0(1), p1(1)], [p0(2), p1(2)], 'r--');
            end


            p0_tail = R1 * [L; 0] + [xi; yi];
            p1_tail = R1 * [L; yM_local(end) + epsilon] + [xi; yi];
            plot([p0_tail(1), p1_tail(1)], [p0_tail(2), p1_tail(2)], 'r--');

            % === SF  ===
            V_local = VA + w(e) * x_local;
            yV_local = V_local * scale;
            coords_V = R1 * [x_local; yV_local] + [xi; yi];

            figure(2);
            plot([xi, xj], [yi, yj], 'k-', 'LineWidth', 1.5);
            plot(coords_V(1,:), coords_V(2,:), 'b-', 'LineWidth', 1.5);

            for k = 1:5:length(x_local)
                xk = x_local(k);
                yk = yV_local(k);
                p0 = R1 * [xk; 0] + [xi; yi];
                p1 = R1 * [xk; yk] + [xi; yi];
                plot([p0(1), p1(1)], [p0(2), p1(2)], 'b--');
            end


            p0_tail_V = R1 * [L; 0] + [xi; yi];
            p1_tail_V = R1 * [L; yV_local(end) + epsilon] + [xi; yi];
            plot([p0_tail_V(1), p1_tail_V(1)], [p0_tail_V(2), p1_tail_V(2)], 'b--');
        end

        figure(1); legend('Structure', 'BM Curve');
        figure(2); legend('Structure', 'SF Curve');
    end

    function BM_matrix = plot_BM_only(X, connectivity, Finternal, w, scale, draw)

        num_elem = size(connectivity, 1);
        BM_matrix = zeros(11, num_elem);

        if draw
            figure(1); clf; hold on;
            title('Bending Moment Diagram');
            xlabel('X (m)'); ylabel('Y (m)');
            axis equal; grid on;
        end

        %
        annotated_values = containers.Map('KeyType','char','ValueType','any');

        for e = 1:num_elem
            n1 = connectivity(e,1); n2 = connectivity(e,2);
            xn1 = X(n1,1); yn1 = X(n1,2);
            xn2 = X(n2,1); yn2 = X(n2,2);
            L_elem = norm([xn2 - xn1, yn2 - yn1]);

            cx = (xn2 - xn1)/L_elem; cy = (yn2 - yn1)/L_elem;
            R_bm = [cx, -cy; cy, cx];

            MA = Finternal(e,1,3);
            VA = Finternal(e,1,2);

            x_bm_local = linspace(0, L_elem, 11);
            M_bm = -MA + VA * x_bm_local + 0.5 * w(e) * x_bm_local.^2;
            offset_bm = -M_bm * scale;

            BM_matrix(:, e) = M_bm(:);
            coords_bm = R_bm * [x_bm_local; offset_bm] + [xn1; yn1];

            if draw
                plot([xn1, xn2], [yn1, yn2], 'k-', 'LineWidth', 1.5);
                plot(coords_bm(1,:), coords_bm(2,:), 'r-', 'LineWidth', 1.5);

                epsilon = 1e-4 * (-1)^e;
                for k = 1:length(x_bm_local)
                    xk = x_bm_local(k); yk = offset_bm(k);
                    p0 = R_bm * [xk; 0] + [xn1; yn1];
                    p1 = R_bm * [xk; yk] + [xn1; yn1];
                    plot([p0(1), p1(1)], [p0(2), p1(2)], 'r--');
                end
                p0_tail = R_bm * [L_elem; 0] + [xn1; yn1];
                p1_tail = R_bm * [L_elem; offset_bm(end) + epsilon] + [xn1; yn1];
                plot([p0_tail(1), p1_tail(1)], [p0_tail(2), p1_tail(2)], 'r--');

                % 选中 3 个位置：起点、中点、终点
                idxs = [1, 6, 11];
                node_ids = [n1, -1, n2];  % 中点没有 node id
                pts_local = [x_bm_local(idxs); offset_bm(idxs)];
                pts_global = R_bm * pts_local + [xn1; yn1];

                for j = 1:3
                    BM_val = M_bm(idxs(j));
                    if abs(BM_val) < 1e-6
                        continue;  %
                    end
                    pos = pts_global(:,j);
                    if j == 2
                        % 中点 — 始终绘制圆圈和文本
                        plot(pos(1), pos(2), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
                        text(pos(1), pos(2), sprintf('%.2f kN·m', BM_val), ...
                            'FontSize', 12, 'Color', 'r', 'VerticalAlignment','top');
                    else
                        node_id_bm = node_ids(j);
                        key = sprintf('N%d', node_id_bm);
                        if ~isKey(annotated_values, key)
                            annotated_values(key) = [];
                        end
                        if ~ismember(round(BM_val,6), round(annotated_values(key),6))
                            text(pos(1), pos(2), sprintf('%.2f kN·m', BM_val), ...
                                'FontSize', 12, 'Color', 'r', 'VerticalAlignment','bottom');
                            annotated_values(key) = [annotated_values(key), BM_val];
                        end
                    end
                end
            end
        end

        if draw
            xlimits = xlim;
            ylimits = ylim;
            xpad = 0.02 * (xlimits(2) - xlimits(1));
            ypad = 0.02 * (ylimits(2) - ylimits(1));
            xlim([xlimits(1) - xpad, xlimits(2) + xpad]);
            ylim([ylimits(1) - ypad, ylimits(2) + ypad]);
            legend('Structure', 'BM Curve');
        end
    end

    function SF_matrix = plot_SF_only(X, connectivity, Finternal, w, scale, draw)

        num_elem = size(connectivity, 1);
        SF_matrix = zeros(11, num_elem);

        if draw
            figure(2); clf; hold on;
            title('Shear Force Diagram');
            xlabel('X (m)'); ylabel('Y (m)');
            axis equal; grid on;
        end

        % Dictionary to track annotated values at nodes
        annotated_values = containers.Map('KeyType','char','ValueType','any');

        for e = 1:num_elem
            n1 = connectivity(e,1); n2 = connectivity(e,2);
            xn1 = X(n1,1); yn1 = X(n1,2);
            xn2 = X(n2,1); yn2 = X(n2,2);
            L_elem = norm([xn2 - xn1, yn2 - yn1]);

            cx = (xn2 - xn1)/L_elem;
            cy = (yn2 - yn1)/L_elem;
            R_sf = [cx, -cy; cy, cx];  % Rotation matrix

            VA = Finternal(e,1,2);  % Shear force at node A

            x_sf_local = linspace(0, L_elem, 11);
            V_sf = VA + w(e) * x_sf_local;  % Shear force at 11 points
            offset_sf = V_sf * scale;

            SF_matrix(:, e) = V_sf(:);
            coords_sf = R_sf * [x_sf_local; offset_sf] + [xn1; yn1];

            if draw
                % Draw element line and SF curve
                plot([xn1, xn2], [yn1, yn2], 'k-', 'LineWidth', 1.5);
                plot(coords_sf(1,:), coords_sf(2,:), 'b-', 'LineWidth', 1.5);

                epsilon = 1e-4 * (-1)^e;
                for k = 1:length(x_sf_local)
                    xk = x_sf_local(k); yk = offset_sf(k);
                    p0 = R_sf * [xk; 0] + [xn1; yn1];
                    p1 = R_sf * [xk; yk] + [xn1; yn1];
                    plot([p0(1), p1(1)], [p0(2), p1(2)], 'b--');  % Dotted shear lines
                end

                % Tail-end vertical line
                p0_tail = R_sf * [L_elem; 0] + [xn1; yn1];
                p1_tail = R_sf * [L_elem; offset_sf(end) + epsilon] + [xn1; yn1];
                plot([p0_tail(1), p1_tail(1)], [p0_tail(2), p1_tail(2)], 'b--');

                % Annotate values at start, mid, and end of each element
                idxs = [1, 6, 11];  % start, middle, end
                node_ids = [n1, -1, n2];  % node IDs for start and end, -1 for mid
                pts_local = [x_sf_local(idxs); offset_sf(idxs)];
                pts_global = R_sf * pts_local + [xn1; yn1];

                for j = 1:3
                    SF_val = V_sf(idxs(j));
                    if abs(SF_val) < 1e-6
                        continue;  % Skip near-zero values
                    end
                    pos = pts_global(:,j);
                    if j == 2
                        % Middle point: draw circle and label
                        plot(pos(1), pos(2), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
                        text(pos(1), pos(2), sprintf('%.2f kN', SF_val), ...
                            'FontSize', 12, 'Color', 'b', 'VerticalAlignment','top');
                    else
                        node_id_sf = node_ids(j);
                        key = sprintf('N%d', node_id_sf);
                        if ~isKey(annotated_values, key)
                            annotated_values(key) = [];
                        end
                        if ~ismember(round(SF_val,6), round(annotated_values(key),6))
                            text(pos(1), pos(2), sprintf('%.2f kN', SF_val), ...
                                'FontSize', 12, 'Color', 'b', 'VerticalAlignment','bottom');
                            annotated_values(key) = [annotated_values(key), SF_val];
                        end
                    end
                end
            end
        end

        if draw
            % Auto-expand axes boundaries
            xlimits = xlim;
            ylimits = ylim;
            xpad = 0.02 * (xlimits(2) - xlimits(1));
            ypad = 0.02 * (ylimits(2) - ylimits(1));
            xlim([xlimits(1) - xpad, xlimits(2) + xpad]);
            ylim([ylimits(1) - ypad, ylimits(2) + ypad]);
            legend('Structure', 'SF Curve');
        end
    end

    function [TotalYd, TotalXd, TotalGlobal] = plot_deformation_local(X, connectivity, U, Ud, Nelem, w, E, I, A, Finternal, fixdof, UFINAL)

        num_elem = size(connectivity,1);
        TotalYd = zeros(11, num_elem);        % Local y-direction (deflection)
        TotalXd = zeros(11, num_elem);        % Local x-direction (axial displacement)
        TotalGlobal = zeros(11, 2, num_elem); % Global [ΔX, ΔY] displacement at 11 points per element

        for ielem = 1:num_elem

            % === Node coordinates and direction cosines ===
            i = connectivity(ielem,1);
            j = connectivity(ielem,2);
            xi = X(i,1); yi = X(i,2);
            xj = X(j,1); yj = X(j,2);
            L = norm([xj - xi, yj - yi]);
            c = (xj - xi)/L;
            s = (yj - yi)/L;

            % === Degrees of freedom ===
            idof = Nelem(:,1:3,ielem);
            jdof = Nelem(:,4:6,ielem);

            % === Transform global displacements to local (axial) ===
            U1_global = [U(idof(1)); U(idof(2))];
            U2_global = [U(jdof(1)); U(jdof(2))];
            R2 = [c s; -s c];  % Rotation matrix
            U1_local = R2 * U1_global;
            U2_local = R2 * U2_global;
            u1 = U1_local(1);
            u2 = U2_local(1);

            % === Rotate global displacements to obtain local deflection ===
            T = [c s 0; -s c 0; 0 0 1];
            T2 = [T, T; T, T];
            Udx = U;
            x_start = Udx(idof(1));
            x_end = Udx(jdof(1));

            Ud(idof) = T2(1:3,1:3) * U(idof);
            Ud(jdof) = T2(4:6,4:6) * U(jdof);

            y1 = Ud(idof(2));
            r1 = Ud(idof(3));
            y2 = Ud(jdof(2));
            r2 = Ud(jdof(3));

            % === Beam shape function parameter computation ===
            ya = y1; ra = r1;
            yb = y2; rb = r2;
            d = ya;
            c1 = ra;

            MA = [L^3/6, L^2/2; L^2/2, L];
            Mb = [yb - ya - ra*L - w(ielem)*L^4/(24*E(ielem)*I(ielem));
                rb - ra - w(ielem)*L^3/(6*E(ielem)*I(ielem))];
            Rab = MA \ Mb;
            a = Rab(1);
            b = Rab(2);

            % === Displacement interpolation along the element ===
            Yd = zeros(11,1);
            Xd = zeros(11,1);

            for k = 1:11
                xk = (k-1)/10 * L;

                % Local deflection using beam theory
                Yd(k) = w(ielem)*xk^4/(24*E(ielem)*I(ielem)) + a*xk^3/6 + b*xk^2/2 + c1*xk + d;

                % Local axial displacement using linear interpolation
                Xd(k) = (u2 - u1)/10*(k-1);

                % === Transform local displacement to global coordinates ===
                local_disp = [Xd(k); Yd(k)];
                global_disp = [c, -s; s, c] * local_disp;
                TotalGlobal(k,:,ielem) = global_disp';  % Store [ΔX, ΔY]
            end

            % === Overwrite first and last points for higher accuracy ===
            TotalGlobal(1,1,ielem)  = x_start;
            TotalGlobal(end,1,ielem) = x_end;
            Xd(1)  = u1;
            Xd(end) = u2;
            Yd(1)  = ya;
            Yd(end) = yb;

            % === Store results ===
            TotalYd(:,ielem) = Yd;
            TotalXd(:,ielem) = Xd;
        end
    end

    function plot_deformation_global(X, connectivity, UFINAL, Ndof, TotalGlobal, scale)
        figure; hold on; axis equal; grid on;
        title('Deformed Structure (Global Coordinates)');
        xlabel('X (m)'); ylabel('Y (m)');

        num_elem = size(connectivity, 1);
        nnode = size(X, 1);

        % === Compute the deformed coordinates of each node (used only to correct endpoints) ===
        X_def = X;
        for n = 1:nnode
            X_def(n,1) = X(n,1) + scale * UFINAL(Ndof(n,1));  % x
            X_def(n,2) = X(n,2) + scale * UFINAL(Ndof(n,2));  % y
        end

        % === Plot the original structure using undeformed coordinates ===
        for e = 1:num_elem
            n1 = connectivity(e,1);
            n2 = connectivity(e,2);
            x1 = X(n1,1); y1 = X(n1,2);
            x2 = X(n2,1); y2 = X(n2,2);
            if e == 1
                plot([x1, x2], [y1, y2], 'k-', 'LineWidth', 1.5, 'DisplayName', 'Original Structure');
            else
                plot([x1, x2], [y1, y2], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
        end

        % === Plot the actual deformed curves using TotalGlobal ===
        for e = 1:num_elem
            n1 = connectivity(e,1); n2 = connectivity(e,2);
            x1 = X(n1,1); y1 = X(n1,2);
            x2 = X(n2,1); y2 = X(n2,2);

            x_orig = linspace(x1, x2, 11);
            y_orig = linspace(y1, y2, 11);

            dx = TotalGlobal(:,1,e) * scale;
            dy = TotalGlobal(:,2,e) * scale;

            x_def = x_orig(:) + dx;
            y_def = y_orig(:) + dy;

            % === Replace endpoints with true deformed coordinates (prevent node tearing) ===
            x_def(1)  = X_def(n1,1);
            y_def(1)  = X_def(n1,2);
            x_def(11) = X_def(n2,1);
            y_def(11) = X_def(n2,2);

            if e == 1
                plot(x_def, y_def, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Deformed Shape');
            else
                plot(x_def, y_def, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
        end

        % === Automatically adjust axis limits with padding ===
        xlimits = xlim; ylimits = ylim;
        xpad = 0.02 * (xlimits(2) - xlimits(1));
        ypad = 0.02 * (ylimits(2) - ylimits(1));
        xlim([xlimits(1) - xpad, xlimits(2) + xpad]);
        ylim([ylimits(1) - ypad, ylimits(2) + ypad]);

        legend show
    end



    function Finternal_local = global2local_Finternal(Finternal, connectivity, X)
        % Convert global Finternal (Axial, Shear) into local coordinates for each element
        Finternal_local = Finternal;

        for e = 1:size(connectivity, 1)
            n1 = connectivity(e, 1);
            n2 = connectivity(e, 2);
            X1 = X(n1, :); X2 = X(n2, :);

            L = norm(X2 - X1);
            cx = (X2(1) - X1(1)) / L;
            cy = (X2(2) - X1(2)) / L;

            % Rotation matrix: from global to local
            Rl = [cx, cy; -cy, cx];

            for node = 1:2
                % Extract global Axial + Shear
                Fg = squeeze(Finternal(e, node, 1:2));  % (2x1)
                Fl = Rl * Fg;                            % (2x1)
                Finternal_local(e, node, 1) = Fl(1);    % Local axial
                Finternal_local(e, node, 2) = Fl(2);    % Local shear
                % Moment (index 3) is assumed unchanged (in-plane)
            end
        end
    end

end