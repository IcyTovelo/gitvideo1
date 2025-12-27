function  [ndof,X,connectivity,Ndof,Nelem,totaldof,E,A,I,H,F,U,alldof,Elastic,Bnode,w,forceunits,lengthunits,UDLFinternal]=Read_input(datafilename)

%Structural Mechanics
%Input data Function of Matrix analysis of frame
%Program created by Dr.Oubay Hassan and Shen Ma

% Open the file
fid = fopen(datafilename, 'r');

% Check if the file was successfully opened
if fid == -1
    error('Cannot open file: %s', datafilename);
end

% Read and display the title
AnalysisTitle = fscanf(fid, 'TITLE = %[^\n]', 1);


% Move the file pointer to the next line
fgetl(fid);

% Read the subtitle

AnalysisStuff = fscanf(fid, 'SUBTITLE = %[^\n]', 1);

% Read the input list
stufflist = fscanf(fid, '\n%d %d %d %d %d %d %d', [7,1]);


% =========================================================================
% Data input phase. Read information from input data file 'datafilename.txt'
% =========================================================================


ndof  = 3;
% Unit section
lengthunits = fscanf(fid, '\nLENGTHUNITS = %s', 1);
forceunits  = fscanf(fid, '\nAXIALUNITS = %s', 1);



%-------------------------------------------------------
% Table of node coordinates
%
% 1.X  2.HN  3.HCT  4.Ndof
%-------------------------------------------------------

npoin = fscanf(fid, '\nNODE_COORDINATES = %d', 1);

if npoin ~= stufflist(1,1)
    disp('error: Number of node is wrong')
end

coord = fscanf(fid, '\n%d %f %f', [3,npoin]);
coord = coord';
coordinate = coord(:,1);

coordRows = size(coord,1);
%X: node coordinate
X = coord(:,2:3);





%-------------------------------------------------------
% Table of connectivities
%
% 1.emp  2.A  3.I  4.connectivity
%-------------------------------------------------------


% Total number of elements in the mesh
nelem = fscanf(fid, '\nELEMENTS = %d', 1);

if nelem ~= stufflist(2,1)

    disp('error: Number of elements is wrong')
end

lnods = fscanf(fid, '\n%d %d %d %d %f %f %d %d', [8,nelem]);
lnods = lnods';

%element metrials type
emp = lnods(:,4);
%cross-sectional areas
A = lnods(:,5);
%cross-sectional areas
I = lnods(:,6);

H = lnods(:,7:8);
%store only nodal connectivities in lnods
connectivity = lnods(:,2:3);
% Initialize an array to store the DOF for each node (3 DOFs per node)
nodeDOF = zeros(size(X,1), 3);
Ndof    = zeros(size(X,1), 3);
% Initialize next available DOF and a counter for the number of DOFs assigned
nextDOF = 1;
Nelem = zeros(1,6,size(nelem,1));
ncnode = zeros(size(X,1),3);



% Loop through all elements
for i = 1:nelem
    node1 = lnods(i, 2);  % First node of the element
    node2 = lnods(i, 3);  % Second node of the element
    modifyNode1 = lnods(i, 7);  % Modification flag for node1 (second-to-last column)
    modifyNode2 = lnods(i, 8);

    % Check and assign DOFs for node1


    if nodeDOF(node1, 1) == 0  % Node1 has not been assigned DOFs
        % Assign the first DOF to node1
        nodeDOF(node1, 1) = nextDOF;  % Assign current nextDOF
        nextDOF = nextDOF + 1;        % Increment nextDOF after assigning the first DOF

        % Assign the second DOF to node1
        nodeDOF(node1, 2) = nextDOF;  % Assign updated nextDOF
        nextDOF = nextDOF + 1;        % Increment nextDOF after assigning the second DOF

        % Assign the third DOF to node1
        nodeDOF(node1, 3) = nextDOF;  % Assign updated nextDOF
        nextDOF = nextDOF + 1;        % Increment nextDOF after assigning the third DOF

        Ndof(node1, 1:3) = nodeDOF(node1, 1:3);
        Nelem(1,1:3,i) = nodeDOF(node1, :);



        if modifyNode1 == 1  % Modify node1's third DOF

            ncnode(node1,1:2) =  nodeDOF(node1,1:2);
            ncnode(node1,3)   = -nodeDOF(node1,3);

            Nelem(1,1:3,i)  = abs(ncnode(node1,:));
            Ndof(node1, :)  = ncnode(node1,:);
        else
            ncnode(node1,1:3) =  nodeDOF(node1,1:3);
            Nelem(1,1:3,i)  = nodeDOF(node1, :);
            Ndof(node1, :)  = ncnode(node1,:);
        end

    elseif nodeDOF(node1, 1) ~= 0

        if modifyNode1 == 0 && Ndof(node1, 3) <= 0 % Modify node1's third DOF

            nodeDOF(node1, 3) = nextDOF;  % Assign updated nextDOF
            nextDOF = nextDOF + 1;
            ncnode(node1,:) =  nodeDOF(node1,:);
            Ndof(node1, 3) = ncnode(node1, 3);
            Nelem(1,1:3,i)  = ncnode(node1,:);
        elseif modifyNode1 == 0 && Ndof(node1, 3) >= 0
            ncnode(node1,:) =  nodeDOF(node1,:);
            Nelem(1,1:3,i)  = Ndof(node1, :);
        elseif modifyNode1 == 1
            nodeDOF(node1, 3) = nextDOF;  % Assign updated nextDOF
            ncnode(node1,:) =  nodeDOF(node1,:);
            nextDOF = nextDOF + 1;
            Nelem(1,1:3,i)  = ncnode(node1,:);
        end


    end

    if nodeDOF(node2, 1) == 0  % Node1 has not been assigned DOFs
        % Assign the first DOF to node1
        nodeDOF(node2, 1) = nextDOF;  % Assign current nextDOF
        nextDOF = nextDOF + 1;        % Increment nextDOF after assigning the first DOF

        % Assign the second DOF to node1
        nodeDOF(node2, 2) = nextDOF;  % Assign updated nextDOF
        nextDOF = nextDOF + 1;        % Increment nextDOF after assigning the second DOF

        % Assign the third DOF to node1
        nodeDOF(node2, 3) = nextDOF;  % Assign updated nextDOF
        nextDOF = nextDOF + 1;
        Ndof(node2, 1:3) = nodeDOF(node2, 1:3);
        Nelem(1,4:6,i) = nodeDOF(node2, :);

        if modifyNode2 == 1  % Modify node1's third DOF

            ncnode(node2,1:2) =  nodeDOF(node2,1:2);
            ncnode(node2,3)   = -nodeDOF(node2,3);

            Nelem(1,4:6,i)  = abs(ncnode(node2,:));
            Ndof(node2, :)  = ncnode(node2,:);
        else
            ncnode(node2,1:3) =  nodeDOF(node2,1:3);
            Nelem(1,4:6,i)  = nodeDOF(node2, :);
            Ndof(node2, :)  = ncnode(node2,:);
        end


    elseif nodeDOF(node2, 1) ~= 0


        if modifyNode2 == 0 && Ndof(node2, 3) <= 0 % Modify node1's third DOF

            nodeDOF(node2, 3) = nextDOF;  % Assign updated nextDOF
            nextDOF = nextDOF + 1;
            ncnode(node2,:) =  nodeDOF(node2,:);
            Ndof(node2, 3) = ncnode(node2, 3);
            Nelem(1,4:6,i)  = ncnode(node2,:);
        elseif modifyNode2 == 0 && Ndof(node2, 3) >= 0
            ncnode(node2,:) =  nodeDOF(node2,:);

            Nelem(1,4:6,i)  = Ndof(node2, :);
        elseif modifyNode2 == 1
            nodeDOF(node2, 3) = nextDOF;  % Assign updated nextDOF
            ncnode(node2,:) =  nodeDOF(node2,:);
            nextDOF = nextDOF + 1;
            Nelem(1,4:6,i)  = ncnode(node2,:);
        end

        % else
        %       ncnode(node2,1:3) =  nodeDOF(node2,1:3);
        %
        %       Nelem(1,4:6,i)  = nodeDOF(node2, :);
        %       Ndof(node2, :)  = ncnode(node2,:);
    end

end

% disp(Nelem)
% disp(Ndof)



% disp(Nelem)
totaldof = (nextDOF-1);


alldof  = zeros(3*size(Ndof,1),1);
U       = zeros(3*size(Ndof,1),1);
Elastic = zeros(3*size(Ndof,1),1);
% disp(size(alldof))



%-------------------------------------------------------
% Table of material properties
%
% 1.E  2.alpha  3.rho
%-------------------------------------------------------


totalemp = fscanf(fid, '\nMATERIAL PROPERTIES = %d', 1);

if totalemp ~= stufflist(3,1)
    disp('error: Number of material properties is wrong')
end

%empdata是MATERIAL PROPERTIES的data矩阵
empdata = fscanf(fid, '\n%d %e %e %f',[4,totalemp]);

empdata = empdata';

%element material properties

%read all type of material properties
Elist     = empdata(:,2);
alphalist = empdata(:,3);
rholist   = empdata(:,4);

%initial allocation
E     = zeros(size(emp));
alpha = zeros(size(emp));
rho   = zeros(size(emp));
% loop over elements
for i = 1:length(emp)
    if emp(i) == 1
        E(i) = Elist(1);
        alpha(i) = alphalist(1);
        rho(i) = rholist(1);
    elseif emp(i) == 2
        E(i) = Elist(2);
        alpha(i) = alphalist(2);
        rho(i) = rholist(2);
    elseif emp(i) == 3
        E(i) = Elist(3);
        alpha(i) = alphalist(3);
        rho(i) = rholist(3);
    end
end








%-------------------------------------------------------
% Table of boundary conditions
%
% 1.fix boundary conditions 2.prescribed displacements
%-------------------------------------------------------

nnodefix = fscanf(fid,'\nNODES_WITH_PRESCRIBED_DISPLACEMENTS = %d',1);
if nnodefix ~= stufflist(4,1)
    disp('error: Number of BC is wrong')
end
ipresc   = fscanf(fid, '\n%d %d %d %d %f %f %f %f %f %f', [10,nnodefix]);
ipresc   = ipresc';
Bnode    = ipresc(:,1);


for ibc=1:size(ipresc,1)

    bnode = Ndof(Bnode(ibc),:);
    alldof(bnode)  = ipresc(ibc, 2:4);
    U(bnode)       = ipresc(ibc, 5:7);
    Elastic(bnode) = ipresc(ibc, 8:10);

end





%-------------------------------------------------------
% Table of load
%
% 1.P1  2.P2  3.UDLFinternal 4.F
%-------------------------------------------------------

%----  ----  ----  ----
% 1. Table of point load
%----  ----  ----  ----
nnodP = fscanf(fid,'\nNODES_WITH_POINT_LOAD = %d',1);

%initial allocation
F  = zeros(totaldof,1);
P1 = zeros(totaldof,1);
P2 = zeros(totaldof,1);

Pndof = zeros(nnodP,3);
Pnode = zeros(nnodP,1);

UDLFinternal = zeros(size(connectivity, 1), 2, 3);


if nnodP ~= stufflist(5,1)
    disp('error: Number of point load is wrong')
end
if nnodP == 0
    disp('No point load')
else
    Pdata1 = fscanf(fid,'\n%d %f %f %f', [4, nnodP]);
    Pdata1 = Pdata1';
    dnode = Pdata1(:,2);

    for iP = 1:nnodP


        Pnode(iP) = Pdata1(iP,1);

        P1(Ndof(Pnode(iP),1:2)) = Pdata1(iP,2:3);
        if Ndof(Pnode(iP),3) >= 0
            P1(Ndof(Pnode(iP),3)) = Pdata1(iP,4);
        end

    end
    F=P1;
end


%---- Read UDL Load Table ----
nelep = fscanf(fid,'\nELEMENT_WITH_UDL_LOAD = %d',1);
if nelep ~= stufflist(6,1)
    disp('Error: Number of UDL loads does not match stufflist.');
end

% Initialize w vector: local y-direction UDL for each element
w = zeros(size(connectivity,1),1);

if nelep ~= 0
    Pdata2 = fscanf(fid,'\n%d %d %d %f', [4, nelep]);
    Pdata2 = Pdata2';  % Transpose for row-wise access

    elegather = Pdata2(:,1);      % Element numbers
    udl_type  = Pdata2(:,2);      % Type (1–4)
    udl_dir   = Pdata2(:,3);      % Direction flag
    udl_value = Pdata2(:,4);      % UDL magnitude

    for i = 1:nelep
        e   = elegather(i);
        typ = udl_type(i);
        dir = udl_dir(i);
        q   = udl_value(i);

        % Determine the sign for local y-direction
        sign = 1;

        % For types where downward UDL → negative in local y

        
        if (typ == 2 && dir == 2) || ...       % Global Y, downward
           (typ == 3 && dir == 2) || ...       % Diagonal x+/y−
           (typ == 1 && dir == 1) || ...       
           (typ == 4 && dir == 2)              % Diagonal x−/y−
            sign = -1;
        end

        w(e) = w(e) + sign * q;
    end
end



% === Process UDL into nodal forces ===
for i = 1:nelep
    elem = elegather(i);
    type = udl_type(i);
    flag = udl_dir(i);
    val  = udl_value(i);

    % Get element geometry
    nodei = connectivity(elem,1);
    nodej = connectivity(elem,2);
    idof = Nelem(1,1:3,elem);
    jdof = Nelem(1,4:6,elem);

    Xi = X(nodei,:); Xj = X(nodej,:);
    Xij = Xj - Xi;
    L = norm(Xij);
    c = Xij(1)/L;
    s = Xij(2)/L;

    % UDL converted to nodal equivalent force and moment
    Ppoint = val * L / 2;
    Pbending = val * L^2 / 12;

    % === CASE 1: Global X+, type=1, flag=1 ===
    if type == 1 && flag == 1
        Px = Ppoint;  % Only in x
        Py = 0;
        Pz = Pbending;

        Ptotal1 = [Px; Py; -Pz];
        Ptotal2 = [Px; Py;  Pz];

        P2(idof,1) = P2(idof,1) + Ptotal1;
        P2(jdof,1) = P2(jdof,1) + Ptotal2;
        F = P1 + P2;

        % === CASE 2: Global X-, type=1, flag=2 ===
    elseif type == 1 && flag == 2
        Px = -Ppoint;  % Only in -x direction
        Py = 0;
        Pz = Pbending;

        Ptotal1 = [Px; Py; -Pz];
        Ptotal2 = [Px; Py;  Pz];

        P2(idof,1) = P2(idof,1) + Ptotal1;
        P2(jdof,1) = P2(jdof,1) + Ptotal2;
        F = P1 + P2;

        % === CASE 3: Global Y+, type=2, flag=1 ===
    elseif type == 2 && flag == 1
        Px = 0;
        Py = Ppoint;
        Pz = Pbending;

        Ptotal1 = [Px; Py; -Pz];
        Ptotal2 = [Px; Py;  Pz];

        P2(idof,1) = P2(idof,1) + Ptotal1;
        P2(jdof,1) = P2(jdof,1) + Ptotal2;
        F = P1 + P2;
        % === CASE 4: Global Y−, type=2, flag=2 ===
    elseif type == 2 && flag == 2
        Px = 0;
        Py = -Ppoint;
        Pz = Pbending;

        Ptotal1 = [Px; Py; -Pz];
        Ptotal2 = [Px; Py;  Pz];

        P2(idof,1) = P2(idof,1) + Ptotal1;
        P2(jdof,1) = P2(jdof,1) + Ptotal2;
        F = P1 + P2;
        % === CASE 5: UDL perpendicular to element (left side), type=3, flag=1 ===
    elseif type == 3 && flag == 1
        nx = -s;  % s = sin(theta)
        ny =  c;  % c = cos(theta)

        Px = nx * Ppoint;
        Py = ny * Ppoint;
        Pz = Pbending;

        Ptotal1 = [Px; Py; -Pz];
        Ptotal2 = [Px; Py;  Pz];

        P2(idof,1) = P2(idof,1) + Ptotal1;
        P2(jdof,1) = P2(jdof,1) + Ptotal2;
        F = P1 + P2;

        % === CASE 6: UDL perpendicular to element (right side), type=3, flag=2 ===
    elseif type == 3 && flag == 2
        nx =  s;  % s = sin(theta)
        ny = -c;  % c = cos(theta)

        Px = nx * Ppoint;
        Py = ny * Ppoint;
        Pz = Pbending;

        Ptotal1 = [Px; Py; -Pz];
        Ptotal2 = [Px; Py;  Pz];

        P2(idof,1) = P2(idof,1) + Ptotal1;
        P2(jdof,1) = P2(jdof,1) + Ptotal2;
        F = P1 + P2;


        % === CASE 7: UDL perpendicular to element (left side), type=4, flag=1 ===
    elseif type == 4 && flag == 1
        nx = -s;
        ny =  c;

        Px = nx * Ppoint;
        Py = ny * Ppoint;
        Pz = Pbending;

        Ptotal1 = [Px; Py; -Pz];
        Ptotal2 = [Px; Py;  Pz];

        P2(idof,1) = P2(idof,1) + Ptotal1;
        P2(jdof,1) = P2(jdof,1) + Ptotal2;
        F = P1 + P2;



        % === CASE 8: UDL perpendicular to element (right side), type=4, flag=2 ===
    elseif type == 4 && flag == 2
        nx =  s;
        ny = -c;

        Px = nx * Ppoint;
        Py = ny * Ppoint;
        Pz = Pbending;

        Ptotal1 = [Px; Py; -Pz];
        Ptotal2 = [Px; Py;  Pz];

        P2(idof,1) = P2(idof,1) + Ptotal1;
        P2(jdof,1) = P2(jdof,1) + Ptotal2;
        F = P1 + P2;










    end
end

%----  ----  ----  ----  ----
% 3. Add UDL Finternal effects (Full 8-type support with local coordinate correction)
%----  ----  ----  ----  ----

for ie = 1:size(elegather, 1)
    elem_id = elegather(ie);
    type    = udl_type(ie);      % Type (1–4)
    flag    = udl_dir(ie);       % Direction flag (1 or 2)
    q       = udl_value(ie);     % UDL magnitude (kN/m)

    % === Element geometry ===
    ni = connectivity(elem_id, 1);   % Start node
    nj = connectivity(elem_id, 2);   % End node
    Xi = X(ni, :); Xj = X(nj, :);    % Coordinates

    L = norm(Xj - Xi);               % Element length
    cx = (Xj(1) - Xi(1)) / L;        % Cosine X
    cy = (Xj(2) - Xi(2)) / L;        % Cosine Y

    % === Equivalent nodal force/moment in global ===
    P = q * L / 2;         % Equivalent nodal force
    M = q * L^2 / 12;      % Equivalent nodal moment

    % === Direction-based global Px / Py ===
    switch type
        case 1  % Along global X axis
            Px = (flag == 1)*P + (flag == 2)*(-P);
            Py = 0;

        case 2  % Along global Y axis
            Px = 0;
            Py = (flag == 1)*P + (flag == 2)*(-P);
            M  = (flag == 1)*(-M) + (flag == 2)*M;

        case 3  % Diagonal (x+, y±)
            Px = P / sqrt(2);
            Py = Px * ((flag == 1)*1 + (flag == 2)*-1);

        case 4  % Diagonal (x−, y±)
            Px = -P / sqrt(2);
            Py = abs(Px) * ((flag == 1)*1 + (flag == 2)*-1);

        otherwise
            Px = 0;
            Py = 0;
    end

    % === Global force → Local (axial / shear) ===
    T_force = [cx, cy; -cy, cx];      % Rotation matrix
    F_local = T_force * [Px; Py];     % [Axial; Shear]

    Plocal_x = F_local(1);  % Axial force in local x
    Plocal_y = F_local(2);  % Shear force in local y

    % === Final moment signs ===
    M_local = -M;  % M: left node -M, right node +M (counterclockwise convention)

    % === Fill UDLFinternal matrix ===
    UDLFinternal(elem_id, 1, 1) = -Plocal_x;
    UDLFinternal(elem_id, 1, 2) = -Plocal_y;
    UDLFinternal(elem_id, 1, 3) = M_local;

    UDLFinternal(elem_id, 2, 1) = -Plocal_x;
    UDLFinternal(elem_id, 2, 2) = -Plocal_y;
    UDLFinternal(elem_id, 2, 3) = -M_local;
end


end
