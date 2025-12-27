[M] = compute_mass_matrix(X, connectivity, Nelem, A, rho, I, totaldof);

%...reduce global mass matrix and compute A-matrix
M_mod = M(freedof,freedof);
K_mod2 = K_mod(freedof,freedof);
A = M_mod \ K_mod2;
%...solve eigenvalue problem
[qbar,w2]=eig(A);
%...get natural angular frequencies
w=diag(sqrt(w2));
%...convert angular frequencies into frequencies
freq=w/(2*pi);
%...reorder array of frequencies in increasing order
[freq,idx]=sort(freq);
%...and reorder fundamental modes matrix accordingly
qbar=qbar(:,idx);

