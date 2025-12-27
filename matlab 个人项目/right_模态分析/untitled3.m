%...reduce global mass matrix and compute A-matrix
Mreduced = M(freedoftable,freedoftable);
A = Mreduced \ Kreduced;
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