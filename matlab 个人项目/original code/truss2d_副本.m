function truss2d
%
%EG-225: Structural Mechanics IIb
%Matrix analysis of trusses
%Program created by Dr. A.J.Gil
%
clc;close all;
%input data
[ndof,X,connectivity,E,A,rho,freedof,P,fixdof,Uf,axialunits,lengthunits]=readinputdata;
%matrix allocation
F=zeros(ndof*size(X,1),1);F(freedof)=P;U=zeros(ndof*size(X,1),1);U(fixdof)=Uf;
%assemble global stiffness matrix
[K]=compute_stiffness_matrix(ndof,X,connectivity,E,A);

%eigenvalue analysis (prior to BCs)
eigenvalue_analysis=1;
if eigenvalue_analysis==1
   [V,D] = eig(K);
   figure(1);factor=1e3;plotmodes(connectivity,X,factor*V,D,lengthunits);
end
%eigenvalue analysis (after BCs)
eigenvalue_analysis=1;
if eigenvalue_analysis==1
   [V,D] = eig(K(freedof,freedof));
   figure(2);factor=1e3;
   plotmodes_2(connectivity,X,ndof,freedof,fixdof,Uf,factor*V,D,lengthunits);
end


%assemble global stiffness matrix
[M]=compute_mass_matrix(ndof,X,connectivity,rho,A);
%eigenvalue analysis (after BCs and dynamic)
eigenvalue_analysis=1;
if eigenvalue_analysis==1
   [V,D] = eig(K(freedof,freedof),M(freedof,freedof));
   figure(3);factor=1e3;
   plotmodes_2(connectivity,X,ndof,freedof,fixdof,Uf,factor*V,D,lengthunits);
end




%computing unknown displacements
U(freedof)=K(freedof,freedof)\(F(freedof)-K(freedof,fixdof)*U(fixdof));

%computing unknown reactions
F(fixdof)=K(fixdof,freedof)*U(freedof)+K(fixdof,fixdof)*U(fixdof);
%computing internal forces
Finternal=compute_internal_forces(ndof,X,connectivity,E,A,U);
%print results in a tabular format
printresults(X,U,F,Finternal(:,2,1),axialunits,lengthunits)
%plot axial forces
figure(4); plotaxialforces(connectivity,X,A,Finternal(:,2,1),axialunits,lengthunits)
%plot displacements
figure(5); plotdisplacements_truss(connectivity,X,U,lengthunits)
%
end


function [K]=compute_stiffness_matrix(ndof,X,connectivity,E,A)
%initial allocation
K=zeros(ndof*size(X,1),ndof*size(X,1)); 
%loop over elements
for ielem=1:size(connectivity,1)
    %extract connectivity information
    nodei=connectivity(ielem,1);idof=ndof*(nodei-1)+1:ndof*(nodei-1)+ndof;
    nodej=connectivity(ielem,2);jdof=ndof*(nodej-1)+1:ndof*(nodej-1)+ndof;
    %compute geometric information
    Xij=X(nodej,:)-X(nodei,:);
    L=norm(Xij,2); c=Xij(1)/L; s=Xij(2)/L; 
    T=[c s];
    %compute element global stiffness block (diagonal)
    Ke=E(ielem)*A(ielem)/L*(T'*T);
    %assemble global stiffness matrix
    K(idof,idof)=K(idof,idof)+Ke; K(idof,jdof)=K(idof,jdof)-Ke;
    K(jdof,idof)=K(jdof,idof)-Ke; K(jdof,jdof)=K(jdof,jdof)+Ke;
end
end

function [M]=compute_mass_matrix(ndof,X,connectivity,rho,A)
%initial allocation
M=zeros(ndof*size(X,1),ndof*size(X,1)); 
%loop over elements
for ielem=1:size(connectivity,1)
    %extract connectivity information
    nodei=connectivity(ielem,1);idof=ndof*(nodei-1)+1:ndof*(nodei-1)+ndof;
    nodej=connectivity(ielem,2);jdof=ndof*(nodej-1)+1:ndof*(nodej-1)+ndof;
    %compute geometric information
    Xij=X(nodej,:)-X(nodei,:);
    L=norm(Xij,2); 
    %compute element global stiffness block (diagonal)
    m=rho(ielem)*A(ielem)/L/6;
    %assemble global stiffness matrix
    M(idof,idof)=M(idof,idof)+2*m*eye(ndof); M(idof,jdof)=M(idof,jdof)+m*eye(ndof);
    M(jdof,idof)=M(jdof,idof)+m*eye(ndof);   M(jdof,jdof)=M(jdof,jdof)+2*m*eye(ndof);
end
end

function Finternal=compute_internal_forces(ndof,X,connectivity,E,A,U)
%initial allocation
Finternal=zeros(size(connectivity,1),2,1);
%loop over elements
for ielem=1:size(connectivity,1)
    %extract connectivity information
    nodei=connectivity(ielem,1);idof=ndof*(nodei-1)+1:ndof*(nodei-1)+ndof;
    nodej=connectivity(ielem,2);jdof=ndof*(nodej-1)+1:ndof*(nodej-1)+ndof;
    %compute geometric information
    Xij=X(nodej,:)-X(nodei,:);
    L=norm(Xij,2); c=Xij(1)/L; s=Xij(2)/L; 
    T=[c s];
    %compute element global stiffness block (diagonal) 
    Ke=E(ielem)*A(ielem)/L*(T'*T);
    %Node 1
    Finternal(ielem,1,1)=-T*( Ke*U(idof,1) -Ke*U(jdof,1));
    %node 2
    Finternal(ielem,2,1)= T*(-Ke*U(idof,1) +Ke*U(jdof,1));
end
end

function printresults(X,U,F,N,axialunits,lengthunits)

fprintf('      Table results \n');
fprintf('-------------------------   \n');
fprintf('   \n');

fprintf(['  Node             U' lengthunits '       V' lengthunits '\n']);
for i=1:size(X,1)
    idof=2*(i-1)+1:2*(i-1)+2;
    fprintf('   %2g       %10.5g  %10.5g\n', i, U(idof,:));
end
fprintf(['  Node             FX' axialunits '       FY' axialunits '\n']);
for i=1:size(X,1)
    idof=2*(i-1)+1:2*(i-1)+2;
    fprintf('   %2g       %10.5g  %10.5g\n', i, F(idof,:));
end
fprintf(['  Element   Axial Force' axialunits '\n']);
for ielem=1:size(N,1)
    fprintf('   %2g       %10.5g\n', ielem, N(ielem));
end

end

function plotdisplacements_truss(connectivity,X,U,units)
%
nelem=size(connectivity,1);ndof=2;
for ielem=1:nelem
    nodei=connectivity(ielem,1);idof=ndof*(nodei-1)+1:ndof*(nodei-1)+ndof;
    nodej=connectivity(ielem,2);jdof=ndof*(nodej-1)+1:ndof*(nodej-1)+ndof;
    Xi=X(nodei,:);  Xj=X(nodej,:);
    xi=Xi+U(idof)'; xj=Xj+U(jdof)';
    plot([Xi(1),Xj(1)],[Xi(2),Xj(2)],'b',[xi(1),xj(1)],[xi(2),xj(2)],'r') 
    hold on;
end
grid on;xlabel(['OX axis ' units]);ylabel(['OY axis ' units]); 
title(['2D truss: displacements result ' units]);
%
end

function plotaxialforces(connectivity,X,A,N,axialunits,lengthunits)
%
n=10;maxwidth=6;maxcolor=100;map=colormap(jet(maxcolor));
nelem=size(connectivity,1);
for ielem=1:nelem
    nodei=connectivity(ielem,1);nodej=connectivity(ielem,2);
    Xi=X(nodei,:); Xj=X(nodej,:);
    x=linspace(Xi(1),Xj(1),n);y=linspace(Xi(2),Xj(2),n);
    width=A(ielem)/max(A)*maxwidth;
    color=1+floor((N(ielem)-min(N))/(max(N)-min(N))*(maxcolor-1));
    plot(x,y,'LineWidth',width,'Color',map(color,:));
    hold on;
end
grid on;
xlabel(['OX axis ' lengthunits]);ylabel(['OY axis ' lengthunits]); 
caxis([min(N) max(N)]);colorbar('location','eastoutside');
title(['2D truss: axial force result ' axialunits]);
%
end

function plotmodes(connectivity,X,V,D,lengthunits)
%
nelem=size(connectivity,1);
for mode=1:size(V,2)
    U=V(:,mode);
    subplot(size(V,2),1,mode);
    for ielem=1:nelem
        nodei=connectivity(ielem,1);idof=[2*(nodei-1)+1:2*(nodei-1)+2];
        nodej=connectivity(ielem,2);jdof=[2*(nodej-1)+1:2*(nodej-1)+2];
        Xi=X(nodei,:);  Xj=X(nodej,:);
        xi=Xi+U(idof)'; xj=Xj+U(jdof)';
        plot([Xi(1),Xj(1)],[Xi(2),Xj(2)],'b',[xi(1),xj(1)],[xi(2),xj(2)],'r')
        hold on;
    end
    grid on;axis equal;xlabel(['OX axis ' lengthunits]);ylabel(['OY axis ' lengthunits]);
    title(['Mode ' num2str(mode) '   \lambda= ' num2str(D(mode,mode))]);
end
%
end

function plotmodes_2(connectivity,X,ndof,freedof,fixdof,Uf,V,D,lengthunits)
%
nelem=size(connectivity,1);
for mode=1:size(V,2)
    U=zeros(ndof*size(X,1),1);
    U(fixdof)=Uf;
    U(freedof)=V(:,mode);
    subplot(size(V,2),1,mode);
    for ielem=1:nelem
        nodei=connectivity(ielem,1);idof=[2*(nodei-1)+1:2*(nodei-1)+2];
        nodej=connectivity(ielem,2);jdof=[2*(nodej-1)+1:2*(nodej-1)+2];
        Xi=X(nodei,:);  Xj=X(nodej,:);
        xi=Xi+U(idof)'; xj=Xj+U(jdof)';
        plot([Xi(1),Xj(1)],[Xi(2),Xj(2)],'b',[xi(1),xj(1)],[xi(2),xj(2)],'r')
        hold on;
    end
    grid on;axis equal;xlabel(['OX axis ' lengthunits]);ylabel(['OY axis ' lengthunits]);
    title(['Mode ' num2str(mode) '   \lambda= ' num2str(D(mode,mode))]);
end
%
end


function [ndof,X,connectivity,E,A,rho,freedof,P,fixdof,Uf,axialunits,lengthunits]=readinputdata
%
% ndof=2;
% X=[0 0; 1e4 0; 1e4 1e4; 0 1e4];
% connectivity=[1 2; 1 3; 1 4; 3 4; 2 4; 2 3];
% E=2e1*[1; 1; 1; 1; 1; 1];
% A=484*[1; 1; 1; 1; 1; 1];
% freedof=[1 2 7 8]; P=[1000 0 0 0];
% fixdof=[3 4 5 6]; Uf=[0 0 0 0];
% axialunits='(kN)';
% lengthunits='(mm)';
%

%
% ndof=2;
% X=[0 5*sqrt(3); 0 5*(sqrt(3)-1); 0 0; 5 5*sqrt(3)]*1e3;
% connectivity=[1 4; 2 4; 3 4];
% E=[1 1 1];
% A=[1 2 1]*1e4;
% freedof=[7 8]; P=[0 -1];
% fixdof=[1 2 3 4 5 6]; Uf=[0 0 0 0 0 0];
% axialunits='(kN)';
% lengthunits='(mm)';
%

ndof=2;
X=[0 6; 0 0; 3*(1+sqrt(3)) 0; 3*(1+sqrt(3)) 6; 3*sqrt(3) 3]*1e3;
connectivity=[1 5; 2 5; 4 5; 3 5];
E=2e2*[1 1 1 1];
A=[5 5 3 3 ]*1e3;
rho=1*[1 1 1 1];
%freedof=[9 10]; P=[1 0]*1e3;
%fixdof=[1 2 3 4 5 6 7 8]; Uf=[0 0 0 0 0 0 0 0];
freedof=[5 6 8 9 10]; P=[-1 1 1 1 0]*1e3;
fixdof=[1 2 3 4 7]; Uf=[0 0 0 0 0];

axialunits='(kN)';
lengthunits='(mm)';
%


end





