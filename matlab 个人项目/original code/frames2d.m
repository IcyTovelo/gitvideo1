function frames2d
%
clc;close all;
%input data
[ndof,X,connectivity,E,A,I,freedof,P,fixdof,Uf,forceunits,lengthunits]=readinputdata;
%matrix allocation
F=zeros(ndof*size(X,1),1);F(freedof)=P;U=zeros(ndof*size(X,1),1);U(fixdof)=Uf;
%assemble global stiffness matrix
[K]=compute_stiffness_matrix(ndof,X,connectivity,E,A,I);
%computing unknown displacements
U(freedof)=K(freedof,freedof)\(F(freedof)-K(freedof,fixdof)*U(fixdof));
%computing unknown reactions
F(fixdof)=K(fixdof,freedof)*U(freedof)+K(fixdof,fixdof)*U(fixdof);
%computing internal forces
Finternal=compute_internal_forces(ndof,X,connectivity,E,A,I,U);
%print results in a tabular format
printresults(X,U,F,Finternal,forceunits,lengthunits)
%plot axial forces
figure(1); factor=1;offset=0;plot_internal_forces(X,connectivity,Finternal,factor,offset);
%plot displacements
figure(2); factor=1;plotdisplacements_frame(X,connectivity,U,factor,0);
%
end


function [K]=compute_stiffness_matrix(ndof,X,connectivity,E,A,I)
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
end
end

function Finternal=compute_internal_forces(ndof,X,connectivity,E,A,I,U)
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
    %frame transformation matrix
    T=[c s 0; -s c 0; 0 0 1];
    %local displacements
    ui=T*U(idof,1);uj=T*U(jdof,1);
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
  
    %compute local internal forces          
    fi=kii*ui+kij*uj;fj=kji*ui+kjj*uj;
    %extract axial, shear and bending moment for Node 1
    Finternal(ielem,1,1)=-fi(1,1);
    Finternal(ielem,1,2)= fi(2,1);
    Finternal(ielem,1,3)=-fi(3,1);
    %extract axial, shear and bending moment for Node 2
    Finternal(ielem,2,1)= fj(1,1);
    Finternal(ielem,2,2)=-fj(2,1);
    Finternal(ielem,2,3)= fj(3,1);
end
end

function printresults(X,U,F,Finternal,forceunits,lengthunits)
%
fprintf('      Table results \n');
fprintf('-------------------------   \n');
fprintf('   \n');

fprintf(['  Node             U' lengthunits '       V' lengthunits '   Theta' ' \n']);
for i=1:size(X,1)
    idof=3*(i-1)+1:3*(i-1)+3;
    fprintf('   %2g       %10.5g  %10.5g  %10.5g\n', i, U(idof,:));
end
fprintf(['\n  Node             FX' forceunits '      FY' forceunits '    M' forceunits '*' lengthunits '\n']);
for i=1:size(X,1)
    idof=3*(i-1)+1:3*(i-1)+3;
    fprintf('   %2g       %10.5g  %10.5g  %10.5g\n', i, F(idof,:));
end
fprintf(['\n  Internal Forces' forceunits '\n\n']);

%Axial force results
k=1;
fprintf(['  Axial Force' forceunits '\n']);
fprintf('  Element        Node1        Node2 \n');
for ielem=1:size(Finternal,1)
    fprintf('   %2g       %10.5g    %10.5g\n', ielem, Finternal(ielem,1,k), Finternal(ielem,2,k));
end

%Shear force results
k=2;
fprintf(['  Shear Force' forceunits '\n']);
fprintf('  Element        Node1        Node2 \n');
for ielem=1:size(Finternal,1)
    fprintf('   %2g       %10.5g    %10.5g\n', ielem, Finternal(ielem,1,k), Finternal(ielem,2,k));
end

%Bending moment results
k=3;
fprintf(['  Bending moment' forceunits '*' lengthunits '\n']);
fprintf('  Element        Node1        Node2 \n');
for ielem=1:size(Finternal,1)
    fprintf('   %2g       %10.5g    %10.5g\n', ielem, Finternal(ielem,1,k), Finternal(ielem,2,k));
end
 %  
end

function plotdisplacements_frame(X,connectivity,U,factor,truss)
%
for ielem=1:size(connectivity,1)
    %initial shape
    nodei=connectivity(ielem,1);nodej=connectivity(ielem,2);
    Xij=X(nodej,:)-X(nodei,:);
    L=norm(Xij);n=Xij/L;
    plot([X(nodei,1) X(nodej,1)],[X(nodei,2) X(nodej,2)],'k-','LineWidth',1);
    hold on;
    %deformed shape
    Ui=factor*U(3*(nodei-1)+1:3*nodei);Uj=factor*U(3*(nodej-1)+1:3*nodej);
    T=[n(1) n(2) 0;-n(2) n(1) 0; 0 0 1];
    ui=T*Ui;uj=T*Uj;
    if truss
       theta=(uj(2)-ui(2))/(L+uj(1)-ui(1));
    else
       theta=0; 
    end
    xlocal=[ui(1) uj(1)]+[0 L]; ylocal=[ui(2) uj(2)];   
    cubic=spline(xlocal,[ui(3)+theta ylocal uj(3)+theta]);    
    x_vec = linspace(xlocal(1),xlocal(2),25)';
    y_vec = ppval(cubic,x_vec);
    x=X(nodei,1)+n(1)*x_vec-n(2)*y_vec;
    y=X(nodei,2)+n(2)*x_vec+n(1)*y_vec;
    plot(x,y,'-b','LineWidth',2);grid on;
    hold on;      
end
title('Displaced shape');xlabel('OX axis');ylabel('OY axis');
%
end

function plot_internal_forces(X,connectivity,Int_forces,factor,tol)
%
for ielem=1:size(connectivity,1)
    nodei=connectivity(ielem,1);nodej=connectivity(ielem,2);
    Xij=X(nodej,:)-X(nodei,:);
    n=Xij/norm(Xij);v=[-n(2) n(1)];
    for k=1:3
        if k==1
            style='b-';graph_name='Axial force diagram';
        elseif k==2
            style='r-';graph_name='Shear force diagram';
        else
            style='g-';graph_name='Bending moment diagram';
        end
        subplot(1,3,k);
        
        plot([X(nodei,1)+factor*Int_forces(ielem,1,k)*v(1) X(nodej,1)+factor*Int_forces(ielem,2,k)*v(1)],...
             [X(nodei,2)+factor*Int_forces(ielem,1,k)*v(2) X(nodej,2)+factor*Int_forces(ielem,2,k)*v(2)],style,'LineWidth',2);
        hold on;
        plot([X(nodei,1) X(nodei,1)+factor*Int_forces(ielem,1,k)*v(1) ],...
            [X(nodei,2) X(nodei,2)+factor*Int_forces(ielem,1,k)*v(2) ],style,'LineWidth',2);
        hold on;
        plot([X(nodej,1) X(nodej,1)+factor*Int_forces(ielem,2,k)*v(1)],...
            [X(nodej,2) X(nodej,2)+factor*Int_forces(ielem,2,k)*v(2)],style,'LineWidth',2);
        hold on;
        plot([X(nodei,1) X(nodej,1)],[X(nodei,2) X(nodej,2)],'k-','LineWidth',1);
        hold on;
        
        text(X(nodei,1)+factor*Int_forces(ielem,1,k)*v(1),X(nodei,2)+factor*Int_forces(ielem,1,k)*v(2),num2str(Int_forces(ielem,1,k)),'FontSize',10);
        text(X(nodej,1)+factor*Int_forces(ielem,2,k)*v(1),X(nodej,2)+factor*Int_forces(ielem,2,k)*v(2),num2str(Int_forces(ielem,2,k)),'FontSize',10);
         
          
        title(graph_name);
        grid on;
        axis equal;
        Xmin=min(X(:,1));Xmax=max(X(:,1));dX=Xmax-Xmin;Ymin=min(X(:,2));Ymax=max(X(:,2));dY=Ymax-Ymin;
        if abs(dX)>1.0e-4 && abs(dY)>1.0e-4
           axis([Xmin-tol*dX Xmax+tol*dX Ymin-tol*dY Ymax+tol*dY]);
        end
        
    end 
end
%
end

function [ndof,X,connectivity,E,A,I,freedof,P,fixdof,Uf,forceunits,lengthunits]=readinputdata
%

%frame problem
% ndof=3;
% X=[0 0; 0 1; 1 1; 1 0];
% connectivity=[1 2; 2 3; 3 4];
% E=1*[1; 1; 1];
% A=1e4*[1; 1; 1];
% I=1*[1; 1; 1];
% freedof=[4 5 6 7 8 9]; P=[0 0 -1 0 0 1];
% fixdof=[1 2 3 10 11 12]; Uf=[0 0 0 0 0 0];
% forceunits='(kN)';
% lengthunits='(m)';
%
%beam problem: exercise 6
ndof=3;
X=[0 0; 4 0; 7 0; 10 0];
connectivity=[1 2; 2 3; 3 4];
E=1*[1; 1; 1];
A=1e4*[1; 1; 1];
I=1e5*[1; 1; 1];
freedof=[4 6 7 8 9]; P=[0 0 0 -10 0];
fixdof=[1 2 3 5 10 11 12]; Uf=[0 0 0 0 0 0 0];
forceunits='(kN)';
lengthunits='(m)';
%
% ndof=3;
% X=[0 0; 1 0; 2 0; 0 1; 0 -1];
% connectivity=[1 2; 2 3; 1 4; 1 5];
% E=1*[1; 1; 1; 1];
% A=1e4*[1; 1; 1; 1];
% I=1*[1; 1; 1; 1];
% freedof=[1 2 3 6 9]; P=[0 0 -1 -1 2];
% fixdof=[4 5 7 8 10 11 12 13 14 15]; Uf=[0 0 0 0 0 0 0 0 0 0];
% forceunits='(kN)';
% lengthunits='(m)';
%

% ndof=3;
% X=[0 0; 0.5 0; 1 0; 1.5 0; 2 0; 0 1; 0 -1];
% connectivity=[1 2; 2 3; 3 4; 4 5;1 6; 1 7];
% E=1*[1; 1; 1; 1; 1; 1];
% A=1e4*[1; 1; 1; 1; 1; 1];
% I=1*[1; 1; 1; 1; 1; 1];
% freedof=[1 2 3 4 5 6 9 10 11 12 15]; P=[0 0 0 0 -5 0 0 0 -10 0 0];
% fixdof=[7 8 13 14 16 17 18 19 20 21]; Uf=0*fixdof;
% forceunits='(kN)';
% lengthunits='(m)';

%beam problem: exercise 4
ndof=3;
X=[0 0; 2 0; 4 0];
connectivity=[1 2; 2 3];
E=1*[1; 1];
A=1e4*[1; 1];
I=1e4*[1; 1];
freedof=[3 4 6 9]; P=[0 0 -8  5];
fixdof=[1 2 5 7 8]; Uf=[0 0 0 0 0];
forceunits='(kN)';
lengthunits='(m)';
% % %

%frame problem
% ndof=3;
% X=[0 0; 0 1; 1 1; 1 0];
% connectivity=[1 2; 2 3; 3 4];
% E=1*[1; 1; 1];
% A=1*[1; 1; 1];
% I=10*[1; 1; 1];
% freedof=[4 5 6 7 8 9]; P=[0 0 0 0 0 10];
% fixdof=[1 2 3 10 11 12]; Uf=[0 0 0 0 0 0];
% forceunits='(kN)';
% lengthunits='(m)';

%frame problem
%  ndof=3;
%  X=[0 0; 2 0; 4 0];
%  connectivity=[1 2; 2 3];
%  E=1*[1; 1];
%  A=0*[1; 1];
%  I=1e4*[1; 1];
%  freedof=[6 8 9]; P=[0 -5 0];
%  fixdof=[1 2 3 4 5 7]; Uf=[0 0 0 0 0 0];
%  forceunits='(kN)';
%  lengthunits='(m)';
 
 
 %frame problem
%  ndof=3;
%  %X=1e3*[0 0; sqrt(3) 3; sqrt(3)+3 3; sqrt(3)+3 0];
%  X=1e3*[0 0; 1.732 3; 1.732+3 3; 1.732+3 0];
%  connectivity=[1 2; 2 3; 3 4];
%  E=200*[1; 1; 1];
%  A=1e2*[32; 32; 66];
%  I=1e4*[2356; 2356; 5245];
%  freedof=[4 5 6 7 8 9 ]; P=[0 -20 0 10 0 0];
%  fixdof=[1 2 3 10 11 12]; Uf=[0 0 0 0 0 0];
%  forceunits='(kN)';
%  lengthunits='(m)';


end





