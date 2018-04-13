function [N,L,C,D]=modelschur(nx)
%This program solves a model of the form
%[xt+1][xt]
%B[]=A[]+G[et]
%[Etyt+1][yt]
%using a Schur decomposition of the matrices B and A
%nx is the number of expectational variables in Etyt+1
%if plotcode=1, impulse response is plotted
%solution is yt=-N xt-L et
%and xt+1 = C xt + D et

[kbar,ybar,rbar,cbar,hbar]=steady();
[delta, theta, beta,gamma]=params();

B=[kbar 0 -ybar 0 0
    0 1 0 0 0
    0 -1 theta 0 0 
    0 0 1 0 0
    0 0 0 1 -rbar*beta];
A=[(1-delta)*kbar 0 0 -cbar 0
    0 gamma 0 0 0
    theta 0 0 -(1-theta) 0
    1 0 0 0 1
    0 0 0 1 0
    ];
G=[0 1 0 0 0]';

[AA,BB,Q,Z]=qz(A,B);
[AA,BB,Q,Z] = ordqz(AA,BB,Q,Z,'udi');
%ordeig(AA,BB)   to check eigenvalues 
Zp=Z;
Qp=Q';
[a,b]=size(Z);
N=inv(Zp(a-nx+1:a,a-nx+1:a))*Zp(a-nx+1:a,1:a-nx)
L=inv(Zp(a-nx+1:a,a-nx+1:a))*inv(AA(a-nx+1:a,a-nx+1:a));
L=L*(Qp(a-nx+1:a,1:a-nx)*G(1:a-nx,1)+Qp(a-nx+1:a,a-nx+1:a)*G(a-nx+1:a,1))
invBBN=inv(B(1:a-nx,1:a-nx)-B(1:a-nx,a-nx+1:a)*N);
C=invBBN*(A(1:a-nx,1:a-nx)-A(1:a-nx,a-nx+1:a)*N)
D=invBBN*(G(1:a-nx,1)-A(1:a-nx,a-nx+1:a)*L)


