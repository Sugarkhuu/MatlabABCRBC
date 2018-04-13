%stationary state values are found in another program
[kbar,ybar,rbar,cbar,hbar]=steady();
[delta, theta, beta]=params();
A=[0 -kbar 0 0]';
B=[0 (1-delta)*kbar theta -1]';
C=[1 -1 -1 0
    ybar -cbar 0 0
    -1 0 1-theta 0
    1 0 0 -1];
D=[0 0 1 0]';
F=[0];
G=F;
H=F;
J=[0 -1 0 beta*rbar];
K=[0 1 0 0];
L=F;
M=F;
N=[.95];

Cinv=inv(C);
a=F-J*Cinv*A;
b=-(J*Cinv*B-G+K*Cinv*A);
c=-K*Cinv*B+H;
P1=(-b+sqrt(b^2-4*a*c))/(2*a);
P2=(-b-sqrt(b^2-4*a*c))/(2*a);
if abs(P1)<1
    P=P1
else
    P=P2
end
R=-Cinv*(A*P+B)
Q=(J*Cinv*D-L)*N+K*Cinv*D-M;
QD=kron(N',(F-J*Cinv*A))+(J*R+F*P+G-K*Cinv*A);
Q=Q/QD
S=-Cinv*(A*Q+D)

