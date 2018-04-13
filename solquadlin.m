%this program find the value function (P matrix) for the linear quadratic
%optimal regulator problem and also calculates the resulting policy
%function
theta=.36;
beta=.99;
delta=.025;
A=1.72;
kbar=12.6695;
hbar=.3335;
ybar=kbar^theta*hbar^(1-theta);
cbar=ybar-delta*kbar;
aa=(theta*ybar/kbar+1-delta);
a(1,1)=-1/(2*cbar*cbar)*aa*aa-1/(2*cbar)*theta*(1-theta)*ybar/(kbar*kbar);
a(1,2)=1/(2*cbar*cbar)*aa;
a(2,1)=1/(2*cbar*cbar)*aa;
a(1,3)=-1/(2*cbar*cbar)*aa*(1-theta)*ybar/hbar;
a(1,3)=a(1,3)+1/(2*cbar)*theta*(1-theta)*ybar/(kbar*hbar);
a(3,1)=a(1,3);
a(2,2)=-1/(2*cbar*cbar);
a(2,3)=1/(2*cbar*cbar)*(1-theta)*ybar/hbar;
a(3,2)=a(2,3);
a(3,3)=-1/(2*cbar*cbar)*(1-theta)*ybar/hbar*(1-theta)*ybar/hbar;
a(3,3)=a(3,3)-1/(2*cbar)*theta*(1-theta)*ybar/(hbar*hbar);
a(3,3)=a(3,3)-A/(2*(1-hbar)*(1-hbar));
x=[kbar kbar hbar]';
m(1,1)=log(kbar^theta*hbar^(1-theta)-delta*kbar)+A*log(1-hbar);
mm1=1/cbar*(theta*ybar/kbar+1-delta);
mm2=(1-theta)*ybar/(cbar*hbar)-A/(1-hbar);
m(1,1)=m(1,1)-mm1*kbar+kbar/cbar-mm2*hbar;
m(1,1)=m(1,1)+(x')*a*x;
m(1,2)=mm1/2-1*a(1:3,1)'*x;
m(2,1)=m(1,2);
m(1,3)=-1/(2*cbar)-1*a(1:3,2)'*x;
m(3,1)=m(1,3);
m(1,4)=mm2/2-1*a(1:3,3)'*x;
m(4,1)=m(1,4);
m(2:4,2:4)=a;
AA=[1 0
    0 0];
B=[0 0
    1 0];
R=m(1:2,1:2);
Q=m(3:4,3:4);
W=m(1:2,3:4)';
P=[1 0
    0 1];
%iterating the Ricotti equation
for i=1:1000
    zinv=inv(Q+beta*B'*P*B);
    z2=beta*AA'*P*B+W';
    P=R+beta*AA'*P*AA-z2*zinv*z2';
end

%finding the policy function
P=P
F=-zinv*(W+beta*B'*P*AA)



