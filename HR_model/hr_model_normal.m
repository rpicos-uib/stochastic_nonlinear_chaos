function hr_model_normal()

%    Original equations
% dx = y + phi(x) - z + I
% dy = psi(x) - y
% dz = r * ( s* (x - xr) - z)
% phi(x) = -ax^3 + bx^2
% psi(x) = c - d*x^2
% usual parameters:
% a = 1;
% b = 3;
% c = 1;
% d = 5;
% r = 1e-3;
% s = 4;
% xr = -8/5;
% I = -10:10;

a = 1;
b = 3;
c = 1;
d = 5;
r = 1e-3;
s = 4;
xr = -8/5;
%xr = -1.618;





A=12;

B=1/max([1,A*b, a*A*A,c/A,d*A,r*s,r*s*xr/A]);


% re-scaled variables:
% ai correspond to the dx
% bi correspond to the dy
% ci correspond to the dz

a1=B;
a2=b.*A.*B;
a3=a.*A.^2.*B;
a4=B./A;
a5 = B;
b1=c.*B./A;
b2 = d.*A.*B;
b3 = B;
c1=r.*s.*B;
c2 = r.*s.*B.*xr./A;
c3=r.*B;

% We change the B value to deal with the 1/2 additions
B=B/8;

x1=a2
x2=a3
x3=a1/2
x4=a4/2
x5=a5/2
y1=b2/2
y2=b3/2
y3=b1/4
z1=c2/2
z2=c3/2
z3=c1/4



dt=0.01;
dt = dt/B;

tf=1000/B;
ll=2000000;
ll=ceil(tf/dt)
t=0;
xx=zeros(1,ll);
yy=xx;
zz=xx;

x0 = [0.1,-0.1,3];
x0 = x0/A;

xx(1)=x0(1);
x=xx(1);
yy(1)=x0(2);
y=yy(1);
zz(1)=x0(3);
z=zz(1);


I=3.0;

for iii=1:1


for i=2:ll

    
%     dx = y - a*x*x*x + b*x*x - z + I;
%     dy = c - d*x*x - y;
%     dz = r * ( s* (x - xr) - z);

    F=hindmarsh_rose([x,y,z],I);
    
    x = x + dt*F(1);
    y = y + dt*F(2);
    z = z + dt*F(3);
    
    xx(i)=x;
    yy(i)=y;
    zz(i)=z;
    
end

figure(1)
subplot(1,2,1)
tt=dt*(0:(ll-1))*B;
plot(tt,xx);
hold on
xlabel('time[s]');
ylabel('magnitude');
title('Hindmarsh–Rose neuron model Stimulation')
grid on;
plot(tt,yy)
plot(tt,zz)


% % This is the usual MATLAB integration
% u = 3
% 
% tspan= [0 2000];
% [t,x] = ode45(@(t,x) hindmarsh_rose(x,u), tspan, x0);
% 
% 
% figure(2)
% plot(t,x);
% %plot(t,x(:,1));
% xlabel('time[s]');
% ylabel('magnitude');
% title('Hindmarsh–Rose neuron model Stimulation')
% grid on;



end

function F=hindmarsh_rose(x,u)

X=x(1);
Y=x(2);
Z=x(3);


% Original set of equations
% f1=Y-X.^3+3*X.^2-Z+u.';
% f2=1-5*X.^2-Y;
% f3=0.0005*(4*(X-(-1.618))-Z);


% Original set of equations with the variables scaled to a maximum value of
% 1 and using B to keep all the constants also below 1.
%f1=Y.*B+(b.*A.*B - a.*A.^2.*B.*X).*X.^2-Z.*B+(u'.*B./A);
%f2=c.*B./A-d.*A.*B.*X.^2-B.*Y;
%f3=r.*s.*B*(X-(xr)./A)-r.*B.*Z;

% The same than the above set, but using transformed constants:
% f1 = Y.*a1 + (a2 - a3.*X).*X.^2 - Z.*a5 + u'.*a4;
% f2 = b1 - b2.*X.^2 - b3.*Y;
% f3 = c1*X - c2 - c3.*Z;


% We rewrite the diff eq system in terms of simple
% multiplication and addition functions.
% We also incorporate the transformation of numbers into
% stochastic numbers.
% Notice that, in the current case, since we are with real numbers,
% there is no difference.

% These are the equations using the firts transformed constants
% If you wish to use them, comment out the B=B/8 line before defining the 
% final constants.

% f11=suma(8*a2,mult(-8*a3,X));
% f12=suma(mult(4*a1,Y), mult(f11,mult(X,X)));
% f13=suma(mult(u,a4*4),mult(-4*a5,Z));
% f1=suma(f12,f13);
% f21 = suma(mult(4*b2,mult(X,X)), mult(4*b3,Y));
% f2= suma(2*b1, - f21);
% f31 = suma(4*c2, mult(4*c3,Z));
% f3 = suma(mult(2*c1,X),-f31);

% These are the equations using the FINAL transformation
% Notice that, in order to use them, we have re-scaled also B.
f11=suma(n2s(x1),mult(neg(n2s(x2)),X));
f12=suma(mult(n2s(x3),Y), mult(f11,mult(X,X)));
f13=suma(mult(n2s(u),n2s(x4)),mult(neg(n2s(x5)),Z));
f1=suma(f12,f13);
f21 = suma(mult(n2s(y1),mult(X,X)), mult(n2s(y2),Y));
f2= suma(n2s(y3), neg(f21));
f31 = suma(n2s(z1), mult(n2s(z2),Z));
f3 = suma(mult(n2s(z3),X),neg(f31));



F=[f1;f2;f3];

end
    function s=neg(x)
        s=-x;
    end
    function s=suma(x,y)
        s=0.5*(x+y);
    end
    function s=mult(x,y)
        s=x*y;
    end
    function s=n2s(x)
        s=x;
    end
end