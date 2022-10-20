function [tt,xx,yy,zz]=hr_model_normal_1(varargin)
% Parameters:
% tf: total simulation time
% noise: noise level, in % over the signal (see below)

% This functions calculates the Hindmarsh-Rose model
% The variables are re-scaled to be inside the [0,1] range
% The constants are re-scaled to be inside the [0,1] range
% The equations are re-arranged to use the half-addition of Stochastic L.
Nb=30;
lb=2^Nb;

% Let's define a noise to be included in the dx,dy,dz after each
% calculation. This value defines the peak amplitude of the noise as 
noise_d= 0.005;

% Noise mode set multiplicative or additive:
% in the following formula: dx=dx*( (1-noise_d/2) + rand()*noise_d)
% noise_mode=0
% in the following formula: dx=dx + (0.5-noise_d)*rand())
noise_mode=1;

%https://www.mathworks.com/help/matlab/matlab_prog/support-variable-number-of-inputs.html
     if(nargin>0)
        tf=varargin{1};
        if(nargin>1)
            noise_d=varargin{2};
        end
     else
         tf=50;
     end





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

% Original values of the constants:
a = 1;
b = 3;
c = 1;
d = 5;
r = 1e-3;
s = 4;
xr = -8/5;


% Scaling parameters: (X=(x-xa)/xm), where:
% X (uppercase) is the normalized variable and,
% x (lowercase) is the original variable.
xm=6;
ym=14;
zm=0.6;
xa=-2;
ya=-12;
za=2.6;

% These are the initial conditions (at zero time, and so the underscript):
% They will be normalized later on...
xZ=0.1;
yZ=0.1;
zZ=3.0;

% Integration parameters:
%tf=100;    % Final time, in seconds
dt=0.01;    % time interval

% The value of the forcing variable (current):
I=3.0;

% Calculation of the variables for the pre-final form of the equations,
% where the variables are normalized in [0,1].
x0 = 8*(ya/xm-a*xm^2*(xa/xm).^3 + b*xm*(xa/xm).^2-za/xm);
x1= 8*ym/xm;
x2 = 4*a*xm^2;
x3 = 8*(b*xm - 3*a*xa*xm);
x4 = 8*(2*b*xa - 3*a*xa^2);
x5 = 8*zm/xm;
x6 = 8/xm;

y0 = 4*(c/ym -ya/ym - d*(xm^2/ym)*(xa^2/xm^2));
y1 = 4*d*(xm^2/ym);
y2 = 8*d*xm/ym*xa;
y3 = 4;

z1 = 4*r*s*xm/zm;
z0 = 2*(r*s*(xa-xr)/zm - r*za/zm);
z2 = 4*r;

% Now, re-normalize even the constants into the [0,1] interval
% Obviously, we'll need to tweak the time to do that...

Bx = max([x0 x1 x2 x3 x4 x5 x6]);
By = max([y0 y1 y2 y3]);
Bz = max([z0 z1 z2]);

% Just a margin factor of 20%....
B= 1.2*max([Bx By Bz]);

% and we normalize the variables:
x0=x0/B;
x1=x1/B;
x2=x2/B;
x3=x3/B;
x4=-x4/B;
x5=x5/B;
x6=x6/B;
y0=-y0/B;
y1=y1/B;
y2=-y2/B;
y3=y3/B;
z0=-z0/B;
z1=z1/B;
z2=z2/B;


% Now, we absorve the B constant into the time scale:
dt = B*dt;
tf=tf*B;

% We create vectors to store the results:
ll = ceil(tf/dt);
xx=zeros(1,ll);
yy=xx;
zz=xx;

% These are the initial conditions, scaled:
X0 = [(xZ-xa)/xm,(yZ-ya)/ym,(zZ-za)/zm];

xx(1)=X0(1);
x=xx(1);
yy(1)=X0(2);
y=yy(1);
zz(1)=X0(3);
z=zz(1);


for i=2:ll
    
   
    F=hindmarsh_rose([x,y,z],x6*I);
    
%    F=F.*( (1-noise_d/2) + rand(3,1)*noise_d);
    
    if(noise_mode<1)
        x = x*( (1-noise_d/2) + rand()*noise_d) + dt*F(1);
        y = y*( (1-noise_d/2) + rand()*noise_d) + dt*F(2);
        z = z*( (1-noise_d/2) + rand()*noise_d) + dt*F(3);
    else
        x = x + noise_d*(0.5-rand()) + dt*F(1);
        y = y + noise_d*(0.5-rand()) + dt*F(2);
        z = z + noise_d*(0.5-rand()) + dt*F(3);
    end
    
    xx(i)=x;
    yy(i)=y;
    zz(i)=z;
% Just to see the evolution:
%100*i/ll
    
end

figure(1)
subplot(2,1,1)
tt=dt*(0:(ll-1))/B;
plot(tt,xx);
hold on
xlabel('time[s]');
ylabel('magnitude');
title('Hindmarsh–Rose neuron model Simulation, normal')
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



    function F=hindmarsh_rose(x,u)
        
        X=get_sn(x(1));
        Xx=get_sn(x(1));
        Xxx=get_sn(x(1));
        Xx2=get_sn(x(1));
        X2=mult(get_sn(x(1)),get_sn(x(1)));
        X22=mult(get_sn(x(1)),get_sn(x(1)));
        X3=mult(get_sn(x(1)),mult(get_sn(x(1)),get_sn(x(1))));
        
        Y=get_sn(x(2));
        Y1=get_sn(x(2));
        Z=get_sn(x(3));
        Z1=get_sn(x(3));
        
        % Original set of equations with the variables scaled to a maximum 
        % value of 1 and re-ordered to keep up with the semi-additions. 
        % Notice that in order to use them, we have re-scaled also B. 
        % See the beginning ofthe function for the values.
                
%         f11 = suma(suma(x0, mult(x1,Y)), neg(mult(x2,X.^3)));
%         f12 = suma(suma(mult(x3,X.^2), -mult(x4,X)), suma(neg(mult(x5,Z)), mult(x6,u)));
%         f1 = suma(f11,f12);
%         f2 = suma(neg(suma(y0, y1*X.^2)), suma(y2*X, neg(y3*Y)));
%         f3 = suma(neg(z0), suma(z1*X, neg(z2*Z)));

        f11 = suma(suma(get_sn(x0), mult(get_sn(x1),Y)), neg(mult(get_sn(x2),X3)));
        f12 = suma(mult(get_sn(x3),X2), neg(mult(get_sn(x4),Xx)));
        f13 = suma(neg(mult(get_sn(x5),Z)), get_sn(u) );
        f14 = suma(f12,f13);
        f1 = suma(f11,f14);
        f21 = neg(suma(get_sn(y0), mult(get_sn(y1),X22)));
        f22 = suma(mult(get_sn(y2),Xxx), neg(mult(get_sn(y3),Y1)));
        f2 = suma(f21,f22);
        f3 = suma(neg(get_sn(z0)), suma(mult(get_sn(z1),Xx2), neg(mult(get_sn(z2),Z))));
        
        f1=get_n(f1);
        f2=get_n(f2);
        f3=get_n(f3);

        F=[f1;f2;f3];
        
    end
    function s=neg(x)
        s=-x;
    end
    function s=suma(x,y)
        % Just in case x or y are scalars...
%         s=y.*ones(1,lb);
%         sx=x.*ones(1,lb);
%         ii=rand(1,lb);
%         s(ii>0.5)=sx(ii>0.5);
%        s=(x.*du + y.*(1-du));
        s=0.5*(x+y);
    end
    function s=mult(x,y)
         s=x.*y;
%         s=2*mean(s0)-1;
%       s=x.*y;
    end
    function s=n2sn(x)
        s=x;
    end
    function s=get_n(x)
        s=ceil(mean(x)*(2^(Nb-1)))/2^(Nb-1);
    end
    function s=get_sn(x)
%         s=zeros(1,lb);
%         du=rand(1,lb);
%         s(du<=0.5+0.5*x)=1;
        s=x;
    end
end