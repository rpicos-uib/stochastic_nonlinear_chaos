function shimizu_morioka()
dt=0.1;

x0=0.51;
y0=0.51;
z0=0.51;

%x0=0.875;
%y0=0.875;
%z0=0.875;


w00=0;
w10=0.1;


npunts=1000000;

npunts=284078;

so_x=zeros(npunts,1);
so_y=zeros(npunts,1);
so_z=zeros(npunts,1);

so_w0=zeros(npunts,1);
so_w1=zeros(npunts,1);


omega=1.0;

mu=0.81;
alpha=3/8;

x=x0;
y=y0;
z=z0;
t=0;


w0=w00;
w1=w10;

Vmax=1;

figure(2)
hold off;

%plot(x,y,'.r');
%plot(t,x,'.r');
hold on;

for i=1:npunts
%    Original equations
%    dx = y;
%    dy = x - mu*y - x*z;
%    dz = -alpha*z + x*x;
    

%   Normalized to (0--1)
% x'=(x+xmin)/xmax
%    dx = y - 1/2;
%    dy = 4*x -2 +mu/2 - mu*y +3*z - 6*x*z;
%    dz = -alpha*z + alpha/2 + 8/3*x*x +2/3 - 8/3*x;
    
%   Normalized with t to have everything <1
% t'=8t
%    dx = y/8 - 1/16;
%    dy = x/2 -1/4 + mu/16 -mu*y/8 + 3*z/8 -3/4*x*z;
%    dz = -alpha*z/8 + alpha/16 + 1/3*x*x + 1/12 - x/3;

% t'=32t
%    dx = y/32 - 1/64;
%    dy = 1/2*(1/2*(1/2*(x-1/2)+ 1/2*(mu/8 + 3/4*z) )) - 1/2*(1/2*(3/4*x*z + mu*y/8));
%    dz = -alpha*z/32 + alpha/64 + 1/12*x*x + 1/48 - x/12;
% dx


%    dx = y/8 - 1/16; = 1/2(y/4 - 1/8)
sn_x = n2sn(x/Vmax);
sn_x2 = n2sn(x/Vmax);
sn_y = n2sn(y/Vmax);
sn_z = n2sn(z/Vmax);
c1 = n2sn(1/16);
c2 = n2sn(1/32);
dx = resta(mult(sn_y,c1), c2 );


%    dy = x/2 -1/4 + mu/16 + 3*z/8 -3/4*x*z -mu*y/8;
%       = 1/2(x-1/2) + 1/2(mu/8 + z*3/4) - 1/2( 3/2*x*z + mu*y/4)
sn_x = n2sn(x/Vmax);
sn_x2 = n2sn(x/Vmax);
sn_y = n2sn(y/Vmax);
sn_z = n2sn(z/Vmax);
c1=n2sn(1/2);
c2=n2sn(mu/8);
c3=n2sn(3/4);
c4=n2sn(3/4);
c5=n2sn(mu/8);

dy1= suma(resta(sn_x, c1), suma(c2, mult(c3, sn_z)));
sn_x = n2sn(x/Vmax);
sn_x2 = n2sn(x/Vmax);
sn_y = n2sn(y/Vmax);
sn_z = n2sn(z/Vmax);
dy2 = suma(mult(c4,mult(sn_x,sn_z)), mult(c5,sn_y));
dy = resta(dy1, dy2);

%    dz = -x/3 -alpha*z/8 + alpha/16 + 1/3*x*x + 1/12;
c1 = n2sn(1/2);
c2 = n2sn(alpha/4);
c3 = n2sn(alpha/8);
c4 = n2sn(2/3);
c5 = n2sn(1/24);
c6 = n2sn(2/3);

sn_x = n2sn(x/Vmax);
sn_x2 = n2sn(x/Vmax);
sn_y = n2sn(y/Vmax);
sn_z = n2sn(z/Vmax);

dz1 = suma(mult(c4,sn_x), mult(c2, sn_z));
sn_x = n2sn(x/Vmax);
sn_x2 = n2sn(x/Vmax);
sn_y = n2sn(y/Vmax);
sn_z = n2sn(z/Vmax);
dz2 = suma(c3, mult(c6, mult(sn_x, sn_x2) ) );
dz3 = resta(dz2, dz1);
dz = suma(dz3, c5);



    x = x + dx*dt;
    y = y + dy*dt;
    z = z + dz*dt;
    t = t+dt;
    
    
    if(x>Vmax) 
        x=0;
    end
    if(x<0) 
        x=Vmax;
    end
        
    if(y>Vmax) 
        y=0;
    end
    if(y<0) 
        y=Vmax;
    end
    if(z>Vmax) 
        z=0;
    end
    if(z<0) 
        z=Vmax;
    end

    
%    plot((x+2)/4,(y+2)/4,'.r');

    so_x(i)=x;
    so_y(i)=y;
    so_z(i)=z;

    
% Now, let's also implement an oscillator, just to compare...


    dw0 = mult(1,w1);
    dw1 = -w0;

    
    w0 = w0 + dw0*dt;
    w1 = w1 + dw1*dt;

    so_w0(i)=w0;
    so_w1(i)=w1;

    if(w0>Vmax) 
        w0=Vmax;
    end
    if(w0<-Vmax) 
        w0=-Vmax;
    end
        
    if(w1>Vmax) 
        w1=Vmax;
    end
    if(w1<-Vmax) 
        w1=-Vmax;
    end

    
end

plot(so_x,'r')
hold on;
plot(so_y,'g')
plot(so_z,'b')

legend('x','y','z');
xlabel('Iteration');

%plot(so_w0,'k')
%plot(so_w1,'k.-')



figure(3)
plot3(so_x,so_y,so_z)
box on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Shimizu-Morioka');










xx2 = so_x;
xf2 = fft(xx2);

L=length(xf2);

P2 = abs(xf2/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

Fs=1;
f = Fs*(0:(L/2))/L;

figure(5)

m=max(P1);

hold on;
semilogx(f,P1/m,'b') 

legend('SL','Usual Integration')
xlabel('Frequency (a.u.)');
ylabel('Amplitude (a.u.)')



    function suma=suma(x,y)
        suma=0.5*(x+y);
    end
    
    function resta=resta(x,y)
        resta=suma(x,-y);
    end

    function mult=mult(x,y)
        mult=x*y;
    end
    
    function n=n2sn(x)
        n=x;
    end
    
    end


