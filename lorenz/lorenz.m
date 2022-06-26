function lorenz()
close all;
dt=0.5;

x0=0.1;
y0=0.1;
z0=0.1;

x01=0.1;
y01=0.1;
z01=0.1;

%x0=0.875;
%y0=0.875;
%z0=0.875;


w00=0;
w10=0.1;


npunts=50000;

so_x=zeros(npunts,1);
so_y=zeros(npunts,1);
so_z=zeros(npunts,1);

so_w0=zeros(npunts,1);
so_w1=zeros(npunts,1);

so_x1=zeros(npunts,1);
so_y1=zeros(npunts,1);
so_z1=zeros(npunts,1);

so_w01=zeros(npunts,1);
so_w11=zeros(npunts,1);


omega=1.0;

mu=0.81;
alpha=3/8;

x=x0;
y=y0;
z=z0;
t=0;

xb=x01;
yb=y01;
zb=z01;


w0=w00;
w1=w10;

Vmax=1000;

%figure(2)
%hold off;

%plot(x,y,'.r');
%plot(t,x,'.r');
%hold on;

sigma=10;
beta=8/3;
rho=28;
A=50;
D=max([ A rho beta sigma]);
E=4;
%D=1;

p32=2.0^32;
p24=2.0^24;
p16=2.0^16;
p14=2.0^14;
p12=2.0^12;
p10=2.0^10;
p8=2.0^8;
p5=2.0^5;

p1=p24;
p2=p16;
sync=0.0;

for i=1:npunts
%    Original equations
%    dx = y;
%    dy = x - mu*y - x*z;
%    dz = -alpha*z + x*x;
    

precision=p1;

%dx = sigma/D*(y-x);
%dy = x*(rho/D - A/D*z) - y/D;
%dz = A/D*x*y - beta/D*z;

sigmap=n2sn(2*sigma/D/E);
rhop=n2sn(rho/D);
c2 = n2sn(A/D);
c3 = n2sn(2*A/D/E);
c4=n2sn(2/D/E);
c5 = n2sn(2*beta/D/E);
c6=n2sn(4/E);

dx = mult(sigmap, resta(y,x));

dy1 = resta(rhop, mult(c2,z));

dy = resta(mult(mult(c6,x),dy1), mult(c4,y));

dz = resta(mult(c3,mult(x,y)), mult(c5,z));



    x = x + dx*dt;
    y = y + dy*dt;
    z = z + dz*dt;
    t = t+dt;
    
    

    
%    plot((x+2)/4,(y+2)/4,'.r');

    so_x(i)=x;
    so_y(i)=y;
    so_z(i)=z;

    
% Now, let's also implement an oscillator, just to compare...







sigmap=n2sn(2*sigma/D/E);
rhop=n2sn(rho/D);
c2 = n2sn(A/D);
c3 = n2sn(2*A/D/E);
c4=n2sn(2/D/E);
c5 = n2sn(2*beta/D/E);
c6=n2sn(4/E);

precision=p2;

yb=n2sn(sync*y + (1-sync)*yb);

dx = mult(sigmap, resta(yb,xb));

dy1 = resta(rhop, mult(c2,zb));

dy = resta(mult(mult(c6,xb),dy1), mult(c4,yb));

dz = resta(mult(c3,mult(xb,yb)), mult(c5,zb));



    xb = xb + dx*dt;
    yb = yb + dy*dt;
    zb = zb + dz*dt;
    
    

    
%    plot((x+2)/4,(y+2)/4,'.r');

    so_x1(i)=xb;
    so_y1(i)=yb;
    so_z1(i)=zb;









    
end

plot(so_x,'r')
hold on;
plot(so_y,'g')
plot(so_z,'b')

plot(so_x1,'.r')
plot(so_y1,'.g')
plot(so_z1,'.b')

figure(2)
plot(so_y1,so_y,'g')

figure(4)
plot(so_y1-so_y,'.g')

figure(5)
plot3(so_x1-so_x,so_y1-so_y,so_z1-so_z,'r')



figure(3)
plot3(so_x,so_y,so_z,'b')
hold on
plot3(so_x1,so_y1,so_z1,'r')



    function suma=suma(x,y)
        dummy=0.5*(x+y);
        suma= ceil(precision*dummy)/precision;
    end
    
    function resta=resta(x,y)
        resta=suma(x,-y);
    end

    function mult=mult(x,y)
        dummy=x*y;
        mult = ceil(precision*dummy)/precision;
    end
    
    function n=n2sn(x)
        dummy=x;
        n = ceil(precision*dummy)/precision;
    end
    
    end


