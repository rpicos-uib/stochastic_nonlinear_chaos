function zhang_tang_s(dibuixa2)

if (nargin==0)
  dibuixa=1;
else
  dibuixa=dibuixa2;
end

llindar = 0.5

fitxer='zhang_tang.txt';


fileID=fopen(fitxer,'a');

%iteracions1 = 2^28;
%iteracions2 = 2^26;
%Vmax = 2 ^26;
%NNN=100;


% Test
iteracions1 = 2^20;
iteracions2 = 2^24;
Vmax = 2 ^24;
NNN=iteracions1/50;
increment=2^0;

% Thus, we draw 100 points each time..
% We iterate iteracions2 and we get ONE point.
% We iterate iteracions1/NNN to get an array of points
% We iterate NNN times to get all the iterations
% increment is related to dt= increment/Vmax;

% Condicions inicials per continuar:
% Llegim el fitxer, i continuam allÃ  on erem....

%if(0<1)
    try
        comanda=['tail -n 1 ' fitxer]
        [status result] = system(comanda);
        if isstrprop(result(end), 'cntrl'), result(end) = []; end
        n=sscanf(result,'%f %f %f %f %f')
        n0=n(1)
        x=n(2)
        y=n(3)
        z=n(4)
        q=n(4)
        string="Result gotten"
    catch
        x = Vmax *0.51;
        y = Vmax *0.51;
        z = Vmax *0.51;
        q = Vmax *0.51;
        n0=1;
    end
%end

x_o=zeros(iteracions1,1);
y_o=x_o;
z_o=x_o;
q_o=x_o;






dt=0.01;
t=0;


npunts=5000000;

so_x=zeros(npunts,1);
so_y=zeros(npunts,1);
so_z=zeros(npunts,1);
so_q=zeros(npunts,1);

sortida1=0;
sortida2=0;

inc=1;
t1=1;

%Vmax=1;

close all;
%figure(2)
%hold off;

%plot(t,x1,'.r');
%hold on;

% Provarem es coeficients de sa Fig3
a1=-0.3;
a2=-0.5;
a3=-0.6;
a4=-0.1;
a5=-0.1;
a6=-0.65;
a7=-0.1;

b1=0.8;
b2=1.5;
b3=3;
b4=0.6;



% Provarem es coeficients de sa Fig4
a1=-0.3;
a2=-0.5;
a3=-0.6;
a4=-0.1;
a5=-0.1;
a6=-0.65;
a7=-0.1;

b1=0.8;
b2=1.5;
b3=1;
b4=0.6;

% Initial conditions
x0=0.0;
y0=-0.0010;
z0=-0.0001;
q0=0.00001;




x1 = x0;
x2 = y0;
x3 = z0;
x4 = q0;

% Formules equilibri article

A=(a1*a7*b1) + (a2*a3*b4);

B=(a2*a2*b4) + (a2*a6*b1) + (a3*a7);

A=50;
B=6*A;



for cc=1:1000


for c= 0:(NNN-1)


for a = 1:iteracions1/NNN
    
for b = 1:iteracions2


%for i=1:npunts

    %    Original equations
%    x1=a1x1 + a2x4 - x2x3;
%    x2=-a3x1 + a4x2 - b1x1x3;
%    x3=a5x3 + b2x1x2 + b3x1x4;
%    x4=a6x2 + a7x4 + b4x1x3;
   
%    x1 = 0.5 * [0.5 * [4/B] * [a1 * z1 + a2 * z4] - [2/B] * A * z2 * z3]

dx1 = a1*x1 + a2*x4 - x2*x3;
dx2 = -a3*x1 + a4*x2 + b1*x1*x3;
dx3 = a5*x3 + b2*x1*x2 + b3*x1*x4;
dx4 = a6*x2 + a7*x4 - b4*x1*x3;


x=x1;
y=x2;
z=x3;
q=x4;






%%%%%%%%%%%%%%%%%%%%%%%

sn_x = n2sn(x/Vmax);
sn_y = n2sn(y/Vmax);
sn_z = n2sn(z/Vmax);
sn_q = n2sn(q/Vmax);

c1 = n2sn(4*a1/B);
c2 = n2sn(4*a2/B);
c3 = n2sn(2*A/B);

x_1=mult(c1,sn_x);
x_2=mult(c2,sn_q);
x_3=mult(sn_y,sn_z);

dx=resta(suma(x_1,x_2),mult(c3,x_3));

%    x2 = 0.5 * [0.5 * [4/B] * [a4 * z2 - a3 * z1] + [2/B] * b1 * A * z1 * z3]
sn_x = n2sn(x/Vmax);
sn_y = n2sn(y/Vmax);
sn_z = n2sn(z/Vmax);
sn_q = n2sn(q/Vmax);

c1 = n2sn(4*a4/B);
c2 = n2sn(4*a3/B);
c3 = n2sn(2*A*b1/B);

y_1=mult(c1,sn_y);
y_2=mult(c2,sn_x);
y_3=mult(sn_x,sn_z);

dy=suma(resta(y_1,y_2),mult(c3,y_3));


%    x3 = 0.5 * [0.5 * [4/B] * [a5 * z3 + b1 * A * z1 * z2] + [2/B] * b3 * A * z1 * z4]
sn_x = n2sn(x/Vmax);
sn_y = n2sn(y/Vmax);
sn_z = n2sn(z/Vmax);
sn_q = n2sn(q/Vmax);

c1 = n2sn(4*a5/B);
c2 = n2sn(4*b2*A/B);
c3 = n2sn(2*A*b3/B);

z_1=mult(c1,sn_z);

z_2_1=mult(sn_x,sn_y);
z_2=mult(c2,z_2_1);

z_3=mult(sn_x,sn_q);

dz=suma(suma(z_1,z_2),mult(c3,z_3));


%    x4 = 0.5 * [0.5 * [4/B] * [a6 * z2 + a7 * z4] - [2/B] * b4 * A * z1 * z3]
sn_x = n2sn(x/Vmax);
sn_y = n2sn(y/Vmax);
sn_z = n2sn(z/Vmax);
sn_q = n2sn(q/Vmax);

c1 = n2sn(4*a6/B);
c2 = n2sn(4*a7/B);
c3 = n2sn(2*A*b4/B);

q_1=mult(c1,sn_y);
q_2=mult(c2,sn_q);
q_3=mult(sn_x,sn_z);

dq=resta(suma(q_1,q_2),mult(c3,q_3));



dx1=dx;
dx2=dy;
dx3=dz;
dx4=dq;


%%%%%%%%%%%%%%%%%%%%%%%

    
    
if(dx1>llindar) 
    x1 = x1+increment;
else
    x1 = x1-increment;
end

if(dx2>llindar) 
    x2 = x2+increment;
else
    x2 = x2-increment;
end

if(dx3>llindar) 
    x3 = x3+increment;
else
    x3 = x3-increment;
end

if(dx4>llindar) 
    x4 = x4+increment;
else
    x4 = x4-increment;
end


if(x1>Vmax)
   x1=Vmax; 
end

if(x2>Vmax)
   x2=Vmax; 
end

if(x3>Vmax)
   x3=Vmax; 
end

if(x4>Vmax)
   x4=Vmax; 
end

    
    
end

t= (c*iteracions1/NNN)+a;
(t)/iteracions1

x_o(t)=x1;
y_o(t)=x2;
z_o(t)=x3;
q_o(t)=x4;


end    
    
    
% Let's print to file every NN iterations....
fprintf(fileID, '%f %f  %f  %f  %f \n',n0+t, x,y,z,q);

%dibuixa=1;
if(dibuixa==0)
    figure(1);
    %hold off;
    plot(n0+(t1:t),x_o(t1:t)/Vmax,'r.');
    hold on
    %figure();
    plot(n0+(t1:t),y_o(t1:t)/Vmax,'g.');
    plot(n0+(t1:t),z_o(t1:t)/Vmax,'b.');
    plot(n0+(t1:t),q_o(t1:t)/Vmax,'y.');
end

t1=t;    
    
    
    
    
    
    
    
    
    
end



end   % cc

fclose(fileID);



plot(so_x,'r')
hold on;
plot(so_y,'g')
plot(so_z,'b')
plot(so_q,'y')

figure(3)
plot3(so_x,so_y,so_z)  %Figura 4_1 Article
title('Capell de bruixa')
xlabel('x1')
ylabel('x2')
zlabel('x3')
campos([0.8 2.1 1.8])

figure(4)
plot3(so_x,so_y,so_q)  %Figura 4_2 Article
title('Discòbol')
xlabel('x1')
ylabel('x2')
zlabel('x4')
campos([-1.3 -1.6 2.1])

figure(5)
plot3(so_y,so_z,so_q)  %Figura 4_3 Article
title('Cor trencat')
xlabel('x2')
ylabel('x3')
zlabel('x4')
campos([-0.02 0.08 -3.5])

figure(6)
plot3(so_z,so_x,so_q)  %Figura 4_4 Article
title('Peixet')
xlabel('x3')
ylabel('x1')
zlabel('x4')
campos([-2 -1 2])










%     function suma=suma(x,y)
%         suma=0.5*(x+y);
%     end
%    
%     function resta=resta(x,y)
%         resta=suma(x,-y);
%     end
% 
%     function mult=mult(x,y)
%         mult=x*y;
%     end
%    
%     function n=n2sn(x)
%         n=x;
%     end






    function n2sn=n2sn(x)
%            n2sn = (Vmax*0.5+x/2 > Vmax*rand());
            n2sn = (0.5+x/2 > rand());
    end

    function suma=suma(x,y)
               suma= (x&y) | ((~x)&y&(1/2>rand()) | ((~y)&x&(1/2>rand()) ) )  ;
    end
    
    function resta=resta(x,y)
                resta = suma(x,~y);
    end

    function mult=mult(x,y)
            mult = ((x&y) | ((~x)&(~y)));
    end












    end