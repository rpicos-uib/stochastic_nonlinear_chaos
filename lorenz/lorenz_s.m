%function [x_o, y_o, z_o]=shimizu_morioka_s()
function lorenz_s(dibuixa2)


if (nargin==0)
  dibuixa=1;
else
  dibuixa=dibuixa2;
end

llindar = 0.5
%llindar1 = 0.508; %(si no, estam integrant "0", i el resultat és zero...)

% iter2 = 2^8, Vmax=2^22
% 1000 iters en 1.6414s
% 10000 iters en 15.3314s

fitxer='lorenz.txt';


iteracions1 = 1e8;

% Aquest funciona.... lento, però funciona
iteracions2 = 2^24;
Vmax = 2 ^22;
NNN=2000
x = Vmax *0.1;
y = Vmax *0.1;
z = Vmax *0.1;
n0=1;

% Condicions inicials per continuar:
% Llegim el fitxer, i continuam allà on erem....

if(1<1)
    fileID=fopen(fitxer,"a");
    comanda=['tail -n 1 ' fitxer]
    [status result] = system(comanda);
    if isstrprop(result(end), 'cntrl'), result(end) = []; end
    n=sscanf(result,'%f %f %f %f');
    n0=n(1)
    x=n(2)
    y=n(3)
    z=n(4)
else
    comanda=['rm ' fitxer]
    [status result] = system(comanda);
    fileID=fopen(fitxer,"w");    
end

% Test
iteracions2 = 2^12;
Vmax = 2 ^28;
NNN=iteracions1/100;



%   Normalized with t to have everything <1

%    dx = y/8 - 1/16;
%    dy = x/2 -1/4 + mu/16 -mu*y/8 + 3*z/8 -3/4*x*z;
%    dz = -alpha*z/8 + alpha/16 + 1/3*x*x + 1/12 - x/3;

sigma=10*Vmax;
beta=8/3*Vmax;
rho=28*Vmax;
A=50*Vmax;
D=max([ A rho beta sigma]);
E=4*Vmax;
%D=1;


x_o=zeros(iteracions1,1);
y_o=x_o;
z_o=x_o;

t=0;

sortida1=0;
sortida2=0;

inc=1;

t1=1;

for c= 0:(NNN-1)


for a = 1:iteracions1/NNN
    
for b = 1:iteracions2

    

%    dx = y/8 - 1/16; = 1/2(y/4 - 1/8)
snx_y = n2sn(y);
snx_x = n2sn(x);

sigmap=n2sn(2*sigma/D/E);
rhop=n2sn(rho/D);
c2 = n2sn(A/D);
c3 = n2sn(2*A/D/E);
c4=n2sn(2/D/E);
c5 = n2sn(2*beta/D/E);
c6=n2sn(4/E);

dx = mult(sigmap, resta(snx_y,snx_x));


sny_z = n2sn(z);
sny_x = n2sn(x);
sny_y = n2sn(y);

dy1 = resta(rhop, mult(c2,sny_z));
dy = resta(mult(mult(c6,sny_x),dy1), mult(c4,sny_y));

snz_z = n2sn(z);
snz_x = n2sn(x);
snz_y = n2sn(y);


dz = resta(mult(c3,mult(snz_x,snz_y)), mult(c5,snz_z));



if(dx>llindar) 
    x = x+1;
else
    x = x-1;
end

if(dy>llindar)
    y = y+1;
else
    y = y-1;
end

if(dz>llindar)
    z = z+1;
else
    z = z-1;
end



if(x>Vmax)
   x=Vmax; 
end

if(y>Vmax)
   y=Vmax; 
end

if(z>Vmax)
   z=Vmax; 
end



end

t= (c*iteracions1/NNN)+a;
(t)/iteracions1

x_o(t)=x;
y_o(t)=y;
z_o(t)=z;



end

% Let's print to file every NN iterations....
fprintf(fileID, '%d %d  %d  %d \n',n0+t, x,y,z);

%dibuixa=1;
if(dibuixa==0)
    figure(1);
    %hold off;
    plot(n0+(t1:t),x_o(t1:t)/Vmax,'r.');
    hold on
    %figure();
    plot(n0+(t1:t),y_o(t1:t)/Vmax,'g.');
    plot(n0+(t1:t),z_o(t1:t)/Vmax,'b.');
end

t1=t;

end

fclose(fileID);

figure(4);
plot3(x_o/Vmax,y_o/Vmax,z_o/Vmax);


    function n2sn=n2sn(x)
            n2sn = (Vmax*0.5+x/2 > Vmax*rand());
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