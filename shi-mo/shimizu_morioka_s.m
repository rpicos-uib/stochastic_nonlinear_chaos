%function [x_o, y_o, z_o]=shimizu_morioka_s()
function shimizu_morioka_s()

llindar = 0.5
%llindar1 = 0.508; %(si no, estam integrant "0", i el resultat és zero...)

% iter2 = 2^8, Vmax=2^22
% 1000 iters en 1.6414s
% 10000 iters en 15.3314s

fitxer='shimizu_morioka.txt';

fileID=fopen(fitxer,"a");

iteracions1 = 1e6;

% Aquest funciona.... lento, però funciona
iteracions2 = 2^12;
Vmax = 2 ^22;
NNN=100
x = Vmax *0.51;
y = Vmax *0.51;
z = Vmax *0.51;

% Condicions inicials per continuar:
% Llegim el fitxer, i continuam allà on erem....
comanda=['tail -n 1 ' fitxer]
[status result] = system(comanda);
if isstrprop(result(end), 'cntrl'), result(end) = []; end
n=sscanf(result,'%f %f %f %f');
n0=n(1)
x=n(2)
y=n(3)
z=n(4)

% Test
iteracions2 = 2^12;
Vmax = 2 ^22;
NNN=iteracions1/100;



%   Normalized with t to have everything <1

%    dx = y/8 - 1/16;
%    dy = x/2 -1/4 + mu/16 -mu*y/8 + 3*z/8 -3/4*x*z;
%    dz = -alpha*z/8 + alpha/16 + 1/3*x*x + 1/12 - x/3;

alpha=3/8;
mu=0.81;

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
snx_y = n2sn(y/Vmax);
cx1 = n2sn(1/16);
cx2 = n2sn(1/32);
dx = resta(mult(snx_y,cx1), cx2 );


%    dy = x/2 -1/4 + mu/16 + 3*z/8 -3/4*x*z -mu*y/8;
%       = 1/2(x-1/2) + 1/2(mu/8 + z*3/4) - 1/2( 3/2*x*z + mu*y/4)
sny_x = n2sn(x/Vmax);
sny_z = n2sn(z/Vmax);
sny_x2 = n2sn(x/Vmax);
sny_y = n2sn(y/Vmax);
sny_z2 = n2sn(z/Vmax);
cy1=n2sn(1/2);
cy2=n2sn(mu/8);
cy3=n2sn(3/4);
cy4=n2sn(3/4);
cy5=n2sn(mu/8);

dy1= suma(resta(sny_x, cy1), suma(cy2, mult(cy3, sny_z)));
dy2 = suma(mult(cy4,mult(sny_x2,sny_z2)), mult(cy5,sny_y));
dy = resta(dy1, dy2);

%    dz = -x/3 -alpha*z/8 + alpha/16 + 1/3*x*x + 1/12;
cz2 = n2sn(alpha/4);
cz3 = n2sn(alpha/8);
cz4 = n2sn(2/3);
cz5 = n2sn(1/24);
cz6 = n2sn(2/3);

snz_x = n2sn(x/Vmax);
snz_z = n2sn(z/Vmax);
dz1 = suma(mult(cz4,snz_x), mult(cz2, snz_z));
snz_x2 = n2sn(x/Vmax);
snz_x3 = n2sn(x/Vmax);
dz2 = suma(cz3, mult(cz6, mult(snz_x3, snz_x2) ) );
dz3 = resta(dz2, dz1);
dz = suma(dz3, cz5);


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

dibuixa=0;
if(dibuixa==1)
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