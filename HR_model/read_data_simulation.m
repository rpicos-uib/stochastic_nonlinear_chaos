function [ex,x]=read_data_simulation()
close all

% This functions plots the Hindmarsh-Rose model as calculated by
% the function hr_model_normal_1_s.m
%
% We need first to downsample the exact solution, so we read another one..
Nb=19;
[t, x, y, z] = read_from_file(Nb,'_3');


tf=max(t);

[xt,xx,xy,xz]=hr_model_normal_1(tf,0.0);

% Just let's get the same length for the 'exact' solution....
% so we downsample it....
t1=length(t);
t2=length(xt);
ratio=(t2/t1)
et=xt(1:ratio:end)';
ex=xx(1:ratio:end)';
ey=xy(1:ratio:end)';
ez=xz(1:ratio:end)';

displayname='Exact';

figure(20)
subplot(4,1,1)
plot_read_data(et,ex,ey,ez,displayname);
figure(20)
subplot(4,1,1)
xlim([0 100]);

[xt,xx,xy,xz]=hr_model_normal_1(10*tf);
plot_read_map(xx,xy,xz,displayname,2,2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%
% Nb=24;
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

Nb=20;

[t, x, y, z] = read_from_file(Nb,'_1');
displayname=[num2str(Nb) ' bits'];

figure(20)
subplot(4,1,2)
plot_read_data(t,x,y,z,displayname);
figure(20)
subplot(4,1,2)
xlim([0 100]);

t=t(1:1000);
x=x(1:1000);
y=y(1:1000);
z=z(1:1000);

%let's check the SNR:
get_SNR(ex,x, Nb);

plot_read_map(x,y,z,displayname,2,2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nb=19;

[t, x, y, z] = read_from_file(Nb,'');
displayname=[num2str(Nb) ' bits'];
figure(20)
subplot(4,1,3)
plot_read_data(t,x,y,z,displayname);
figure(20)
subplot(4,1,3)
xlim([0 100]);


t=t(1:min(end,1000));
x=x(1:min(end,1000));
y=y(1:min(end,1000));
z=z(1:min(end,1000));


%let's check the SNR:
get_SNR(ex,x, Nb);

plot_read_map(x,y,z,displayname,2,2,3);
%plot_read_map(x(x<0.15),y(x<0.15),z(x<0.15),displayname,2,2,2);
displayname=[num2str(22) ' bits'];
plot_read_map(x(y<0.8 & x<0.4),y(y<0.8 & x<0.4),z(y<0.8 & x<0.4),displayname,2,2,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Nb=16;

[t, x, y, z] = read_from_file(Nb,'_2');
displayname=[num2str(Nb) ' bits'];
figure(20)
subplot(4,1,4)
plot_read_data(t,x,y,z,displayname);
figure(20)
subplot(4,1,4)
xlim([0 100]);

t=t(1:1000);
x=x(1:1000);
y=y(1:1000);
z=z(1:1000);
%let's check the SNR:
get_SNR(ex,x, Nb);

plot_read_map(x,y,z,displayname,2,2,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

Nb=12;

% [t, x, y, z] = read_from_file(Nb);
% displayname=[num2str(Nb) ' bits'];
% plot_read_data(t,x,y,z,displayname);
% 
% t=t(1:1000);
% x=x(1:1000);
% y=y(1:1000);
% z=z(1:1000);
% %let's check the SNR:
% get_SNR(ex,x, Nb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


    function plot_read_map(x,y,z,displayname,n,m,p)
        figure(50)
        subplot(n,m,p);
%        plot3(ex,ey,ez,'DisplayName',displayname);
        plot(x,y,'b.','DisplayName',displayname);
        hold on
        title(displayname);
        xlabel('X (a.u.)');
        ylabel('Y (a.u.)');
        xlim([0.1 0.8]);
        ylim([0.2 1.0]);
        
        plot_derivative(x,y,z,displayname,m,n,p);
        
    end

    function get_SNR(ex,x, Nb)
        noise = (abs(ex)-abs(x));
        SNR = mean(abs(ex))./mean(abs(noise));  % DEFINITION 1
        % However, I think (on 7/Oct/2022), that it makes more sense using
        % the energy (as by Parceval's theorem):
        % (notice that we're using the same bin width, so no problems from
        % this side).
        SNR = sum(abs(ex))./sum(abs(noise));  % DEFINITION 2

        % It seems that the result is the same... (as expected?)
        
        SNRdB=20*log10(SNR)
        
        figure(40)
        plot(Nb, SNRdB,'ro');
        hold on;
        xlabel('# of bits');
        ylabel('Energy Signal/Noise (a.u.)');
        
    end




    function plot_derivative(X,Y,Z,displayname,m,n,p)
        
        x0=X*0;
        y0=x0;
        z0=y0;
        i0=25;
        
        for i=i0+1:length(x0)-i0
            x0(i)=sum(X(i-i0:i+i0))/(2*i0+1);
            y0(i)=sum(Y(i-i0:i+i0))/(2*i0+1);
            z0(i)=sum(Z(i-i0:i+i0))/(2*i0+1);
        end
        
        x0=x0(i0+1:end-i0);
        y0=y0(i0+1:end-i0);
        z0=z0(i0+1:end-i0);
        
        dx=diff(x0);
        dy=diff(y0);
        dz=diff(z0);
        
        figure(31);
        subplot(n,m,p);
        plot3(dx,dy,dz);
        title(displayname);
        xlabel('dx');
        ylabel('dy');
        zlabel('dz');
        
    end






end