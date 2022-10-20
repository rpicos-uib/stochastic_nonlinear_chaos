function read_data_simulation_vs_SC

% Just closes the windows we are using here....
try
close(101)
close(100)
end

Nb=19;

tf=100;

noise_level=0.02;
noise_level=2^(-Nb/2+4);
displayname=['Noisy D.E. (noise =' num2str(100*noise_level) '%)'];
delta_step=5;

[tt, xx1, yy1, zz1]=hr_model_normal_1(tf, noise_level);
tt=tt(1:delta_step:end);
xx1=xx1(1:delta_step:end);
yy1=yy1(1:delta_step:end);
zz1=zz1(1:delta_step:end);


[tt, xx2, yy2, zz2]=hr_model_normal_1(tf, noise_level);
tt=tt(1:delta_step:end);
xx2=xx2(1:delta_step:end);
yy2=yy2(1:delta_step:end);
zz2=zz2(1:delta_step:end);

[tt, xx3, yy3, zz3]=hr_model_normal_1(tf, noise_level);
tt=tt(1:delta_step:end);
xx3=xx3(1:delta_step:end);
yy3=yy3(1:delta_step:end);
zz3=zz3(1:delta_step:end);


figure(101)
subplot(2,1,1)
plot(tt(tt<tf),xx1(tt<tf),'DisplayName','simul_1')
title(displayname);
hold on

subplot(2,1,1)
plot(tt(tt<tf),xx2(tt<tf),'DisplayName','simul_2')
title(displayname);

subplot(2,1,1)
plot(tt(tt<tf),xx3(tt<tf),'DisplayName','simul_3')
title(displayname);
legend

ll=length(tt(tt<tf))



% Let's get the parameters for the Fourier:
% Sampling frequency
Fs = max(tt(tt<tf))/ll
% Sampling period
T = 1/Fs;
% Length of signal
L = ll;
% Time vector
t = tt;

% Let's compare the three FFTs:
fft1 = fft(xx1(tt<tf));
fft2 = fft(xx2(tt<tf));
fft3 = fft(xx3(tt<tf));

freq=Fs*(0:(ll/2))/ll;
% Compute the two-sided spectrum P2. 
%Then compute the single-sided spectrum P1 based on P2 
% and the even-valued signal length L.
% see https://www.mathworks.com/help/matlab/ref/fft.html
P2 = abs(fft1/L);
fft1 = P2(1:L/2+1);
fft1(2:end-1) = 2*fft1(2:end-1);
P2 = abs(fft2/L);
fft2 = P2(1:L/2+1);
fft2(2:end-1) = 2*fft2(2:end-1);
P2 = abs(fft3/L);
fft3 = P2(1:L/2+1);
fft3(2:end-1) = 2*fft3(2:end-1);

delta_f=1;
figure(100)
semilogx(freq, (fft1(1:delta_f:end)),'.','DisplayName',displayname);
hold on;
loglog(freq, (fft2(1:delta_f:end)),'.','DisplayName',displayname);
loglog(freq, (fft3(1:delta_f:end)),'.','DisplayName',displayname);
xlabel('Freq. (Hz)');
ylabel('Magnitude (a.u.)');


%%%%%%%%%%%%%%%%%%%%%%%
% Now, the SC data:
%%%%%%%%%%%%%%%%%%%%%%%

Nb=19;

fes_coses(Nb,'_1');


%%%%%%%%%%%%%%%%%%%%%%%
% Now, the SC data:
%%%%%%%%%%%%%%%%%%%%%%%

Nb=19;

fes_coses(Nb,'_2');

%%%%%%%%%%%%%%%%%%%%%%%
% Now, the SC data:
%%%%%%%%%%%%%%%%%%%%%%%

Nb=19;

fes_coses(Nb,'_3');


    function fes_coses(Nb, cadena)
[ts, xs, ys, zs] = read_from_file(Nb, cadena);
displayname=[num2str(Nb) ' bits'];

figure(101)
subplot(2,1,2)
plot(ts,xs,'DisplayName',displayname);
hold on
title(displayname);
xlim([0 tf]);
xlabel('Time (s)');

% Let's FFT....
figure(100)
fft_sc=fft(xs(ts<tf));
ll_sc=length(xs(ts<tf));
FS_sc=max(ts)/ll_sc
freq_sc=FS_sc*(0:(ll_sc/2))/ll_sc;

P2 = abs(fft_sc/ll_sc);
fft_sc = P2(1:ll_sc/2+1);
fft_sc(2:end-1) = 2*fft_sc(2:end-1);
loglog(freq_sc, fft_sc,'.','DisplayName',displayname);
legend
end
end