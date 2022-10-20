function plot_SNR()

figure(40)
close(40);
% We need first to downsample the exact solution, so we read another one..
Nb=19;
[t, x, y, z] = read_from_file(Nb,'_3');

tf=max(t);
tf=100;

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

length(xx)
figure(40)


Nb=16;
plot_range(Nb,tf,t2,xx,1,5);

Nb=17;
plot_range(Nb,tf,t2,xx,1,1);

Nb=18;
plot_range(Nb,tf,t2,xx,1,5);

Nb=19;
plot_range(Nb,tf,t2,xx,0,3);

Nb=20;
plot_range(Nb,tf,t2,xx,1,2);

Nb=22;
plot_range(Nb,tf,t2,xx,1,1);

Nb=24;
plot_range(Nb,tf,t2,xx,1,1);


i=16:0.1:24;
SNRdB=4*log10(2.^(i/2-4));
SNRdB2=(2*(i-9))*log10(2);


plot(i,SNRdB,'b');
plot(i,SNRdB2,'k');
















xlabel('# of bits');
ylabel('Energy Signal/Noise (a.u.)');


    function plot_range(Nb,tf,t2,xx,iin,ifin)
        counter=0;
        suma=0;
        for i=iin:ifin
            [t, x, y, z] = read_from_file(Nb,['_' num2str(i)]);
            x=x(t<tf);
            t=t(t<tf);
            t1=length(t);
            ratio=(t2/t1);
            ex=xx(1:ratio:end)';
            SNRdB=get_SNR(ex,x(t<tf), Nb);
            suma=suma+SNRdB;
            plot(Nb, SNRdB,'ro');
            hold on;
            counter=counter+1;
        end
        plot(Nb,suma/counter,'xg');
    end


    function SNRdB=get_SNR(ex,x, Nb)
        noise = (abs(ex)-abs(x));
        SNR = mean(abs(ex))./mean(abs(noise));  % DEFINITION 1
        % However, I think (on 7/Oct/2022), that it makes more sense using
        % the energy (as by Parceval's theorem):
        % (notice that we're using the same bin width, so no problems from
        % this side).
        SNR = sum(abs(ex))./sum(abs(noise));  % DEFINITION 2
        
        % It seems that the result is the same... (as expected?)
        
        SNRdB=20*log10(SNR)
        
    end

end

