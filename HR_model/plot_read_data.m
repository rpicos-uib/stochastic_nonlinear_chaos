    function plot_read_data(t,x,y,z,displayname)
        figure(20)
        plot(t,x,'DisplayName',displayname)
        hold on
        % plot(t, y);
        % plot(t, z);
        legend
        xlabel('Time (s)');
%        ylabel('Variable (normalized)');
        ylabel('x (norm.)');
        
        figure(21)
        fx=fft(x);
        % fy=fft(y);
        % fz=fft(z);
        
        
        loglog(abs(fx),'.','DisplayName',displayname);
        hold on;
        legend
        xlabel('Freq. (a.u.)');
        ylabel('Magnitude (a.u.)');
        
        
        figure(22)
        plot3(x,y,z,'DisplayName',displayname);
        hold on;
        legend
        xlabel('X (a.u.)');
        ylabel('Y (a.u.)');
        zlabel('Z (a.u.)');
        
        
    end