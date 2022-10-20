    function [t, x, y, z] = read_from_file(varargin)
%https://www.mathworks.com/help/matlab/matlab_prog/support-variable-number-of-inputs.html
    
    cadena='_1';
    Nb=varargin{1};
    if(nargin>1)
        cadena=varargin{2}
    end


    % Let's open a file so we can interrupt the process...
        fitxer=['hindmarsh-rose-stochastic_' num2str(Nb) 'bits' cadena '.txt'];
        fileID=fopen(fitxer,"r");
        %comanda=['tail -n 1 ' fitxer]
        %[status result] = system(comanda);
        %if isstrprop(result(end), 'cntrl'), result(end) = []; end
        % n=sscanf(result,'%f %f %f %f');
        M=csvread(fitxer);
        
        
        t=M(:,1);
        x=M(:,2);
        y=M(:,3);
        z=M(:,4);
        
        
    end
