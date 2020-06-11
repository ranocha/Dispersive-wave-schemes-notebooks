%%% Solitary wave to the Holm-Hone equation using the Petviashvili
%%% method

function a = fw_traveling_wave(c, l, N)

% Input  : - c,  speed (must be > 1)
%          - l,  half-length of the domain
%          - N,  number of Fourier modes used in the computation
% Output : fw_traveling_wave_init.txt file with the solution (x, u(x))
% Example: fw_traveling_wave(1.2)
    
    if nargin < 3
        N = 512;  % number of Fourier modes
    end
    if nargin < 2
        l = 20.0; % half-length of the domain
    end
    if nargin < 1
        c = 1.2;  % speed
    end

    dx = 2*l/N;             % distance between two points in the real space
    x = (-N/2:N/2-1)'*dx;   % real space discretization

    % vector of wavenumbers
    k = [0:N/2-1 0 1-N/2:-1]'*pi/l;

    j = (N/4+2:N/4*3);      % antialising treatment
    k(j) = 0;               % the frequencies we sacrify
    
    L = c*(1 + k.^2) - 1; % linear operator
    Linv = 1./L;          % inverse linear operator

    u0 = exp(-x.^2);
    

    tol = 1e-15;        % iterations stopping criterium
    gam = 1.5;          % a parameter in the Petviashvili method
    
    err = inf;
    iter = 0;
    while (err > tol)
        iter = iter + 1;
        u0_hat = fft(u0);
        Nl_hat = (1 + k.^2) .* fft(0.5 * u0.^2);
        Nl_hat(j) = 0;
        Lu0_hat = L .* u0_hat;
        nom = sum(conj(u0_hat) .* Nl_hat);
        denom = sum(conj(u0_hat) .* Lu0_hat);
        S = (denom/nom)^gam;
        u1_hat = S * Linv .* Nl_hat;
        u1 = real(ifft(u1_hat));
        err = norm(u1-u0, inf);
        fprintf('iteration %5d; error = %.2e\n', iter, err);
        u0 = u1;
    end
    plot(x,u0)
    
    
    io = fopen('fw_traveling_wave_init.txt','w');
    fprintf(io, '# x u0 \n');
    fprintf(io, '# Generated using the following parameters \n');
    fprintf(io, '# c  = %20.15e \n', c);
    fprintf(io, '# l  = %20.15e \n', l);
    fprintf(io, '# N  = %d \n', N);
    for i=1:N
        fprintf(io,'%20.15f %20.15e\n', x(i), u0(i));
    end
    fprintf(io,'%20.15f %20.15e\n', -x(1), u0(1));
  
    fprintf('min(u0) = %.2e\n', min(u0));
    fprintf('max(u0) = %.2e\n', max(u0));
    fprintf('u0( 1 ) = %.2e\n', u0(1));
    fprintf('u0(end) = %.2e\n', u0(end));
    a = max(u0);
end