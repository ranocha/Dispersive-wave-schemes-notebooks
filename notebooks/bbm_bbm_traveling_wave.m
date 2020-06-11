%%% Solitary wave to the BBM-BBM system using the Petviashvili
%%% method
%%% ----------------------------------------------------------
%%% Author: Dimitrios Mitsotakis
%%% April 6, 2020
%%% Modified by Hendrik Ranocha

function a = bbm_bbm_traveling_wave(c, l, N)

% Input  : - c, speed (must be > 1)
%          - l, half-length of the domain
%          - N, number of Fourier modes used in the computation
% Output : bbm_bbm_traveling_wave_init.txt file with the solution (x,eta(x),u(x))
% Example: bbm_bbm_traveling_wave(1.2)
    
    if nargin < 3
        N = 512;  % number of Fourier modes
    end
    if nargin < 2
        l = 20.0; % half-length of the domain
    end
    if nargin < 1
        c = 1.2;  % speed (must be > 1)
    end

    g = 1.0;                % gravity acceleration
    d = 1.0;                % undisturbed water depth

    dx = 2*l/N;             % distance between two points in the real space
    x = (-N/2:N/2-1)'*dx;   % real space discretization

    % vector of wavenumbers
    k = [0:N/2-1 0 1-N/2:-1]'*pi/l;

    j = (N/4+2:N/4*3);      % antialising treatment
    k(j) = 0;               % the frequencies we sacrify
    
%     L = c*(1 + 1/6*(d*k).^2);   % linear operator
    L = c*(1 + 1*(d*k).^2);   % linear operator

    
    row1 = [L, -d*ones(size(L))];
    row2 = [-g*ones(size(L)), L];
    
    size(row1)
    %%% Initial guess specification:
    amp = c^2/g - d; kk = sqrt(3*amp/(d+amp))/d;
    x0=0;
    eta0 = amp*sech(0.5*kk*(x-x0)).^2;
    u0 = c*eta0./(d+eta0);


    tol = 1e-14;        % iterations stopping criterium
    gam = 2.0;          % a parameter in the Petviashvili method
    
    A = zeros(2,2);
    
    eta1_hat = zeros(size(eta0));
    u1_hat = zeros(size(u0));
    
    err = inf;
    iter = 0;
    while (err > tol)
        iter = iter + 1;
        u0_hat = fft(u0);
        eta0_hat = fft(eta0);
        Nl1 = eta0 .* u0; 
        Nl1_hat = fft(Nl1); Nl1_hat(j) = 0;
        Nl2 = 0.5 * u0.^2;
        Nl2_hat = fft(Nl2); Nl2_hat(j) = 0;
        
        denom1 = sum(conj(eta0_hat) .* Nl1_hat);
        denom2 = sum(conj(u0_hat) .* Nl2_hat);
        denom = denom1 + denom2;

        Leta0_hat = L .* eta0_hat;
        Lu0_hat = L .* u0_hat;
        nom1 = sum(conj(eta0_hat) .* Leta0_hat) - d * sum(conj(u0_hat) .* eta0_hat);
        nom2 = sum(conj(u0_hat) .* Lu0_hat) - g * sum(conj(eta0_hat) .* u0_hat);
        nom = nom1 + nom2;
        
        S = (nom/denom)^gam;       
        
        for ik = 1: length(k)
            A(1,:) = row1(ik,:);
            A(2,:) = row2(ik,:);
            
            w1_hat = S * (A \ [Nl1_hat(ik); Nl2_hat(ik)]);
            eta1_hat(ik) = w1_hat(1);
            u1_hat(ik) = w1_hat(2);
        end
        
        eta1 = real(ifft(eta1_hat));
        u1 = real(ifft(u1_hat));

        R1 = abs(real(ifft(-L.*eta1_hat+d*u1_hat+fft(eta1.*u1))));
        R2 = abs(real(ifft(-L.*u1_hat+g*eta1_hat+0.5*fft(u1.^2))));

%        err = max(norm(u1-u0, inf),norm(eta1-eta0, inf))
        err = max(norm(R1, inf), norm(R2, inf));
        fprintf('iteration %5d; error = %.2e\n', iter, err);

        u0 = u1;
        eta0 = eta1;
    end
    plot(x, eta0, '-', 'DisplayName', 'eta0')
    hold on
    plot(x, u0, '--', 'DisplayName', 'u0')
    legend()
    hold off
    
    io = fopen('bbm_bbm_traveling_wave_init.txt','w');
    fprintf(io, '# x eta0 u0 \n');
    fprintf(io, '# Generated using the following parameters \n');
    fprintf(io, '# c = %20.15e \n', c);
    fprintf(io, '# l = %20.15e \n', l);
    fprintf(io, '# N = %d \n', N);
    for i=1:N
        fprintf(io,'%20.15f %20.15e %20.15e\n', x(i), eta0(i), u0(i));
    end
    
    fprintf('min(eta0) = %.2e\n', min(eta0));
    fprintf('max(eta0) = %.2e\n', max(eta0));
    fprintf('min( u0 ) = %.2e\n', min( u0 ));
    fprintf('max( u0 ) = %.2e\n', max( u0 ));
    fprintf('eta0( 1 ) = %.2e\n', eta0(1));
    fprintf('eta0(end) = %.2e\n', eta0(end));
    fprintf(' u0( 1 )  = %.2e\n', u0(1));
    fprintf(' u0(end)  = %.2e\n', u0(end));
    a = max(eta0);
end
