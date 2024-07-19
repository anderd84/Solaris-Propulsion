% aerospike nozzle utility script
% nozzle geometry defined by G. Angelino's approximation method
clear
close all
% defining geometry
% nozzle driving dimensions
% inputs
Pc = 300; %psia
Patm = 14; %psia

k = 1.4;





r_e = 1;
r_t = .9;
r_b = .4;

% nozzle derived requirements
M_e = 4;

% aerodynamic data

figure
hold on
for M_e = [2.4 3 3.5 4 4.5 5]
    % calculated parameters
    nu_e = PrandtlMeyer(k, M_e)
    A_star = pi*(r_e^2 - r_t^2)*cos(nu_e);
    
    
    M = linspace(1, M_e, 101);
    nu = PrandtlMeyer(k, M);
    mu = MachAngle(M);
    exitAreaRatio = (pi*r_e^2)/A_star;
    
    alpha = (nu_e - nu + mu);
    
    areaRatio = Isentropic1D(k, M, A_star) / A_star;
    eta = r_b/r_e;
    
    zeta = (1 - sqrt(1 - (areaRatio .* (1-eta^2) .* M .* sin(alpha) / exitAreaRatio)));
    zeta = zeta ./ sin(alpha);
    
    
    x = zeta .* cos(alpha);
    y = zeta .* sin(alpha);
    
    % yyaxis right
    % plot(linspace(0,max(x),101), zeta)
    % plot(linspace(0,max(x),101), alpha)
    % yyaxis left
    plot(x,y)
    % plot(x(y < r_e-r_b & x>0),r_e - y(y < r_e-r_b & x>0))
    % plot(x(y < r_e-r_b & x>0),y(y < r_e-r_b & x>0) - r_e)
end
% set(gca, 'YDir', 'reverse')
axis equal
grid on

function A = Isentropic1D(k, M, A_star)
    f1 = (k-1)/2;
    f2 = (k+1)/(2*(k-1));
    A = A_star*(1./M).*((1+f1*M.^2)/(1+f1)).^f2;
end

function nu = PrandtlMeyer(k, M)
    f1 = (k+1)./(k-1);
    f2 = (M.^2 - 1);
    nu = sqrt(f1).*atan(sqrt(f2./f1)) - atan(sqrt(f2));
end

function mu = MachAngle(M)
    mu = asin(1./M);
end