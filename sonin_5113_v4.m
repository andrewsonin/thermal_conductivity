% Andrew Sonin, 5113. Variant 4
% Two-layer six-point difference scheme for solving the heat equation.

% u'_t = u''_{xx} + e^t * sin(x), x \in [-\pi/2, \pi/2], t \in [0, T]
% u(t=0) = sin(x), u(x=-\pi/2) = -cosh(t) = -u(x=\pi/2)
% exact solution: u = cosh(t) * sin(x)

x_0 = -pi/2;
x_N = -x_0;
t_0 = 0;
t_N = input('Input the edge time: ');
ksi = input('Input the method parameter (ksi): ');
P = input('Input partition parameter: ');

h = (x_N - x_0) / P;
tau = (t_N - t_0) / P;

fprintf('Parabolic Courant number = %f\n', tau / h^2);

U = zeros(P + 1, P + 1); % Numerical solution
f = [];

for i = 0 : P
    U(1, i + 1) = sin(x_0 + h * i);
end

A = zeros(P + 1, P + 1);
A(1, 1) = 1; A(P + 1, P + 1) = A(1, 1);

f(1) = -cosh(tau);
f(P + 1) = -f(1);

boundary = -ksi / h^2;
central = 1 / tau - 2 * boundary;
c_1 = (1 - ksi) / h^2;
c_2 = 1 / tau - 2 * c_1;
for i = 2 : P  % Filling the first time layer
    A(i, i) = central;
    A(i, i - 1) = boundary;
    A(i, i + 1) = boundary;
    f(i) = c_2 * U(1, i) + c_1 * (U(1, i + 1) + U(1, i - 1)) + sin(x_0 + (i - 1) * h);
end

for i = 2 : P + 1  % Time layers filling
    U(i, :) = A \ f';
    for j = 2 : P  % Spatial layers filling
        f(j) = c_2 * U(i, j) + c_1 * (U(i, j + 1) + U(i, j - 1)) + sin(x_0 + (j - 1) * h) * exp(tau * (i - 1));
    end
    f(1) = -cosh(tau * i);
    f(P + 1) = -f(1);
end

[X, Y] = meshgrid(t_0 : tau : t_N, x_0 : h : x_N); % create the grid
ES = sin(Y).*(cosh(X)); %exact solution

% plot3(X, Y, abs(ES - U'), 'r'); % Error plot
plot3(X, Y, ES, 'y')
hold on;
plot3(X, Y, U', 'b')
