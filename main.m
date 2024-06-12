%% time of arrival position triangulation
% This script solves a overdetermined time of arrival problem with a least
% square estimator in 2D, just like a GNSS systems does in 3D.
%
% input:    sat1, sat2, sat3 - sattelite positions
%           rho - measured sattelite pseudoranges
% settings: P0 - receiver position (initial guess)
%           b  - clock bias (initial guess
%           corr_norm - norm of solution vec correction at which the solver
%           stops.
%

% settings:
P0 = [1;0];         % initial position [x,y]
b = 1;              % initial clock bias [b]
corr_norm = 1e-5;   % solution correction norm

% input: positions and ranges
sat1 = [0; 15];
sat2 = [12; 3];
sat3 = [3; -6];
rho = [12; 7; 15];

iterations = 0;       % iteration counter
P = [P0; b];          % solution vector [x,y,b]
delta_P = [1; 1; 1];  % correction for solution vector

while norm(delta_P) > corr_norm
    % calc. range estimations
    rho_est = [norm(P0-sat1); norm(P0-sat2); norm(P0-sat3)];
    % build system matrix
    A = [(P0-sat1)'/rho_est(1), 1;
         (P0-sat2)'/rho_est(2), 1;
         (P0-sat3)'/rho_est(3), 1;];
    % calc. range corrections
    delta_rho = rho - rho_est - P(3);
    % solve for solution vec corrections
    delta_P = inv(A'*A)*A'*delta_rho;
    % update solution
    P = P +delta_P;
    P0 = P(1:2);
    % increment counter
    iterations = iterations + 1;
end

% draw satellites and solution
figure()
circle(sat1(1), sat1(2), rho(1))
circle(sat2(1), sat2(2), rho(2))
circle(sat3(1), sat3(2), rho(3))
scatter(P(1), P(2), 'red', 'X')

iterations
P
D = inv(A'*A);
GDOP = sqrt(sum(diag(D)))

function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp); hold on;
end