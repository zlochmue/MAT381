% Task 1: Solve Example 2
% The algorithm is from the paper:
% A third order Newton-type method to solve systems
% of nonlinear equations, by Darvishi and Barati, (2007).
% Settles on a minimum of [2.4177, 2.8768, 3.5508]
% iteration 1: [1.3275, 2.9197, 1.5280]
% iteration 2: [2.8468, 4.0922, 3.0468]
% iteration 8: [2.4177, 2.8768, 3.5508]

clear all;
initial = [-1; 1; -1];
xvect = initial
for iter = 1:10
 x = xvect(1);
 y = xvect(2);
 z = xvect(3);
 Fvect_x = [(4*(x-1)^3) + (6*(x-y)^5) + (6*(x-z)^5); ...
 (4*(y-3)^3) - (6*(x-y)^5) + (6*(y-z)^5); ...
 (4*(z-5)^3) - (6*(y-z)^5) - (6*(x-z)^5)];
 F_prime_x = [ ...
 (12*(x-1)^2) + (30*(x-y)^4) + (30*(x-z)^4), (-30*(x-y)^4), (-30*(x-z)^4); ...
 (-30*(x-y)^4), (12*(y-3)^2) + (30*(x-y)^4) + (30*(y-z)^4), (-30*(y-z)^4); ...
 (-30*(x-z)^4), (-30*(y-z)^4), (12*(z-5)^2) + (30*(y-z)^4) + (30*(x-z)^4)];
 yvect = xvect - inv(F_prime_x)*Fvect_x ;
 x2 = yvect(1);
 y2 = yvect(2);
 z2 = yvect(3);
 Fvect_y = [(4*(x2-1)^3) + (6*(x2-y2)^5) + (6*(x2-z2)^5); ...
 (4*(y2-3)^3) - (6*(x2-y2)^5) - (6*(y2-z2)^5); ...
 (4*(z2-5)^3) - (6*(y2-z2)^5) - (6*(x2-z2)^5)];
 xvect = xvect - inv(F_prime_x)*( Fvect_x + Fvect_y );
 it = num2str(iter);
 message = ['at iteration ', it, ', we have [x, y, z] ='];
 disp(message)
 disp(round(xvect,10))
end

% Task 2: Solve Example 4
% The algorithm is from the paper:
% A third order Newton-type method to solve systems
% of nonlinear equations, by Darvishi and Barati, (2007).
%
% Solve for x1 and x2 to minimize
% f(x1, x2) =sum_{j=1}^{N} [y_j - exp(x1 t_j) - exp(x2 u_j)]^2
%
% After 10 iterations : [1.7893, 4.1160]
% Greater than 25 iterations : [1.7321, 4.1231] consistently
clear all;
%
% make data
%

 N = 200;
 bvect = zeros(N,1);
 seed = 101;
 rng(seed); % set seed for random number generator
 Amat = rand(N,2);
 for j = 1:N
 t = Amat(j,1);
 u = Amat(j,2);
 bvect(j) = exp(sqrt(3)*t)+exp(sqrt(17)*u);
 end
 scatter3(Amat(:,1),Amat(:,2),bvect,10)
 
%
% Solve the nonlinear system of equations
%
 initial = [0.75; 0.25];
 NumIter = 10; % number of iterations
xvect = initial
for iter = 1:NumIter
 x1 = xvect(1);
 x2 = xvect(2);
 sum1 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);
 u = Amat(j,2);
 sum1 = sum1 + ((-exp(x1*t) - exp(x2*u) + y)*(-exp(x1*t)*t));
 end
 sum2 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);
 u = Amat(j,2);
 sum2 = sum2 + ((-exp(x1*t) - exp(x2*u) + y)*(-exp(x2*u)*u));
 end
 Fvect_x = [sum1; sum2];
 sum3 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);

u = Amat(j,2);
 sum3 = sum3 + ((2*exp(x1*t) - exp(x2*u) + y)*((t*t)*exp(x1*t)));
 end
 sum4 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);
 u = Amat(j,2);
 sum4 = sum4 + ((t * u * exp((t*x1) + (u*x2))));
 end
 sum5 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);
 u = Amat(j,2);
 sum5 = sum5 + ((2*exp(x2*u) - exp(x1*t) + y)*((u*u)*exp(x2*u)));
 end
 F_prime_x = [ sum3 sum4; ...
 sum4 sum5];
 yvect = xvect - inv(F_prime_x)*Fvect_x ;
 y1 = yvect(1);
 y2 = yvect(2);
 sum6 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);
 u = Amat(j,2);
 sum6 = sum6 + ((-exp(y1*t) - exp(y2*u) + y)*(-exp(y1*t)*t));
 end
 sum7 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);
 u = Amat(j,2);
 sum7 = sum7 + ((-exp(y1*t) - exp(y2*u) + y)*(-exp(y2*u)*u));
 end
 Fvect_y = [sum6; sum7];
 zvect = xvect - inv(F_prime_x)*( Fvect_x + Fvect_y );
 xvect = zvect;
 it = num2str(iter);
 message = ['at iteration ', it, ', we have [x, y] ='];
 disp(message)
 disp(round(xvect,10))
end
