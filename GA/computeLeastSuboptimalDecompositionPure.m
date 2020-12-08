clear;
close all;
clc;

%% 

addpath('utils');
load('data/CartPoleSystem.mat');

fun = @(x) computeLQRMeasurePure(sys, reshape(x(1:sys.U_DIMS^2), sys.U_DIMS, sys.U_DIMS), ...
                                  reshape(x((sys.U_DIMS^2 + 1):(sys.U_DIMS^2 + sys.U_DIMS*sys.X_DIMS)), sys.U_DIMS, sys.X_DIMS), ...
                                  x((sys.U_DIMS^2 + sys.U_DIMS*sys.X_DIMS)+1));
nvars = sys.U_DIMS^2 + sys.X_DIMS*sys.U_DIMS + 1;
lb = zeros(1, nvars);
ub = ones(1, nvars);
Aeq = [];
beq = [];
nonlcon = [];
IntCon = linspace(1,nvars,nvars);
options.Display = 'iter';
nconst = 4*(sys.U_DIMS) + sys.U_DIMS + 2*sys.X_DIMS;
A = zeros(nconst, nvars);
b = zeros(nconst, 1);

% Unique order constraint
for ii=1:1:sys.U_DIMS
    % One rank is assigned to only one input
    A(ii, ((ii-1)*sys.U_DIMS + 1):(ii*sys.U_DIMS)) = ones(1, sys.U_DIMS);
    b(ii) = 1;
    A(sys.U_DIMS + ii, ((ii-1)*sys.U_DIMS + 1):(ii*sys.U_DIMS)) = -ones(1, sys.U_DIMS);
    b(sys.U_DIMS + ii) = -1;
    
    % Each input is assigned a rank
    A(2*sys.U_DIMS + ii, ii + sys.U_DIMS*linspace(0,sys.U_DIMS-1,sys.U_DIMS)) = ones(1, sys.U_DIMS);
    b(2*sys.U_DIMS + ii) = 1;
    A(3*sys.U_DIMS + ii, ii + sys.U_DIMS*linspace(0,sys.U_DIMS-1,sys.U_DIMS)) = -ones(1, sys.U_DIMS);
    b(3*sys.U_DIMS + ii) = -1;
    
    % Each input must have at least one state variable
    A(4*sys.U_DIMS + ii, sys.U_DIMS^2 + ii + sys.U_DIMS*linspace(0,sys.X_DIMS-1,sys.X_DIMS)) = -ones(1, sys.X_DIMS);
    b(4*sys.U_DIMS + ii) = -1;
end

for jj=1:1:sys.X_DIMS
    % Each state variable must be assigned to ony one input
    A(5*sys.U_DIMS + jj, ...
      (sys.U_DIMS^2 + ((jj-1)*sys.U_DIMS + 1)):(sys.U_DIMS^2 + jj*sys.U_DIMS)) = ones(1, sys.U_DIMS);
    b(5*sys.U_DIMS + jj) = 1;
    A(5*sys.U_DIMS + sys.X_DIMS + jj, ...
      (sys.U_DIMS^2 + ((jj-1)*sys.U_DIMS + 1)):(sys.U_DIMS^2 + jj*sys.U_DIMS)) = -ones(1, sys.U_DIMS);
    b(5*sys.U_DIMS + sys.X_DIMS + jj) = -1;
end

x = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);

% Extract decomposition
r = reshape(x(1:sys.U_DIMS^2), sys.U_DIMS, sys.U_DIMS);
s = reshape(x((sys.U_DIMS^2 + 1):(sys.U_DIMS^2 + sys.U_DIMS*sys.X_DIMS)), sys.U_DIMS, sys.X_DIMS);
d = x(end);

err_lqr = computeLQRMeasure(sys,r,s,d);