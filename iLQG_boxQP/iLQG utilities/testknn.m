clear;
clc;
close all;

%% 
load('../iLQGDecomposed_test.mat', 'XCartFF', 'UCartFF', 'KCartFF');
Mdl = ExhaustiveSearcher(XCartFF');
% Mdl = KDTreeSearcher(XCartFF');
numPoints = 5000;
Xquery = [-1.5 + 3*rand(numPoints, 1), ...
          -3 + 6*rand(numPoints, 1), ...
          0 + 2*pi*rand(numPoints, 1), ...
          -3 + 6*rand(numPoints, 1)];
% Try knn
tic;
knnclosest = knnsearch(Mdl, Xquery);
time_knn = toc;

% Try direct look-up
dirclosest = zeros(numPoints, 1);
tic;
for ii=1:1:numPoints
    [~, dirclosest(ii)] = min(vecnorm(XCartFF - Xquery(ii,:)'), [], 2);
end
time_dir = toc;

% Try pdist2 look-up
XCartFFT = XCartFF';
tic;
[~, pdistclosestT] = min(pdist2(XCartFFT, Xquery));
time_pdistT = toc;

tic;
[~, pdistclosest] = min(pdist2(XCartFF', Xquery));
time_pdist = toc;