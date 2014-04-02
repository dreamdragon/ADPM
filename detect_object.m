function [ ds ] = detect_object( impath, cls, DISPLAY )
%DETECT_OBJECT Detecting Object using Active Deformable Part Models
%   impath:     input image file
%   cls:        object class of interest
%   DISPLAY:    boolean flag for displaying top detection
%   ds:         detection bounding boxes

% AUTORIGHTS
% -------------------------------------------------------
% Copyright (C) 2013-2014 Menglong Zhu, Nikolay Atanasov
%                         Samarth Brahmbhatt
% 
% This file is part of the Active Deformable Part Models 
% code (http://cis.upenn.edu/~menglong/adpm.html)
% and is available under the terms of an MIT-like license
% provided in COPYING. Please retain this notice and
% COPYING if you use this file (or a portion of it) in
% your project.
% -------------------------------------------------------

    if nargin < 1
        disp('Please specify input image file');
        ds = [];
        return;
    end
    
    if nargin < 2
        disp('Please specify object class');
        ds = [];
        return;
    end
    
    if nargin < 3
        DISPLAY = 1;
    end
    
    % load original DPM model
    load(sprintf('models/VOC2010/%s_final.mat', cls));

    % load pca policy
    pca_pol = load(sprintf('models/%s_2010_pca_policy_%d_%d',cls, 20, 5));
    pca_policy = pca_pol.policy;

    %load full policy
    pol = load(sprintf('models/%s_2010_policy_%d_%d',cls, 50, 5));
    policy  = pol.policy;
    

    im = imread(impath);

    % Display image 
    if DISPLAY
        clf;
        imagesc(im);
        axis image;
        axis off;
    end
    
    % Perform object detection
    ds = detect(im, model, policy, pca_policy);

    % Display detection
    if DISPLAY
        showboxes(im, ds(1,1:4));
        title('ADPM top detections');
    end
end

% main detection function
function [ds] = detect(im, model, policy, pca_policy)

pca = 5;
load('external/pca.mat');
policy_model = grammar2simple(project_model(model, coeff, pca));

fprintf('Building the feature pyramid...');

th = tic();
pyra = featpyramid(double(im), model);
tF = toc(th);

fprintf('done\n');
fprintf('  --> Feature pyramid generation took %f seconds\n', tF);

fprintf('Computing detections with ADPM...');
[dP, bP, tP, nP] = pca_policy_detect(pyra, policy_model, policy, pca_policy);
fprintf('done\n');
fprintf('  --> ADPM detection took %f seconds\n', tP);

ds = getboxes(im, dP, bP);


end

function b = getboxes(image, det, all)
b = [];
if ~isempty(det)
  [det, all] = clipboxes(image, det, all);
  I = nms(det, 0.5);
  det = det(I,:);
  all = all(I,:);
  b = [det(:,1:4) all];
end

end
