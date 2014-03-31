function demo_pca_policy()
% Run policy demo.

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

clear all
close all
clc

startup;
fprintf(['DEMO for Active Deformable Part Models, \n',...
    'The demo include detections on two classes: person and car,\n',...
    'only the top detection will be displayed\n\n']);

fprintf('Compiling mex files...\n');
pca_policy_compile;

classes = {'person','car'};

for i = 1:2
    cls = classes{i};
    fprintf('\nDetecting %s...\n',cls);

    images = dir(sprintf('images/%s/',cls));
    for j = 1:numel(images)
        if length(images(j).name)>3 && strcmp(images(j).name(end-3:end),'.jpg')
            detect_object([sprintf('images/%s/',cls) images(j).name],cls);
            fprintf('----------------------------------------\n\n');
            pause(1);
        end
    end
end

