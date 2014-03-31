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


% compile related DPM code
compile;

% compile policy code
mexcmd = 'mex -outdir bin';

mexcmd = [mexcmd ' -O'];
mexcmd = [mexcmd ' CXXOPTIMFLAGS="-O3 -DNDEBUG -fomit-frame-pointer"'];
mexcmd = [mexcmd ' LDOPTIMFLAGS="-O3"'];

mexcmd = [mexcmd ' CXXFLAGS="\$CXXFLAGS -Wall"'];
mexcmd = [mexcmd ' LDFLAGS="\$LDFLAGS -Wall"'];

mexcmd = [mexcmd ' mex/detect_policy_pca.cc mex/policy.cc mex/model_pca_policy.cc'];

eval(mexcmd);