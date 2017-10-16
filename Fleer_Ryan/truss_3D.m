%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This function will run a 3D analysis on a truss system inputed via an 
% input data file 
%
% Pierce Costello, Braden Barkemeyer - 10/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = truss_3D(inputfile,outputfile)

% read input file
[joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs] = readinput('test3d_1.inp');

% run variables through force analysis
[barforces,reacforces] = forceanalysis3D(joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs);

% write output
writeoutput3D(outputfile,inputfile,barforces,reacforces,joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs);

% plot 3D truss
plottruss3D(joints,connectivity,barforces,reacjoints,3*[0.025,0.04,0.05],[1 1 0 0]);

end


