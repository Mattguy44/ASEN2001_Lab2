%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will run a 3D analysis on an inputed truss system. 
%
% Pierce Costello, Braden Barkemeyer - 10/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = truss_3D(inputfile,outputfile)

% read input file
[joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs] = readinput(inputfile);

% run variables through force analysis
[barforces,reacforces] = forceanalysis3D(joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs);

% write output
writeoutput3D(outputfile,inputfile,barforces,reacforces,joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs);

% define radius of bars
rad_bars = 3*[0.025,0.04,0.05];

% plot 3D truss
plottruss3D(joints,connectivity,barforces,reacjoints,rad_bars,[1 1 0 0]);

end


