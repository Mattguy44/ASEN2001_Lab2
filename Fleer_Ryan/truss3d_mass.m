function truss3d_mass(inputfile,outputfile)
% function truss2d(inputfile,outputfile)
%
% Analysis of 2-D statically determinate truss
%
% Input:  inputfile  - name of input file
%         outputfile - name of output file
%
% Author: Kurt Maute for ASEN 2001, Sept 21 2011
% Modified: Matthew Ryan, Oct. 16, 2017

% read input file
[joints,connectivity,reacjoints,reacvecs,rhobar,jointweight,maxload]=readinput_mass(inputfile);

% compute forces in bars and reactions
[barforces,reacforces,jointloads,jointstatus]=forceanalysis3D_mass(joints,connectivity,reacjoints,reacvecs,rhobar,jointweight,maxload);

% write outputfile
writeoutput3D_mass(outputfile,inputfile,barforces,reacforces,joints,connectivity,reacjoints,reacvecs,rhobar,jointweight,jointloads,jointstatus);

% plot truss
plottruss(joints,connectivity,barforces,reacjoints,3*[0.025,0.04,0.05],[1 1 0 0]);

end
