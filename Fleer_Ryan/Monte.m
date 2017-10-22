function Monte(inputfile)
% function truss3d_mass(inputfile,outputfile)
%
% Stochastic analysis of 2-D statically determinate truss by
% Monte Carlo Simulation. Only positions and strength of joints 
% treated as random variables
%
% Assumption: variation of joint strength and positions described 
%             via Gaussian distributions
% 
%             joint strength : mean = 4.8 N
%                              coefficient of varation = 0.4/4.8 = 0.08
%             joint position : 
%                              coefficient of varation = 0.01
%                              (defined wrt to maximum dimension of truss)
%
%             number of samples is set to 1e5
%
% Input:  inputfile  - name of input file
%
% Author: Kurt Maute for ASEN 2001, Oct 13 2012
% Modified: Ethan Fleer Oct. 21, 2017
% Modified: Matthew Ryan Oct. 22, 2017

% parameters
Joint_stren_mean   = 4.8;    % mean of joint strength
Joint_stren_cov    = 0.08;   % coefficient of variation of joint strength
Joint_posis_cov    = 0.01;  % coefficient of variation of joint position
Number_samples = 1e5;   % number of samples

% read input file
[Joints,connectivity,reacjoints,reacvecs,rhobar,jointweight,maxload]=readinput_mass(inputfile);
% determine extension of truss
ext_x=max(Joints(:,1))-min(Joints(:,1));   % extension in x-direction
ext_y=max(Joints(:,2))-min(Joints(:,2));   % extension in y-direction
ext_z=max(Joints(:,3))-min(Joints(:,3));   % extension in z-direction
ext  =max([ext_x,ext_y,ext_z]);

% loop overall samples
Number_Joints=size(Joints,1);       % number of joints
Max_Forces=zeros(Number_samples,1);  % maximum bar forces for all samples
Max_React=zeros(Number_samples,1);   % maximum support reactions for all samples
Failure=zeros(Number_samples,1);    % failure of truss

for is=1:Number_samples
    
    % This creates a random joint strength limit
    Rand_strength = (Joint_stren_cov*Joint_stren_mean)*randn(1,1);
    J_strength = Joint_stren_mean + Rand_strength;
    
    % This creates random samples
    Rand_joints = (Joint_posis_cov*ext)*randn(Number_Joints,3);    
    
    % This randomizes the joint positions
    randjoints = Joints + Rand_joints;    
    
    % This finds the forces in bars and reactions
    [barforces,reacforces,~,~] = forceanalysis3D_mass(randjoints,connectivity,reacjoints,reacvecs,rhobar,jointweight,maxload);    
    
    % This finds the maximum force in the bars and supports
    Max_Forces(is) = max(abs(barforces));
    Max_React(is)  = max(abs(reacforces));    
    
    % This determines whether or not the truss fails
    Failure(is) = Max_Forces(is) > J_strength || Max_React(is) > J_strength;
end

%These plot the data
figure(1);
subplot(1,2,1);
hist(Max_Forces,30);
title('Maximum Bar Forces');
xlabel('Magnitude of Bar Forces (N)');
ylabel('Frequency');

subplot(1,2,2);
hist(Max_React,30);
title('Histogram of Maximum Support Reactions');
xlabel('Magnitude of Reaction Forces (N)');
ylabel('Frequency');

fprintf('\nFailure probability : %e \n\n',sum(Failure)/Number_samples);

end