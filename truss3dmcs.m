    function truss3dmcs(inputfile)
% function truss2d(inputfile,outputfile)
%
% Stochastic analysis of 2-D statically determinate truss by
% Monte Carlo Simulation. Only positions and strength of joints 
% treated as random variables
%
% Assumption: variation of joint strength and positions described 
%             via Gaussian distributions
% 
%             joint strength : mean = 15
%                              coefficient of varation = 0.1
%             joint position : 
%                              coefficient of varation = 0.01
%                              (defined wrt to maximum dimension of truss)
%
%             number of samples is set to 1e4
%
% Input:  inputfile  - name of input file
%
% Author: Kurt Maute for ASEN 2001, Oct 13 2012

% parameters
jstrmean   = 4.8;    % mean of joint strength
jstrcov    = 0.4;   % coefficient of variation of joint strength
jposcov    = 0.01;  % coefficient of variation of joint position
numsamples = 1e3;   % number of samples                                                  %changed from e5

% read input file
[joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs]=readinput(inputfile);

% determine extension of truss
ext_x=max(joints(:,1))-min(joints(:,1));   % extension in x-direction
ext_y=max(joints(:,2))-min(joints(:,2));   % extension in y-direction
ext_z=max(joints(:,3))-min(joints(:,3));   % extension in z-direction
ext  =max([ext_x,ext_y,ext_z]);

% loop overall samples
numjoints=size(joints,1);       % number of joints
maxforces=zeros(numsamples,1);  % maximum bar forces for all samples
maxreact=zeros(numsamples,1);   % maximum support reactions for all samples
failure=zeros(numsamples,1);    % failure of truss

for is=1:numsamples 
    
    % generate random joint strength limit
    varstrength = (jstrcov*jstrmean)*randn(1,1);
    
    jstrength = jstrmean + varstrength;
    
    % generate random samples
    varjoints = (jposcov*ext)*randn(numjoints,3);
    
    % perturb joint positions
    randjoints = joints + varjoints;
    
    % compute forces in bars and reactions
    [barforces,reacforces] = forceanalysis(randjoints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs);
    
    % determine maximum force magnitude in bars and supports
    maxforces(is) = max(abs(barforces));
    maxreact(is)  = max(abs(reacforces));
    
    % determine whether truss failed
    failure(is) = maxforces(is) > jstrength || maxreact(is) > jstrength;
end

figure(1);
subplot(1,2,1);
hist(maxforces,30);
title('Histogram of maximum bar forces');
xlabel('Magnitude of bar forces');
ylabel('Frequency');

subplot(1,2,2);
hist(maxreact,30);
title('Histogram of maximum support reactions');
xlabel('Magnitude of reaction forces');
ylabel('Frequency');

fprintf('\nFailure probability : %e \n\n',sum(failure)/numsamples);

end

function [joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs]=readinput(inputfile)
% function [joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs]=readinput(inputfile)
%
% read input file
%
% input: inputfile - name of input file
%
% output: joints       - coordinates of joints
%         connectivity - connectivity 
%         reacjoints   - joint id where reaction acts on
%         reacvecs     - unit vector associated with reaction force
%         loadjoints   - joint id where external load acts on
%         loadvecs     - load vector
%
% Author: Kurt Maute, Sept 21 2011

% open inputfile
fid = fopen(inputfile);

if fid<0
    error('inputfile does not exist')
end

% initialze counters and input block id
counter = 0;
input_blk = 1;

% read first line
line=fgetl(fid);

% read input file
while line > 0
    
    % check if comment
    if strcmp(line(1),'#')
        % read next line and continue
        line = fgetl(fid);
        continue;
    end
    
    switch input_blk
        
        case 1 % read number of joints, bars, reactions, and loads
            
            input_par = sscanf(line,'%d%d%d%d%d');
            
            numjoints = input_par(1);
            numbars   = input_par(2);
            numreact  = input_par(3);
            numloads  = input_par(4);
            
            % check for correct number of reaction forces
            if numreact~=6; error('incorrect number of reaction forces');end
            
            % initialize arrays
            joints       = zeros(numjoints,3);
            connectivity = zeros(numbars,2);
            reacjoints   = zeros(numreact,1);
            reacvecs     = zeros(numreact,3);
            loadjoints   = zeros(numloads,1);
            loadvecs     = zeros(numloads,3);
            
            % check whether system satisfies static determiancy condition
            if 3*numjoints ~= 2*numbars
                error('truss is not statically determinate');
            end

            % expect next input block to be joint coordinates
            input_blk = 2;
            
        case 2 % read coordinates of joints
            
            % increment joint id
            counter = counter + 1;
            
            % read joint id and coordinates;
            tmp = sscanf(line,'%d%e%e');
            
            % extract and check joint id
            jointid = tmp(1);
            if jointid > numjoints || jointid<1
                error('joint id number need to be smaller than number of joints and larger than 0');
            end
            
            % store coordinates of joints
            joints(jointid,:) = tmp(2:4);
            
            % expect next input block to be connectivity
            if counter == numjoints
                input_blk  = 3;
                counter = 0;
            end
            
        case 3 % read connectivity of bars
            
            % increment bar id
            counter = counter + 1;
            
            % read connectivity;
            tmp = sscanf(line,'%d%d%d');
            
            % extract bar id number and check
            barid = tmp(1);
            if barid > numbars || barid<0
                error('bar id number needs to be smaller than number of bars and larger than 0');
            end
            
            % check joint ids
            if max(tmp(2:3)) > numjoints || min(tmp(2:3)) < 1
                error('joint id numbers need to be smaller than number of joints and larger than 0');
            end
            
            % store connectivity
            connectivity(barid,:) = tmp(2:3);
            
            % expect next input block to be reaction forces
            if counter == numbars
                input_blk = 4;
                counter = 0;
            end
            
        case 4 % read reaction force information
            
            % increment reaction id
            counter = counter + 1;
            
            % read joint id and unit vector of reaction force;
            tmp = sscanf(line,'%d%e%e');
            
            % extract and check joint id
            jointid=tmp(1);
            if jointid>numjoints || jointid<1
                error('joint id number need to be smaller than number of joints and larger than 0');
            end
            
            % extract untit vector and check length
            uvec = tmp(2:4);
            uvec = uvec/norm(uvec);
            
            % store joint id and unit vector
            reacjoints(counter)  = jointid;
            reacvecs(counter,:)  = uvec;
            
            % expect next input block to be external loads
            if counter == numreact
                input_blk  = 5;
                counter = 0;
            end
            
        case 5 % read external load information
            
            % increment reaction id
            counter = counter + 1;
            
            % read joint id and unit vector of reaction force;
            tmp = sscanf(line,'%d%e%e');
            
            % extract and check joint id
            jointid = tmp(1);
            if jointid > numjoints || jointid < 1
                error('joint id number need to be smaller than number of joints and larger than 0');
            end
            
            % extract force vector
            frcvec = tmp(2:4);
            
            % store joint id and unit vector
            loadjoints(counter) = jointid;
            loadvecs(counter,:) = frcvec;
            
            % expect no additional input block
            if counter == numloads
                input_blk  = 99;
                counter = 0;
            end
            
        otherwise
            %fprintf('warning: unknown input: %s\n',line);
    end
    
    % read next line
    line=fgetl(fid);
end

% close input file
fclose(fid);

end

function [barforces,reacforces]=forceanalysis(joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs)
% function [barforces,reacforces]=forceanalysis(joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs)
%
% compute forces in bars and reaction forces
%
% input:  joints       - coordinates of joints
%         connectivity - connectivity 
%         reacjoints   - joint id where reaction acts on
%         reacvecs     - unit vector associated with reaction force
%         loadjoints   - joint id where external load acts on
%         loadvecs     - load vector
%
% output: barforces    - force magnitude in bars
%         reacforces   - reaction forces
%
% Author: Kurt Maute, Sept 21 2011
% Updated: Braden Barkemeyer, Oct. 9, 2017
%
readinput('test3d_1.inp')
% extract number of joints, bars, reactions, and loads
numjoints = size(joints,1);
numbars   = size(connectivity,1);
numreact  = size(reacjoints,1);
numloads  = size(loadjoints,1);

% number of equations
numeqns = 3 * numjoints;

% allocate arrays for linear system
Amat = zeros(numeqns);
bvec = zeros(numeqns,1);

% build Amat - loop over all joints

for i=1:numjoints
    
   % equation id numbers
   idx = 3*i-2;
   idy = 3*i-1;
   idz = 3*i;
   % get all bars connected to joint
   [ibar,ijt]=find(connectivity==i);
   
   % loop over all bars connected to joint
   for ib=1:length(ibar)
       
       % get bar id
       barid=ibar(ib);
       
       % get coordinates for joints "i" and "j" of bar "barid"
       joint_i = joints(i,:);
       if ijt(ib) == 1
           jid = connectivity(barid,2);
       else
           jid = connectivity(barid,1);
       end
       joint_j = joints(jid,:);
       
       % compute unit vector pointing away from joint i
       vec_ij = joint_j - joint_i;
       uvec   = vec_ij/norm(vec_ij);
       
       % add unit vector into Amat
       Amat([idx idy idz],barid)=uvec;
   end
end

% build contribution of support reactions 
for i=1:numreact
    
    % get joint id at which reaction force acts
    jid=reacjoints(i);

    % equation id numbers
    idx = 3*jid-2;
    idy = 3*jid-1;
    idz = 3*jid;

    % add unit vector into Amat
    Amat([idx idy idz],numbars+i)=reacvecs(i,:);
end

% build load vector
for i=1:numloads
    
    % get joint id at which external force acts
    jid=loadjoints(i);

    % equation id numbers
    idx = 3*jid-2;
    idy = 3*jid-1;
    idz = 3*jid;

    % add unit vector into bvec (sign change)
    bvec([idx idy idz])=-loadvecs(i,:);
end

% check for invertability of Amat
if rank(Amat) ~= numeqns
    error('Amat is rank defficient: %d < %d\n',rank(Amat),numeqns);
end

% solve system
xvec=Amat\bvec;

% extract forces in bars and reaction forces
barforces=xvec(1:numbars);
reacforces=xvec(numbars+1:end);

end

function writeoutput(outputfile,inputfile,barforces,reacforces,joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs)
% writeoutput(outputfile,inputfile,barforces,reacforces,joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs);
%
% output analysis results
%
% input:  outputfile   - name of output file
%         inputfile    - name of input file
%         barforces    - force magnitude in bars
%         reacforces   - reaction forces
%         joints       - coordinates of joints
%         connectivity - connectivity 
%         reacjoints   - joint id where reaction acts on
%         reacvecs     - unit vector associated with reaction force
%         loadjoints   - joint id where external load acts on
%         loadvecs     - load vector
%
%
% Author: Kurt Maute, Sept 21 2011
%Updated: Braden Barkemeyer, Oct. 9, 2017

% open output file
fid=fopen(outputfile,'w');

% write header
fprintf(fid,'3-D Truss analysis\n');
fprintf(fid,'------------------\n\n');
fprintf(fid,'Date: %s\n\n',datestr(now));

% write name of input file
fprintf(fid,'Input file: %s\n\n',inputfile);

% write coordinates of joints
fprintf(fid,'Joints:         Joint-id  x-coord.   y-coord.    z-coord.\n');
for i=1:size(joints,1)
    fprintf(fid,'%17d %12.2f %12.2f %12.2f\n',i,joints(i,1),joints(i,2),joints(i,3));
end
fprintf(fid,'\n\n');

% write external loads
fprintf(fid,'External loads: Joint-id  Force-x      Force-y     Force-z\n');
for i=1:size(loadjoints,1)
    fprintf(fid,'%17d %12.2f %12.2f %12.2f\n',loadjoints(i),loadvecs(i,1),loadvecs(i,2), loadvecs(i,3));
end
fprintf(fid,'\n');
    
% write connectivity and forces
fprintf(fid,'Bars:      Bar-id    Joint-i       Joint-j         Force   (T,C)\n');
for i=1:size(connectivity,1)
    if barforces(i)>0;tc='T';else tc='C';end
    fprintf(fid,'%17d   %7d %12d    %12.3f     (%s)\n',...
        i,connectivity(i,1),connectivity(i,2),abs(barforces(i)),tc);
end
fprintf(fid,'\n');

% write connectivity and forces
fprintf(fid,'Reactions:    Joint-id  Uvec-x    Uvec-y     Uvec-z    Force\n');
for i=1:size(reacjoints,1)
    fprintf(fid,'%17d %10.2f %10.2f %10.2f %10.3f\n',reacjoints(i),reacvecs(i,1),reacvecs(i,2), reacvecs(i,3),reacforces(i));
end

% close output file
fclose(fid);

end