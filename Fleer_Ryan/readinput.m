%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function serves as the primary read input function that will be used
% in the 3-D truss analysis.
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
% Modified By: Ethan Fleer, Oct. 9, 2017
% Modified By: Matthew Ryan,  Oct. 16, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs] = readinput(inputfile)

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
            if 3*numjoints -6 ~= numbars
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

