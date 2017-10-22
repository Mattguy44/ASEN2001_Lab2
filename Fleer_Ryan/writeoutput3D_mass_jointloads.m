%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Modified: Matthew Ryan, Oct. 20, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function writeoutput3D_mass(outputfile,inputfile,barforces,reacforces,joints,connectivity,reacjoints,reacvecs,rhobar,jointweight,jointloads,jointstatus)
%
% open output file
fid=fopen(outputfile,'w');

% write header
fprintf(fid,'3-D Truss analysis with weighted members\n');
fprintf(fid,'----------------------------------------\n\n');
fprintf(fid,'Date: %s\n\n',datestr(now));

% write name of input file
fprintf(fid,'Input file: %s\n\n',inputfile);

% write properties of materials
fprintf(fid, 'Linear weight density of members: %3.5f\n', rhobar);
fprintf(fid, 'Average joint weight: %3.5f\n\n', jointweight);

% write coordinates of joints
fprintf(fid,'Joints:         Joint-id  x-coord.     y-coord.     z-coord.\n');
for i=1:size(joints,1)
    fprintf(fid,'%17d %12.2f %12.2f %12.2f\n',i,joints(i,1),joints(i,2),joints(i,3));
end
fprintf(fid,'\n\n');

% write loads on joints
fprintf(fid, 'Joint loads:    Joint-id  Force-x      Force-y      Force-z\n');
for i = 1:size(joints,1)
    fprintf(fid, '%17d  %12.3f %12.3f %12.3f    %s\n', i, jointloads(i*3-2), jointloads(i*3-1), jointloads(i*3), jointstatus(i));
end
fprintf(fid,'\n');
    
% write connectivity and forces
fprintf(fid,'Bars:           Bar-id    Joint-i      Joint-j      Force    (T,C)\n');
for i=1:size(connectivity,1)
    if barforces(i)>0;tc='T';else tc='C';end
    fprintf(fid,'%17d   %7d %12d     %12.3f     (%s)\n',...
        i,connectivity(i,1),connectivity(i,2),abs(barforces(i)),tc);
end
fprintf(fid,'\n');

% write connectivity and forces
fprintf(fid,'Reactions:      Joint-id  Uvec-x       Uvec-y       Uvec-z      Force\n');
for i=1:size(reacjoints,1)
    fprintf(fid,'%17d %12.2f %12.2f %12.2f %12.3f\n',reacjoints(i),reacvecs(i,1),reacvecs(i,2), reacvecs(i,3),reacforces(i));
end

% close output file
fclose(fid);

end

