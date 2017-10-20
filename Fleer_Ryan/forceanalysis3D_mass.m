function [barforces,reacforces,jointloads,jointfailure]=forceanalysis3D_mass(joints,connectivity,reacjoints,reacvecs,rhobar,jointweight,maxload)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% compute forces in bars and reaction forces
%
% input:  joints       - coordinates of joints (in)
%         connectivity - connectivity 
%         reacjoints   - joint id where reaction acts on
%         reacvecs     - unit vector associated with reaction force (in)
%         rhobar       - linear density of members (lb/in)
%         jointweight  - average weight of a joint (lb)
%
% output: barforces    - force magnitude in bars
%         reacforces   - reaction forces
%
% Author: Kurt Maute, Sept 21 2011
% Modified: Matthew Ryan, Oct. 16, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract number of joints, bars, and reactions
numjoints = size(joints,1);
numbars   = size(connectivity,1);
numreact  = size(reacjoints,1);

% allocate arrays
weightvec = zeros(3*numjoints,1);
%jointfailure = repmat(' ',[numjoints 1]);

% define max design load
k = 1.7;
maxdesign = maxload/k;

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
       vec_ij_mag = norm(vec_ij);
       uvec   = vec_ij/vec_ij_mag;
       
       % add half weight of bar to joint
       weightvec(idz) = weightvec(idz) - 0.5*vec_ij_mag*rhobar;
       
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
for i=1:numjoints

    % equation id numbers
    idz = 3*i;

    % add weight vector into bvec (sign change)
    bvec(idz) = -weightvec(idz) + jointweight;
    
    % check if load on joint is too great
    b_mag = norm([bvec(idx) bvec(idy) bvec(idz)]);
    if b_mag > maxload
        jointfailure(i) = string('Exceeds Maximum Load');
    else
        if b_mag > maxdesign
            jointfailure(i) = string('Exceeds Allowable Load');
        end
    end
    
%     % extract loads on joints
%     jointloads(i) = [bvec(idx) bvec(idy) bvec(idz) jointfailure(i)];
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
jointloads = -bvec;
end