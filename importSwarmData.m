% ************************************************************************
% ********************       importSwarmData.m        ********************
% ************************************************************************
% Created by Kevin Rio, Stephane Bonneaud & Zach Page
% Last modified: 24 July 2014 by Kevin Rio
%   - Fixed a bug caused by inconsistent filenames output from QTM; which
%     sometimes have a 'qtm_' prefix and other times do not.
%   - Fixed the multi-file capability, so it is now possible to run the
%     script on a directory and batch process multiple files at once.
%     Output is separate files (F*.mat, where * is the original filename) 
%     for each trial. 
%   - Fixed bug which computed speed and heading incorrectly, which
%     subsequently caused problems in quant measures like polarization. 


%%
% ========================================================================
%  --- IMPORT RAW POSITION --
% Takes Sayles Swarm data from Qualisys Track Manager (QTM) and returns x,y
% position. These data are stored in a structure that is organized by trial
% and helmet number, e.g. Traj(2,16) is the trajectory data from trial 2
% for helmet 16.
% ========================================================================
clear all; close all;

% Initialize variables. Note that line 20 assumes filenames end in 
% *(3D).mat; this can be changed if desired. Also, filenames cannot contain
% dashes ("-"); this is required by Matlab. 
dataList = dir('*.mat');                                                    % list of Qualisys files in current directory
nTrial = length(dataList);                                                  % number of trials to import 
Traj = struct;                                                              % structure of x,y,t 
fileList = cell(nTrial,1);                                                  % list of filenames for each trial

% Cycle through all trajectories, identify helmet number / tracker letter,
% and import the data into the "Traj" structure. 
for i = 1:nTrial
    
    % Load a trial from the list of trials imported from the directory.
    fileName = dataList(i).name(1:end-4);                                   % remove '.mat' suffix
    dataName = regexprep(fileName,'-','_');                                 % replace hyphens ("-") with underscores ("_")
    load(fileName)
    if exist(dataName) == 0                                                 % if dataName does not exist, add "qtm_" prefix
        dataName = strcat('qtm_',dataName);
    end

    % Import data into "Trajectories" structure. 
    estimates = cell(1,20);
    num_trajectories = eval(strcat(dataName,'.Trajectories.Labeled.Count')); % # trajectories (# helmets x 5)
    for j = 1:num_trajectories
        
        % Determine the helmet number and tracker letter of the trajectory
        % that the loop is currently on. 
        suffix = ['.Trajectories.Labeled.Labels{',num2str(j),'}(2:end-1)'];  % end of path to helmet
        helmetNb = eval(strcat(dataName,suffix));                            % helmet number
        suffix = ['.Trajectories.Labeled.Labels{',num2str(j),'}(end)'];      % end of path to tracker
        trackLetter = eval(strcat(dataName,suffix));                         % tracker letter
        
        % Import each tracker (a-e) for the current helmet. 
        if trackLetter == 'a'                                                
            suffix = ['.Trajectories.Labeled.Data(',num2str(j),',1:3,:)'];
            trackers.a = squeeze(eval(strcat(dataName,suffix)));
        elseif trackLetter == 'b'
            suffix = ['.Trajectories.Labeled.Data(',num2str(j),',1:3,:)'];
            trackers.b = squeeze(eval(strcat(dataName,suffix)));
        elseif trackLetter == 'c'
            suffix = ['.Trajectories.Labeled.Data(',num2str(j),',1:3,:)'];
            trackers.c = squeeze(eval(strcat(dataName,suffix)));
        elseif trackLetter == 'd'
            suffix = ['.Trajectories.Labeled.Data(',num2str(j),',1:3,:)'];
            trackers.d = squeeze(eval(strcat(dataName,suffix)));
        elseif trackLetter == 'e'
            suffix = ['.Trajectories.Labeled.Data(',num2str(j),',1:3,:)'];
            trackers.e = squeeze(eval(strcat(dataName,suffix)));
        end
        
        % If you've reached the 5th, 10th, etc. tracker, you've reached the
        % end of that helmet, so compute head position. 
        if rem(j,5) == 0
            estimates{j/5} = extractHeadTrajectory(helmetNb,trackers);
            H = estimateHead(estimates{j/5});
            Traj(i,j/5).x = H(:,1)/1000;
            Traj(i,j/5).y = H(:,2)/1000;
            Traj(i,j/5).t = (1/60:1/60:length(H)/60)';
        end
    end
    
    % Add NaN's if number of helmets < 20.
    nHelm = (num_trajectories/5);
    nFrame = length(Traj(i,1).x);
    for helm = floor((nHelm+1)):20
        Traj(i,helm).x = NaN(nFrame,1);
        Traj(i,helm).y = NaN(nFrame,1);
        Traj(i,helm).t = (1/60:1/60:nFrame/60)';
    end
    
    % Store filename for each trial.
    fileList{i} = fileName;
    
    % Clear temporary variables from the workspace. 
    clearvars -except Traj fileList dataList i nTrial;
    
end

% Clear unneeded variables and save data. 
clearvars -except Traj fileList;
'Finished 1 of 3'

%%
% ========================================================================
%  --- IMPORT FILTERED POSITION --
% Takes raw Sayles Swarm data (generated from importRawPos) and applies a
% Butterworth low-pass filter. Also computes speed and heading time series.
% ========================================================================

% Apply Butterworth low-pass filter.
nTrial = size(Traj,1);
nHelm = size(Traj,2);
for iTrial = 1:nTrial
    
    nFrame = length(Traj(iTrial).x);
    
    % Divide data into usable 'segments' separated by NaN values.
    clear fragment
    fragment = cell(nHelm,1);
    for jHelm = 1:nHelm
        count = 1;
        if isempty(Traj(iTrial,jHelm).x(1)) == 0 && isnan(Traj(iTrial,jHelm).x(1)) == 0 % first segment begins at 1
            fragment{jHelm}(count,1) = 1;                                   % *unless value is NaN
        end
        for frame = 2:nFrame
            if isnan(Traj(iTrial,jHelm).x(frame)) == 1 && ...               % end of segment
                    isnan(Traj(iTrial,jHelm).x(frame-1)) == 0
                fragment{jHelm}(count,2) = frame-1;
                count = count+1;
            elseif isnan(Traj(iTrial,jHelm).x(frame)) == 0 && ...           % beginning of next segment
                    isnan(Traj(iTrial,jHelm).x(frame-1)) == 1
                fragment{jHelm}(count,1) = frame;
            end
        end
    end
    
    % Remove fragments that have only one data point (i.e. fragments that
    % aren't fragments!).
    for jHelm = 1:nHelm
        nFragment = size(fragment{jHelm},1);
        for iFragment = 1:nFragment
            if length(fragment{jHelm}) == 1
                fragment{jHelm} = [];
            end
            if isempty(fragment{jHelm}) == 0 && fragment{jHelm}(iFragment,2) == 0
                fragment{jHelm}(iFragment,:) = [];
            end
        end
    end
    
    % Filter each segment independently.
    for jHelm = 1:nHelm
        nFragment = size(fragment{jHelm},1);
        Traj(iTrial,jHelm).xF = NaN(length(Traj(iTrial,jHelm).x),1);
        Traj(iTrial,jHelm).yF = NaN(length(Traj(iTrial,jHelm).y),1);
        
        % Use linear extrapolation to 'pad' the data to avoid endpoint
        % error when filtering.
        for iFragment = 1:nFragment
            beg = fragment{jHelm}(iFragment,1);
            fin = fragment{jHelm}(iFragment,2);
            if fin-beg > 12 % 12 = 3 x 4(th order)
                Traj(iTrial,jHelm).xF(beg:fin) = bFiltExtrap(60,Traj(iTrial,jHelm).x(beg:fin));
                Traj(iTrial,jHelm).yF(beg:fin) = bFiltExtrap(60,Traj(iTrial,jHelm).y(beg:fin));
            end
        end     
    end
end
    
% Compute speed and heading, using filtered position.
for iTrial = 1:nTrial
    for jHelm = 1:nHelm
        [s,h] = xyToSpeedHeading(Traj(iTrial,jHelm).xF,Traj(iTrial,jHelm).yF);
        Traj(iTrial,jHelm).spd = s;
        Traj(iTrial,jHelm).hdn = h;
    end
end

% Clear unneeded variables and save data. 
clearvars -except Traj nTrial nHelm fileList;
'Finished 2 of 3'

%%
% ========================================================================
%  --- GENERATE QUANTITATIVE MEASURES --
% Takes filtered Sayles Swarm data (generated from importFiltPos) and
% produces a number of quantitative measures that are useful in subsequent
% analyses. These include a time series of the swarm's density;
% polarization and angular momentum measures of cohesion; and a recording
% of which helmet is nearest to the center of the swarm at each time step. 
% ========================================================================

% Initialize variables.
Quant = struct;

% Compute density by taking the area of the convex hull that encloses all
% (tracked) pedestrians in the swarm.
x = NaN(nHelm,1);
y = NaN(nHelm,1);
for iTrial = 1:nTrial
    for st = 1:length(Traj(iTrial,1).x)
        for jHelm  = 1:nHelm
            x(jHelm,1) = Traj(iTrial,jHelm).x(st);
            y(jHelm,1) = Traj(iTrial,jHelm).y(st);
        end
        x(isnan(x)==1) = []; y(isnan(y)==1) = [];                           % remove NaN values
        Quant(iTrial,1).Nped(st,1) = length(x);                             % number of pedestrians
        if length(x) > 3                                                    % must have 4+ points
            vi = convhull(x,y);                                             % compute convex hull
            area = polyarea(x(vi),y(vi));                                   % compute area
            Quant(iTrial,1).dens(st,1) = length(x)/area;                    % compute density
        else
            Quant(iTrial,1).dens(st,1) = NaN;
        end
    end
end

% Determine which helmet number is closest to the center of the swarm.
for iTrial = 1:nTrial
    
    % find swarm center of mass (COM)
    nFrame = length(Traj(iTrial,1).x);
    x = NaN(nFrame,nHelm);
    y = NaN(nFrame,nHelm);
    for jHelm = 1:nHelm
        x(:,jHelm) = Traj(iTrial,jHelm).x;
        y(:,jHelm) = Traj(iTrial,jHelm).y;
    end
    Mx = NaN(nFrame,1); My = NaN(nFrame,1);
    filtMx = NaN(nFrame,1); filtMy = NaN(nFrame,1);
    for frame = 1:nFrame
        Mx(frame,1) = nanmean(x(frame,:));
        My(frame,1) = nanmean(y(frame,:));
    end
    filtMx = bFiltExtrap(60,Mx(1:end-1));    % filtered COM location
    filtMy = bFiltExtrap(60,My(1:end-1));
    
    % find center helmet
    center_helm = NaN(nFrame,1);
     for frame = 1:nFrame-1
        min_separation = 100;
         for jHelm = 1:20
             separation = pdist([Traj(iTrial,jHelm).x(frame) Traj(iTrial,jHelm).y(frame); ...
                 filtMx(frame,1) filtMy(frame,1)]);
             if separation < min_separation
                 min_separation = separation;
                 center_helm(frame,1) = jHelm;
             end
         end
     end
     center_helm(nFrame,1) = center_helm(nFrame-1,1);
     Quant(iTrial,1).ctrHelm = center_helm;

end

% Compute polarization, a measure of linear cohesion.
for iTrial = 1:nTrial
    nHelm = max(Quant(iTrial).Nped);
    nFrame = length(Traj(iTrial,1).x);
    
    % generate velocity vectors (vx,vy)
    vx = NaN(nFrame,nHelm);
    vy = NaN(nFrame,nHelm);
    v = cell(nHelm,1);
    nv = cell(nHelm,1);
    nvx = NaN(nFrame,nHelm);
    nvy = NaN(nFrame,nHelm);
    for jHelm = 1:nHelm
         vx(:,jHelm) = Traj(iTrial,jHelm).spd.*cos(Traj(iTrial,jHelm).hdn);
         vy(:,jHelm) = Traj(iTrial,jHelm).spd.*sin(Traj(iTrial,jHelm).hdn);
         v{jHelm} = horzcat(vx(:,jHelm),vy(:,jHelm));
         for frame = 1:nFrame
             nv{jHelm}(frame,:) = v{jHelm}(frame,:)/norm(v{jHelm}(frame,:));
         end
         nvx(:,jHelm) = nv{jHelm}(:,1);
         nvy(:,jHelm) = nv{jHelm}(:,2);
    end
     
     % compute polarization
     hcount = NaN(nFrame,1);
     pv = NaN(nFrame,2);
     polarization = NaN(nFrame,1);
     for frame = 1:nFrame
         hcount(frame,1) = nHelm-sum(isnan(nvx(frame,:)));               % number of helmets present on each time step
         if hcount(frame,1) > 0
             pv(frame,1) = nansum(nvx(frame,:));                            % vector sum of nvv's
             pv(frame,2) = nansum(nvy(frame,:));
             polarization(frame,1) = norm(pv(frame,:))/hcount(frame,1);     % polarization
         end
     end
     
     % write to file
     Quant(iTrial).polarization = polarization;
    
end

% Compute angular momentum, a measure of circular cohesion. 
for iTrial = 1:nTrial
    nHelm = max(Quant(iTrial).Nped);
    nFrame = length(Traj(iTrial,1).x);
    
    % find center of mass (COM) of swarm
    fpx = NaN(nFrame,nHelm);                                          % filtered position (x)
    fpy = NaN(nFrame,nHelm);                                          % filtered position (y)
    vx = NaN(nFrame,nHelm);
    vy = NaN(nFrame,nHelm);
    for jHelm = 1:nHelm
        vx(:,jHelm) = Traj(iTrial,jHelm).spd.*cos(Traj(iTrial,jHelm).hdn);
        vy(:,jHelm) = Traj(iTrial,jHelm).spd.*sin(Traj(iTrial,jHelm).hdn);
        fpx(:,jHelm) = Traj(iTrial,jHelm).x;
        fpy(:,jHelm) = Traj(iTrial,jHelm).y;
    end
    ctrUnfilt = NaN(nFrame,2);                                           % center vector
    hcount = NaN(nFrame,1);
    for frame = 1:nFrame
        hcount(frame,1) = nHelm-sum(isnan(fpx(frame,:)));
        if hcount(frame,1) > 0
            ctrUnfilt(frame,1) = nansum(fpx(frame,:))/hcount(frame,1);      % COM (x)
            ctrUnfilt(frame,2) = nansum(fpy(frame,:))/hcount(frame,1);      % COM (y)
        end
    end
    ctr = NaN(nFrame,2);
    ctr(4:nFrame-5,1) = bFiltExtrap(60,ctrUnfilt(4:nFrame-5,1));      % filtered COM (x)
    ctr(4:nFrame-5,2) = bFiltExtrap(60,ctrUnfilt(4:nFrame-5,2));      % filtered COM (y)
    
    % compute vector from COM to each helmet
    rcx = NaN(nFrame,nHelm);                                          % vector from center to each helmet (x)
    rcy = NaN(nFrame,nHelm);                                          % vector from center to each helmet (y)
    rnorm = NaN(nFrame,nHelm);
    for jHelm = 1:nHelm
        rcx(:,jHelm) = Traj(iTrial,jHelm).x-ctr(:,1);
        rcy(:,jHelm) = Traj(iTrial,jHelm).y-ctr(:,2);
    end
    rnorm = sqrt(rcx(:,:).^2 + rcy(:,:).^2);
    rcx(:,:) = rcx(:,:)./rnorm(:,:);
    rcy(:,:) = rcy(:,:)./rnorm(:,:);
    
    % compute angular momentum measure
    rcz = zeros(nFrame,1);
    vz = zeros(nFrame,nHelm);
    cpd = NaN(nFrame,nHelm);
    for frame = 1:nFrame
        if hcount(frame,1) > 0
            for jHelm = 1:nHelm
                v = [vx(frame,jHelm) vy(frame,1) vz(frame,jHelm)];
                r = [rcx(frame,jHelm) rcy(frame,1) rcz(frame,1)];
                tmp = cross(v,r);
                cpd(frame,jHelm) = tmp(3);
            end
        end
    end
    momentum = NaN(nFrame,1);
    for frame = 1:nFrame
        momentum(frame,1) = norm(nansum(cpd(frame,:)))/hcount(frame,1);
    end
    
    % write to file
    Quant(iTrial).momentum = momentum;
    
end

% Clear unneeded variables and save data. 
T = Traj;
Q = Quant;
for i = 1:nTrial
    clear Traj;
    clear Quant;
    Traj = T(i,:);
    Quant = Q(i);
    save( strcat('F',fileList{i}) , 'Traj','Quant' )
end
'Finished 3 of 3!'

