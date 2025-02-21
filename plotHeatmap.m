clear all;
close all;
%load F2_CS_1N_13;
%load F2_CS_1F_7;
load F2_CS_1F_20;
%load F2_CS_1N_8;
%load F2_CS_1S_18;
%load F2_CS_1S_23;
%load F2_CS_2F_12;
%load F2_CS_2F_19;
%load F6_SW_2_2;
%load F6_SW_2_4;
%load F2_SW_1_1;
%load F2_SW_2_2;
%load F4_SW_1_4;
%load F4_SW_2_3;
%load F6_SW_1_3;
%load F2_CFO_90-2_10;
fn = 'Helvetica';

% ========================================================================
% HEAT MAP #1: MEAN ABSOLUTE HEADING DIFFERENCE
% ========================================================================

count = 1;
nTrial = 1;
for iTrial = 1:nTrial
   nHelm = max(Quant(iTrial).Nped);
   nFrame = length(Traj(iTrial,1).x);
    
    % for all non-center helmets, compute position (x,y) and difference in heading
    for jFrame = 1:nFrame-1
        for kHelm = 1:nHelm
            cH = Quant(iTrial).ctrHelm(jFrame);
            if kHelm ~= cH;
                if isnan(Traj(iTrial,kHelm).hdn(jFrame,1)) == 0
                    
                    % difference in heading
                    angularDeviation(count,1) = (angleBn(Traj(iTrial,kHelm).hdn(jFrame,1),...
                        Traj(iTrial,cH).hdn(jFrame,1)));
                                        
                    % position (x,y) in cH-centered coordinates
                    x = Traj(iTrial,kHelm).x(jFrame)-Traj(iTrial,cH).x(jFrame);
                    y = Traj(iTrial,kHelm).y(jFrame)-Traj(iTrial,cH).y(jFrame);
                    th = Traj(iTrial,cH).hdn(jFrame);
                    xR(count,1) = x*cos(th)-y*sin(th);
                    yR(count,1) = x*sin(th)+y*cos(th);
                    dR(count,1) = sqrt(xR(count,1)^2+yR(count,1)^2);
                    
                    % take absolute value of deviation
                    angularDeviation(count,1) = abs(angularDeviation(count,1))*(180/pi);
                    
                    count = count+1;
                    
                end
            end
        end
    end
end

% divide data into grids
grid = 0.5;   % size of grid squares (in m)
gmin = floor(min(xR));  % smallest value
gmax = ceil(max(xR));   % largest value
angularDevList = cell(2*gmax);
count = 1;
gXcount = 1; 
for gX = gmin:grid:gmax
     gYcount = 1;
    for gY = gmin:grid:gmax
        count = 1;
        for j = 1:length(xR)
            if xR(j,1)>gX && xR(j,1)<gX+1
                if yR(j,1)>gY && yR(j,1)<gY+1
                    angularDevList{gXcount,gYcount}(count,1) = (angularDeviation(j,1));
                    count = count+1;
                end
            end
        end
        gYcount = gYcount+1;
    end
    gXcount = gXcount+1;
end
x = gmin+0.5:grid:gmax+0.5;
y = gmin+0.0:grid:gmax+0.0;

% compute heat map values
angularDevMean = NaN(length(angularDevList));
for i = 1:size(angularDevList,1)
    for j = 1:size(angularDevList,2)
        if length(angularDevList{i,j}) > 500
            angularDevMean(i,j) = (mean(angularDevList{i,j}));
        end
    end
end

% plot and aesthetics
h = imagesc(x,y,angularDevMean);
set(h,'alphadata',~isnan(angularDevMean));
colormap(jet);
colorbar;
xL = xlabel('        Left-Right (m)');
yL = ylabel('     Back-Front (m)');
axis([-6 6 -6 6]);
set(gca, ...
    'XTick',-6:3:6, ...
    'YTick',-6:3:6, ...
    'FontName',fn, ...
    'FontSize',16, ...
    'TickDir','out', ...
    'TickLength', [.02 .02]);
set([xL,yL], ...
    'FontName', fn, ...
    'FontSize', 24);
set(gcf, 'PaperPositionMode', 'auto');
% print -depsc2 -painters plotHeatMapHeading.eps


%%
% ========================================================================
% HEAT MAP #2: OCCUPANCY
% ========================================================================

count = 1;
totalFrames = 0;
for iTrial = 1:nTrial
    nHelm = max(Quant(iTrial).Nped);
    nFrame = length(Traj(iTrial,1).x);
    totalFrames = (totalFrames+nFrame*nHelm);
    
    % for all non-center helmets, compute position (x,y) and difference in heading
    for jFrame = 1:nFrame-1
        for kHelm = 1:nHelm
            cH = Quant(iTrial).ctrHelm(jFrame);
            if kHelm ~= cH;
                if isnan(Traj(iTrial,kHelm).hdn(jFrame,1)) == 0
                    
                    % difference in heading
                    angularDeviation(count,1) = (angleBn(Traj(iTrial,kHelm).hdn(jFrame,1),...
                        Traj(iTrial,cH).hdn(jFrame,1)));
                                        
                    % position (x,y) in cH-centered coordinates
                    x = Traj(iTrial,kHelm).x(jFrame)-Traj(iTrial,cH).x(jFrame);
                    y = Traj(iTrial,kHelm).y(jFrame)-Traj(iTrial,cH).y(jFrame);
                    th = Traj(iTrial,cH).hdn(jFrame);
                    xR(count,1) = x*cos(th)-y*sin(th);
                    yR(count,1) = x*sin(th)+y*cos(th);
                    dR(count,1) = sqrt(xR(count,1)^2+yR(count,1)^2);
                    
                    % take absolute value of deviation
                    angularDeviation(count,1) = abs(angularDeviation(count,1))*(180/pi);
                    
                    count = count+1;
                end
            end
        end
    end
end

% generate heat map data
grid = 0.2;   % size of grid squares (in m)
gmax = ceil(max(max(xR),max(yR)))+1;  % smallest value
gmin = -gmax;  % largest value
angularDevList = cell(2*gmax);
count = 1;
gXcount = 1; 
for gX = gmin:grid:gmax-1
     gYcount = 1;
    for gY = gmin:grid:gmax-1
        count = 1;
        for j = 1:length(xR)
            if xR(j,1)>gX && xR(j,1)<gX+1
                if yR(j,1)>gY && yR(j,1)<gY+1
                    angularDevList{gXcount,gYcount}(count,1) = (angularDeviation(j,1));
                    count = count+1;
                end
            end
        end
        gYcount = gYcount+1;
    end
    gXcount = gXcount+1;
end
x = gmin:grid:gmax;
y = gmin:grid:gmax;

% compute heat map values
angularDevCount = NaN(length(angularDevList));
for i = 1:size(angularDevList,1)
    for j = 1:size(angularDevList,2)
        if length(angularDevList{i,j}) > 1
            angularDevCount(i,j) = (length(angularDevList{i,j})/totalFrames)*100;
        end
    end
end

% plot and aesthetics
h = imagesc(x,y,angularDevCount);
% set(h,'alphadata',~isnan(angularDevCount));
colormap(jet);
colorbar;
xL = xlabel('        Left-Right (m)');
yL = ylabel('     Back-Front (m)');
axis([-6 6 -6 6]);
set(gca, ...
    'XTick',-6:3:6, ...
    'YTick',-6:3:6, ...
    'FontName',fn, ...
    'FontSize',16, ...
    'TickDir','out', ...
    'TickLength', [.02 .02]);
set([xL,yL], ...
    'FontName', fn, ...
    'FontSize', 24);
set(gcf, 'PaperPositionMode', 'auto');
% print -depsc2 -painters plotHeatMapOccupancy.eps