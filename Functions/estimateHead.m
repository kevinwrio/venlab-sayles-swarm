function [Hpos] = estimateHead(hpos_estimates)


%% AVERAGE HEAD POSITION
  % THIS FUNCTION CALCULATES THE AVERAGE X, Y & Z POSITIONS FOR EACH TIMESTEP
  % OF ONE SAYLES SWARM TRIAL.
  
  %INOUTS:
    % HPOS_ESTIMATES = 3D MATRIX (TIME, (COL: X,Y,Z), ESTIMATE)
    % THRESHOLD = VALUE OF STANDARD ERROR THRESHOLD FOR UNRELIABLE
    %             ESTIMATE. IF STD ERROR IS > THRESHOLD, FLAG(I) = 1.
  %OUTPUTS:
    % HPOS = MATRIX OF AVG HEAD POSITION ESTIMATES (COL 1 - X, 2 - Y, 3 - Z) 
    % SDEV = STD DEV OF  AT EACH TIMESTEP
    % SDEVTOT = STD DEV ACROSS ENTIRE TRIAL
    % FLAG = PROBLEMATIC TRIALS (1 PROBLEMATIC, 0 OK)


%% CALCULATE AVERAGE HEAD POSITION

for i = 1:length(hpos_estimates(:,1,1));
    
    % Average x,y,z position
    Hpos(i,1) = nanmean(hpos_estimates(i,1,:));
    Hpos(i,2) = nanmean(hpos_estimates(i,2,:));
    Hpos(i,3) = nanmean(hpos_estimates(i,3,:));
    
    % SD of estimates (x-position)
    % SD(i) = std(hpos_estimates(i,1,:));
    
    % if SD > threshold, eliminate estimates until SD < threshold
    % threshold = 10; % cm
    % while SD(i) > threshold
        
        
        
    
    
end
    
end
