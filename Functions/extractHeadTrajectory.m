function [estimates] = extractHeadTrajectory( helmetNb, trackers )
%function [headTrajectory,flag,sd] = extractHeadTrajectory( helmetNb, trackers )
%% EXTRACT HEAD TRAJECTORY GIVEN A HELMET AND 5 TRACKER TRAJECTORIES
% Script to use to extract the location in time of the head of participants
% from the locations in time of the 5 trackers of the helmet used for the
% swarm experiment 2012.
%
% in param helmetNb: the helmet's number as defined in Qualisys data files.
% in param trackers: timeseries of the locations of the 5 trackers for the
%    specific helmet.
%
% out param estimates: returns the estimated trajectories of the head based
% on the different set of 3 trackers (out of the 5), each set of 3 defining
% a coordinate system for the head, hence giving an estimation of the head
% trajectory.
%
% author: Stéphane Bonneaud, 2012
% VENLab, swarm experiment 2012

  % (1) checking trackers contains the timeseries (x,y,z) for the 5
  % trackers.
  if length(fieldnames(trackers)) ~= 5
    error('myApp:argChk', 'Wrong number of columns in input matrix.');
  end
  
  % (2) extract head trajectory estimates.
  combinations=combnk('abcde', 3);    
  estimates=zeros(length(trackers.a(1,:)),3,length(combinations));
  for t = 1:length(combinations)
    X1 = trackers.(combinations(t,1));
    X2 = trackers.(combinations(t,2));
    X3 = trackers.(combinations(t,3));
    
    % Gets the coordinates of X4 in system (X1, X2, X3) using the data file
    % containing all the coordinates for each helmet.
    X4 = getX4inX1X2X3( helmetNb, ...
      combinations(t,1), combinations(t,2), combinations(t,3) );

    for dt=1:length(X1(1,:))
      estimates(dt,:,t) = lc2gb( X1(1:3,dt), X2(1:3,dt), X3(1:3,dt), X4 );
  %      transfMat(X1(1:3,dt), X2(1:3,dt), X3(1:3,dt), X4);
    end
  end
  
  % (3) extract final head trajectory based on estimates and get
  % information on the quality of the data.
  %[Hpos, Sdev, Serr, Flag] = avg_hpos( estimates, 1 );
  
  % (4) extract compact information about the quality of the extracted
  % head trajectory, i.e. 4 percentages and the average of trackers used to
  % calculate the head trajectory.
  %[quality, avgTrackers] = getQuality( estimates, Flag );
      
end

function trackers = fillInInput( trajectoriesArray )
%%
%
  
  max = length(trajectoriesArray{1});
  for i = 2:length(trajectoriesArray)
    max
  end
end


function [quality, avgTrackers] = getQuality( estimates, flag )
%%
%

end
