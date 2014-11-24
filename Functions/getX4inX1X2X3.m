function X4 = getX4inX1X2X3( helmetNb, X1, X2, X3 )
%% Return X4 coordinates in the coordinate system of {X1, X2, X3}
% Script uses the file X4inX1X2X3.csv to extract the coordinate of the head
% tracker in the system coordinate of the 3 trackers which identifiers are
% given in parameters.
%
% in param helmetNb: helmet number as defined in the Qualisys data files.
% in param X1, X2, X3: identifiers of the trackers defining the coordinate
%    system in which the out param X4 coordinates are expressed. Values are
%    e.g. 'a', 'd', 'e'.
%
% author: Stephane Bonneaud, 2012
% VENLab, swarm experiment 2012
% 

  data=importdata('X4inX1X2X3.csv');
  data.textdata=data.textdata(2:end,:);
  result=data.data(strcmp(data.textdata(:,1),num2str(helmetNb)), :);
  result=result(strcmp(data.textdata(strcmp(data.textdata(:,1),num2str(helmetNb)),2),[X1 X2 X3]),:);  
  X4=result(1:3);
end