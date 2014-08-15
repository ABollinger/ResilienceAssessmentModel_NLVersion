% THIS CODE SAVES THE NECESSARY MATRICES TO THE CORRECT NAME
% MATRICES ARE NAMED ACCORDING TO THE SCENARIO

% % save the networkperformancematrix
if fractionofvulnerablebusestoprotect == 0
    networkperformancematrix00 = networkperformancematrix;
elseif fractionofvulnerablebusestoprotect == 0.2
    networkperformancematrix20 = networkperformancematrix;
elseif fractionofvulnerablebusestoprotect == 0.4
    networkperformancematrix40 = networkperformancematrix;
elseif fractionofvulnerablebusestoprotect == 0.6
    networkperformancematrix60 = networkperformancematrix;
elseif fractionofvulnerablebusestoprotect == 0.8
    networkperformancematrix80 = networkperformancematrix;
elseif fractionofvulnerablebusestoprotect == 1
    networkperformancematrix100 = networkperformancematrix;
else
    disp('ERROR: CANNOT SAVE NETWORKPERFORMANCEMATRIX')
end

% ALTERNATIVE CODE
% if fractionofvulnerablebusestoprotect == 0 && flooddefensescenario == 19 
%     networkperformancematrix00a = networkperformancematrix;
%     eventprobabilitymatrix00a = eventprobabilitymatrix;
% 
% elseif fractionofvulnerablebusestoprotect == 0 && flooddefensescenario == 20 
%     networkperformancematrix00b = networkperformancematrix;
%     eventprobabilitymatrix00b = eventprobabilitymatrix;
%         
% elseif fractionofvulnerablebusestoprotect == 0 && flooddefensescenario == 21 
%     networkperformancematrix00c = networkperformancematrix;
%     eventprobabilitymatrix00c = eventprobabilitymatrix;
% 
% elseif fractionofvulnerablebusestoprotect == 0 && flooddefensescenario == 22 
%     networkperformancematrix00d = networkperformancematrix;
%     eventprobabilitymatrix00d = eventprobabilitymatrix;
% 
% 
% elseif fractionofvulnerablebusestoprotect == 0.2 && flooddefensescenario == 19 
%     networkperformancematrix20a = networkperformancematrix;
%     eventprobabilitymatrix20a = eventprobabilitymatrix;
% 
% elseif fractionofvulnerablebusestoprotect == 0.2 && flooddefensescenario == 20 
%     networkperformancematrix20b = networkperformancematrix;
%     eventprobabilitymatrix20b = eventprobabilitymatrix;
%         
% elseif fractionofvulnerablebusestoprotect == 0.2 && flooddefensescenario == 21 
%     networkperformancematrix20c = networkperformancematrix;
%     eventprobabilitymatrix20c = eventprobabilitymatrix;
% 
% elseif fractionofvulnerablebusestoprotect == 0.2 && flooddefensescenario == 22 
%     networkperformancematrix20d = networkperformancematrix;
%     eventprobabilitymatrix20d = eventprobabilitymatrix;
% 
% 
% elseif fractionofvulnerablebusestoprotect == 0.4 && flooddefensescenario == 19 
%     networkperformancematrix40a = networkperformancematrix;
%     eventprobabilitymatrix40a = eventprobabilitymatrix;
% 
% elseif fractionofvulnerablebusestoprotect == 0.4 && flooddefensescenario == 20 
%     networkperformancematrix40b = networkperformancematrix;
%     eventprobabilitymatrix40b = eventprobabilitymatrix;
%         
% elseif fractionofvulnerablebusestoprotect == 0.4 && flooddefensescenario == 21 
%     networkperformancematrix40c = networkperformancematrix;
%     eventprobabilitymatrix40c = eventprobabilitymatrix;
% 
% elseif fractionofvulnerablebusestoprotect == 0.4 && flooddefensescenario == 22 
%     networkperformancematrix40d = networkperformancematrix;
%     eventprobabilitymatrix40d = eventprobabilitymatrix;
% 
% 
% elseif fractionofvulnerablebusestoprotect == 0.6 && flooddefensescenario == 19 
%     networkperformancematrix60a = networkperformancematrix;
%     eventprobabilitymatrix60a = eventprobabilitymatrix;
% 
% elseif fractionofvulnerablebusestoprotect == 0.6 && flooddefensescenario == 20 
%     networkperformancematrix60b = networkperformancematrix;
%     eventprobabilitymatrix60b = eventprobabilitymatrix;
%         
% elseif fractionofvulnerablebusestoprotect == 0.6 && flooddefensescenario == 21 
%     networkperformancematrix60c = networkperformancematrix;
%     eventprobabilitymatrix60c = eventprobabilitymatrix;
% 
% elseif fractionofvulnerablebusestoprotect == 0.6 && flooddefensescenario == 22 
%     networkperformancematrix60d = networkperformancematrix;
%     eventprobabilitymatrix60d = eventprobabilitymatrix;
% 
% 
% elseif fractionofvulnerablebusestoprotect == 0.8 && flooddefensescenario == 19 
%     networkperformancematrix80a = networkperformancematrix;
%     eventprobabilitymatrix80a = eventprobabilitymatrix;
% 
% elseif fractionofvulnerablebusestoprotect == 0.8 && flooddefensescenario == 20 
%     networkperformancematrix80b = networkperformancematrix;
%     eventprobabilitymatrix80b = eventprobabilitymatrix;
%         
% elseif fractionofvulnerablebusestoprotect == 0.8 && flooddefensescenario == 21 
%     networkperformancematrix80c = networkperformancematrix;
%     eventprobabilitymatrix80c = eventprobabilitymatrix;
% 
% elseif fractionofvulnerablebusestoprotect == 0.8 && flooddefensescenario == 22 
%     networkperformancematrix80d = networkperformancematrix;
%     eventprobabilitymatrix80d = eventprobabilitymatrix;
% 
% 
% elseif fractionofvulnerablebusestoprotect == 1 && flooddefensescenario == 19 
%     networkperformancematrix100a = networkperformancematrix;
%     eventprobabilitymatrix100a = eventprobabilitymatrix;
% 
% elseif fractionofvulnerablebusestoprotect == 1 && flooddefensescenario == 20 
%     networkperformancematrix100b = networkperformancematrix;
%     eventprobabilitymatrix100b = eventprobabilitymatrix;
%         
% elseif fractionofvulnerablebusestoprotect == 1 && flooddefensescenario == 21 
%     networkperformancematrix100c = networkperformancematrix;
%     eventprobabilitymatrix100c = eventprobabilitymatrix;
% 
% elseif fractionofvulnerablebusestoprotect == 1 && flooddefensescenario == 22 
%     networkperformancematrix100d = networkperformancematrix;
%     eventprobabilitymatrix100d = eventprobabilitymatrix;
% 
% else
%     disp('ERROR: CANNOT SAVE NETWORKPERFORMANCEMATRIX')
% 
% end

%OLD CODE
% % save the networkperformancematrix
% if fractionofvulnerablebusestoprotect == 0
%     networkperformancematrix00 = networkperformancematrix;
% elseif fractionofvulnerablebusestoprotect == 0.2
%     networkperformancematrix20 = networkperformancematrix;
% elseif fractionofvulnerablebusestoprotect == 0.4
%     networkperformancematrix40 = networkperformancematrix;
% elseif fractionofvulnerablebusestoprotect == 0.6
%     networkperformancematrix60 = networkperformancematrix;
% elseif fractionofvulnerablebusestoprotect == 0.8
%     networkperformancematrix80 = networkperformancematrix;
% elseif fractionofvulnerablebusestoprotect == 1
%     networkperformancematrix100 = networkperformancematrix;
% else
%     disp('ERROR: CANNOT SAVE NETWORKPERFORMANCEMATRIX')
% end
% 
% % save the eventprobabilitymatrix
% if fractionofvulnerablebusestoprotect == 0
%     eventprobabilitymatrix00 = eventprobabilitymatrix;
% elseif fractionofvulnerablebusestoprotect == 0.2
%     eventprobabilitymatrix20 = eventprobabilitymatrix;
% elseif fractionofvulnerablebusestoprotect == 0.4
%     eventprobabilitymatrix40 = eventprobabilitymatrix;
% elseif fractionofvulnerablebusestoprotect == 0.6
%     eventprobabilitymatrix60 = eventprobabilitymatrix;
% elseif fractionofvulnerablebusestoprotect == 0.8
%     eventprobabilitymatrix80 = eventprobabilitymatrix;
% elseif fractionofvulnerablebusestoprotect == 1
%     eventprobabilitymatrix100 = eventprobabilitymatrix;
% else
%     disp('ERROR: CANNOT SAVE EVENTPROBABILITYMATRIX')
% end