
%% IDENTIFY THE CRITICAL SUBSTATIONS

clear

%add paths
addpath('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion');
%addpath('C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\MatlabResilienceModel\\ThesisVersion');

addpath(genpath('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/matpower4.0b4/'));
%addpath(genpath('C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\MatlabResilienceModel\\ThesisVersion\\matpower4.0b4'));

addpath('/home/andrewbollinger/AndrewBollinger/projects/NetlogoModels/NetworkEvolution/modelFiles/outputdata/');
%addpath('C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\NetlogoModels\\NetworkEvolution\\modelFiles\\outputdata')

%set options
printstuff = false;
allowcascadingfailures = true;
allowgeneratorredispatching = true;
redispatchattemptslimit = 100;
mpopt = mpoption('PF_DC', 1, 'OUT_ALL', 0, 'VERBOSE', 0);
iterations = 1000;
plotstuff = 1;

% create a matrix to hold the results
LoadNetwork
networkperformancematrix = [];
performancedropmatrix = [];

i = 1;
while i <= iterations

    % load the network data
    LoadNetwork

    % create some variables to hold the results
    busesremoved = zeros(length(mpc.bus(:,1)),1);
    performancedrop = zeros(length(mpc.bus(:,1)),1);
    networkperformance = zeros(length(mpc.bus(:,1)),1);
    
    disp(i);
    
    % set the bus demand values
    mindemandvector = mpc.bus(:,17);
    maxdemandvector = mpc.bus(:,16);
    demandvector = mindemandvector + (maxdemandvector - mindemandvector) * rand(1);
    mpc.bus(:,3) = demandvector;

    % create a random generation dispatch distribution
    generatoroutputs = random('uniform', 0, 1, length(mpc.gen(:,1)), 1);
    generatoroutputs = generatoroutputs .* mpc.gen(:,9);
    generatoroutputs = generatoroutputs * sum(mpc.bus(:,3)) / sum(generatoroutputs);
    mpc.gen(:,2) = generatoroutputs;
    
    % successively remove buses and calculate network performance
    numberofbusestoremove = length(mpc.bus(:,1));
    busestoremove = mpc.bus;
    lostloadfraction = 0;
    iteration = 1;
    for numberofbusesremoved = 1:numberofbusestoremove
        
        % get the number of branches remaining
        numrows = length(mpc.branch(:,1));
        
        redispatchattempts = 0;
        
        % remove a bus
        x = randsample(1:length(busestoremove(:,1)),1,true);
        y = busestoremove(x,1);
        mpc.branch((mpc.branch(:,1) == y),:) = [];
        mpc.branch((mpc.branch(:,2) == y),:) = [];
        mpc.gen((mpc.gen(:,1) == y),9) = 0;
        mpc.gen((mpc.gen(:,1) == y),2) = 0;
        busestoremove(x,:) = [];
        busesremoved(iteration) = y;

        % calculate the power flows
        CalculatePowerFlows;

        % calculate the network performance
        lostload = sum(mpc.bus(:,3)) - sum(busresults2(:,2));
        lostloadfraction = lostload ./ sum(mpc.bus(:,3));
        networkperformance(iteration) = 1 - lostloadfraction;
        
        % calculate the performance drop
        if iteration == 1
            performancedrop(iteration) = 1 - networkperformance(iteration);
        else
            performancedrop(iteration) = networkperformance(iteration-1) - networkperformance(iteration);
        end
        
        iteration = iteration + 1;
    end
    
    % plot the network performance (optional)
    if plotstuff == 1
        plot(networkperformance,'LineWidth',1,'Color',[0.7 0.7 0.7]);
        hold on
    end

    % update the network performance matrix
    networkperformancematrix = horzcat(networkperformancematrix, networkperformance);
    
    % update the performance drop matrix
    performancedrop = horzcat(busesremoved, performancedrop);
    performancedrop = sortrows(performancedrop, 1);
    performancedrop = performancedrop(:,2);
    performancedropmatrix = horzcat(performancedropmatrix, performancedrop);
    
    % update the iterator
    i = i + 1;

end

% plot mean network performance
networkperformancematrix = networkperformancematrix';
meanperformance = mean(networkperformancematrix);
plot(meanperformance, 'LineWidth',3,'Color','black');
xlabel('Number of substations failed')
ylabel('Network performance')
ylim([0 1.05])
hold off

% calculate the criticality of substations
substationcriticality = mean(performancedropmatrix,2);
csvwrite('SubstationCriticality.csv',substationcriticality);

substationcriticality_sorted = horzcat(mpc.bus(:,1), substationcriticality);
substationcriticality_sorted = sortrows(substationcriticality_sorted,-2);
