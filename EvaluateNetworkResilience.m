
%% TODO

% Assumption: We are ignoring the flood risk of substations in the
% buitendijksgebied I think.

% Verify the CalculatePowerFlows code

% Confirm your assumption for the case of heat waves about ignoring the
% radiation components w/ relation to line capacities.

% Compile risk data for the different types of events and integrate that
% into the code

% Determine adaptations to test and implement adaptations code

% Incorporate code for calculating movement through the state space

% Calculate resilience values, also max and min values, for each case, both
% with the risk factors and without the risk factors.

% Validate the model by comparing a default case with line capacities and
% maybe with load flow data.  If there's not time for this, you can leave
% this blank and do it later.

% OPTIONAL:

% Include power plant maintenance scenarios based on planned unavailability
% data?

% Fill in missing values in distribution grid demand data based on the
% known total peak consumption in 2010?

% Include sudden failures due to line sag/soil movement?

%% SET GLOBAL VARIABLES

%add paths
%addpath('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion');
addpath('C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\MatlabResilienceModel\\ThesisVersion');

%addpath(genpath('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/matpower4.0b4/'));
addpath(genpath('C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\MatlabResilienceModel\\ThesisVersion\\matpower4.0b4'));

%addpath('/home/andrewbollinger/AndrewBollinger/projects/NetlogoModels/NetworkEvolution/modelFiles/outputdata/');
addpath('C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\NetlogoModels\\NetworkEvolution\\modelFiles\\outputdata')

%set options
printstuff = false;
allowcascadingfailures = true;
mpopt = mpoption('PF_DC', 1, 'OUT_ALL', 0, 'VERBOSE', 0);
iterations = 2;
plotnetworkperformance = 1;

%% IMPLEMENT ADAPTATIONS

    
%% EVALUATE WINDSTORM RESILIENCE

% create a matrix to hold the results
LoadNetwork
networkperformancematrix = zeros(length(mpc.branch(:,1)), iterations);

i = 1;
while i <= iterations

    % load the network data
    LoadNetwork

    % create some variables to hold the results
    networkperformance = [];
    lostloadfraction = 0;
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

    % set the wind vulnerability values
    branchlength = mpc.branch(:,17);
    branchlength_normalized = branchlength / max(branchlength);
    branchfractionunderground = mpc.branch(:,18);
    branchwindexposure = mpc.branch(:,19);
    branchwindexposure_normalized = branchwindexposure / max(branchwindexposure);
    branchwindstormvulnerability = branchlength_normalized .* (1 - branchfractionunderground) .* branchwindexposure_normalized;
    
    % add the wind vulnerability vector to the branch matrix
    mpc.branch = horzcat(mpc.branch, branchwindstormvulnerability);

    % successively remove branches and calculate network performance
    numberofbranchestoremove = length(mpc.branch(:,1));
    branchestoremove = mpc.branch;
    for numberofbranchesremoved = 1:numberofbranchestoremove

        % get the number of branches remaining
        numrows = length(mpc.branch(:,1));

        % calculate the power flows
        CalculatePowerFlows;

        % calculate network performance
        lostload = sum(mpc.bus(:,3)) - sum(busresults2(:,2));
        lostloadfraction = lostload ./ sum(mpc.bus(:,3));
        networkperformance = vertcat(networkperformance, 1 - lostloadfraction);
        
        % determine the branch to be removed based on the vulnerability and
        % remove the branch, if it has not already failed;
        if sum(branchestoremove(:,20)) > 0
            x = randsample(1:length(branchestoremove(:,1)),1,true,branchestoremove(:,20));
            branchtoremove_frombusnumber = branchestoremove(x,1);
            branchtoremove_tobusnumber = branchestoremove(x,2);
            branchtoremove_index = find(mpc.branch(:,1) == branchtoremove_frombusnumber & mpc.branch(:,2) == branchtoremove_tobusnumber);
            if ~isempty(branchtoremove_index)
                mpc.branch(branchtoremove_index,:) = [];
            end
            branchestoremove(x,:) = [];
        end
    end
    
    % plot the network performance (optional)
    if plotnetworkperformance == 1
        plot(networkperformance,'LineWidth',1,'Color',[0.7 0.7 0.7]);
        hold on
    end

    % update the network performance matrix
    for z = 1:length(networkperformance(:,1))
        networkperformancematrix(z,i) = networkperformance(z,1);
    end
    
    % update the iterator
    i = i + 1;
end

% plot mean network performance
networkperformancematrix = networkperformancematrix';
meanperformance = mean(networkperformancematrix);

plot(meanperformance, 'LineWidth',3,'Color','black');
ylim([0 1.05])
xlabel('Number of lines failed')
ylabel('Network performance')
saveas(gcf,'output/windstormperformance.fig')
saveas(gcf,'output/windstormperformance.png')
hold off

%% EVALUATE FLOOD RESILIENCE

% create a matrix to hold the results
LoadNetwork
networkperformancematrix = zeros(length(mpc.bus(:,1)), iterations);

i = 1;
while i <= iterations

    % load the network data
    LoadNetwork

    % create some variables to hold the results
    lostloadfraction = 0;
    networkperformance = [];
    disp(i);
    
    % set the bus demand values
    mindemandvector = mpc.bus(:,17);
    maxdemandvector = mpc.bus(:,16);
    demandvector = mindemandvector + (maxdemandvector - mindemandvector) * rand(1);
    mpc.bus(:,3) = demandvector;
    
%     peakpowerdemand = sum(mpc.bus(:,3));
%     electricitydemand2010 = 109351000; %(in MWh, from EnergieInNederland2011 Report)
%     powerdemand2010 = electricitydemand2010 / 8760;
%     demandfactor = powerdemand2010 / peakpowerdemand;
%     meandemandvector = mpc.bus(:,3) * demandfactor;
%     minpowerdemand = powerdemand2010 - (peakpowerdemand - powerdemand2010);
%     randomfactor = ((peakpowerdemand - minpowerdemand) * rand(1) + minpowerdemand) / powerdemand2010;
%     mpc.bus(:,3) = meandemandvector * randomfactor;

    % create a random generation dispatch distribution
    generatoroutputs = random('uniform', 0, 1, length(mpc.gen(:,1)), 1);
    generatoroutputs = generatoroutputs .* mpc.gen(:,9);
    generatoroutputs = generatoroutputs * sum(mpc.bus(:,3)) / sum(generatoroutputs);
    mpc.gen(:,2) = generatoroutputs;

    % create the weight matrix, which will be used to determine the removal of substations
    % The imported vulnerability values vary between 0 and 5, which relate
    % to different maximum flood height levels, as follows: 0 = 0m, 1 =
    % 0-0.2m, 2 = 0.2-0.5m, 3 = 0.5 - 0.8m, 4 = 0.8 - 2m, 5 = 2 - 5m. We
    % assume that substations are protected up to 0.8m (the supposed actual
    % value is 1m), so all values less than 3 are considered to constitute
    % a risk of 0.
    busfloodvulnerabilities = mpc.bus(:,15);
    busfloodvulnerabilities = busfloodvulnerabilities - 3;
    busfloodvulnerabilities_normalized = busfloodvulnerabilities / (5 - 3);
    mpc.bus = horzcat(mpc.bus, busfloodvulnerabilities_normalized);
    
    % successively remove buses and calculate network performance
    numberofbusestoremove = length(mpc.bus(:,1));
    busestoremove = mpc.bus;
    for numberofbusesremoved = 1:numberofbusestoremove
        
        % get the number of branches remaining
        numrows = length(mpc.branch(:,1));

        % calculate the power flows
        CalculatePowerFlows;

        % calculate the network performance
        lostload = sum(mpc.bus(:,3)) - sum(busresults2(:,2));
        lostloadfraction = lostload ./ sum(mpc.bus(:,3));
        networkperformance = vertcat(networkperformance, 1 - lostloadfraction);
        
        % determine the bus to be removed based on the vulnerability and remove the bus
        if sum(busestoremove(:,15)) > 0
            x = randsample(1:length(busestoremove(:,1)),1,true,busestoremove(:,15));
            y = busestoremove(x,1);
            mpc.branch((mpc.branch(:,1) == y),:) = [];
            mpc.branch((mpc.branch(:,2) == y),:) = [];
            mpc.gen((mpc.gen(:,1) == y),9) = 0;
            mpc.gen((mpc.gen(:,1) == y),2) = 0;
            busestoremove(x,:) = [];
        end
    end
    
    % plot the network performance (optional)
    if plotnetworkperformance == 1
        plot(networkperformance,'LineWidth',1,'Color',[0.7 0.7 0.7]);
        hold on
    end

    % update the network performance matrix
    for z = 1:length(networkperformance(:,1))
        networkperformancematrix(z,i) = networkperformance(z,1);
    end
    
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
saveas(gcf,'output/floodperformance.fig')
saveas(gcf,'output/floodperformance.png')

hold off

%% EVALUATE HEAT WAVE RESILIENCE

% create a matrix to hold the results
LoadNetwork

% set the temperature levels
mintemperature = 30; %for a heat wave the peak temperature has to be above 30 deg C
maxtemperature = 41;
meantemperaturesummer2012 = 15.7; %mean daily temp between May and Sept 2012

% create the results matrix
networkperformancematrix = zeros(maxtemperature - mintemperature + 1, iterations);

i = 1;
while i <= iterations

    % load the network data
    LoadNetwork
    fractionunderground = mpc.branch(:,18);

    % create some variables to hold the results
    lostloadfraction = 0;
    networkperformance = [];
    disp(i);
    
    % set the bus demand values
    mindemandvector = mpc.bus(:,17);
    maxdemandvector = mpc.bus(:,16);
    %adjust demand vector to account for lower demand during summer
    %14483 is the maximum May-Sep 2013 demand
    %16101 is the maximum Apr-Dec 2013 demand
    maxdemandvector = maxdemandvector * 14483 / 16101;
    demandvector = mindemandvector + (maxdemandvector - mindemandvector) * rand(1);
    mpc.bus(:,3) = demandvector;
    
    %peakpowerdemand = sum(mpc.bus(:,3));
    %meanpowerdemandsummer2012 = 13616.27; %mean maximum daily power demand between May and Sept 2012, excluding weekends
    %demandfactor = meanpowerdemandsummer2012 / peakpowerdemand;
    %mpc.bus(:,3) = mpc.bus(:,3) * demandfactor;
    %buslist_original = mpc.bus;

    % create a random generation dispatch distribution
    generatoroutputs = random('uniform', 0, 1, length(mpc.gen(:,1)), 1);
    generatoroutputs = generatoroutputs .* mpc.gen(:,9);
    generatoroutputs = generatoroutputs * sum(mpc.bus(:,3)) / sum(generatoroutputs);
    mpc.gen(:,2) = generatoroutputs;
    
    % calculate generator temperature sensitivity factors
    generatorssensitivecapacities = mpc.gen(:,24);
    generatortemperaturesensitivityfactor = zeros(length(mpc.gen(:,1)),1);
    generatorswithnonzerocapacity_indices = find(mpc.gen(:,9));
    generatortemperaturesensitivityfactor(generatorswithnonzerocapacity_indices) = ...
        generatorssensitivecapacities(generatorswithnonzerocapacity_indices) ./ mpc.gen(generatorswithnonzerocapacity_indices,9);

    % create some vectors for use later
    genlist_original = mpc.gen;
    buslist_original = mpc.bus;
    capacitylist_original = capacitylist;
    
    % successively increase the temperature
    for temperature = mintemperature:maxtemperature
        
        % get the number of branches
        numrows = length(mpc.branch(:,1));
        
        % modify line capacities depending on the temperature
        Ta = temperature; %ambient temp
        Taref = 25; %reference temp (rating temp)
        Tcmax = 100; %maximum conductor temp
        capacitylist(:,3) = capacitylist_original(:,3) * sqrt((Tcmax - Ta) / (Tcmax - Taref));
        
        % modify demand magnitudes depending on the temperature
        Ta = temperature; %ambient temp
        Taref = 20.15; %reference temp = mean maximum daily temp at De Bilt between 1/5/2010 and 31/09/2010
        mpc.bus(:,3) = buslist_original(:,3) * (1 + 0.005 * (Ta - Taref));
        
        % modify generation capacity and output depending on the temperature
        Ta = temperature; %ambient temp
        Taref = 25; %reference temp (assumed)
        mpc.gen(:,9) = genlist_original(:,9) .* (1 - 0.006 * generatortemperaturesensitivityfactor * (Ta - Taref));
        
        % modify generation availability depending on the temperature
        % assumption: there are lots of shaky assumptions here
        %coolingwaterrisk = (temperature - 35) * 0.1;
        %mpc.gen(:,9) = mpc.gen(:,9) - generatorssensitivecapacities * coolingwaterrisk;
        
        % change the output of generators whose capacity is output exceeds their capacity
        %generatorswithoutputexceedingcapacity_indices = find(mpc.gen(:,9) < mpc.gen(:,2),2);
        %mpc.gen(generatorswithoutputexceedingcapacity_indices,2) = mpc.gen(generatorswithoutputexceedingcapacity_indices,9); 

        % calculate the power flows
        CalculatePowerFlows;

        % calculate the network performance
        lostload = sum(mpc.bus(:,3)) - sum(busresults2(:,2));
        lostloadfraction = lostload ./ sum(mpc.bus(:,3));
        networkperformance = vertcat(networkperformance, 1 - lostloadfraction);

    end
    
    % plot the network performance (optional)
    if plotnetworkperformance == 1
        plot(networkperformance,'LineWidth',1,'Color',[0.7 0.7 0.7]);
        hold on
    end

    % update the network performance matrix
    for z = 1:length(networkperformance(:,1))
        networkperformancematrix(z,i) = networkperformance(z,1);
    end
    
    % update the iterator
    i = i + 1;

end

% plot mean network performance
networkperformancematrix = networkperformancematrix';
meanperformance = mean(networkperformancematrix);

plot(meanperformance, 'LineWidth',3,'Color','black');
xlabel('Event magnitude (temperature)')
ylabel('Network performance')
ylim([0 1.05])
saveas(gcf,'output/heatwaveperformance.fig')
saveas(gcf,'output/heatwaveperformance.png')

hold off

%% PLOT THE COMBINED RESULTS

plot(scenario1_windstorm_meanperformance,'LineWidth',2,'Color','red');
hold on
plot(scenario2_windstorm_meanperformance,'LineWidth',2,'Color','green');
plot(scenario3_windstorm_meanperformance,'LineWidth',2,'Color','blue');
plot(scenario4_windstorm_meanperformance,'LineWidth',2,'Color','cyan');
plot(scenario5_windstorm_meanperformance,'LineWidth',2,'Color','magenta');
hold off
ylim([0 1.05])
title('Windstorm vulnerability of the infrastructure','FontSize',18);
ylabel('Network performance','Fontsize',16);
xlabel('Number of lines failed','Fontsize',16);
leg = legend('Baseline','Centralized','Distributed','Offshore wind','Import');
set(leg,'FontSize',14);
xtick =  get(gca,'XTickLabel');
set(gca,'XTickLabel',xtick,'FontSize',12)
ytick =  get(gca,'YTickLabel');
set(gca,'YTickLabel',ytick,'FontSize',12)
saveas(gcf,'output/windstormperformance.png')
saveas(gcf,'output/windstormperformance.fig')
saveas(gcf,'output/windstormperformance.png')


plot(scenario1_flood_meanperformance,'LineWidth',2,'Color','red');
hold on
plot(scenario2_flood_meanperformance,'LineWidth',2,'Color','green');
plot(scenario3_flood_meanperformance,'LineWidth',2,'Color','blue');
plot(scenario4_flood_meanperformance,'LineWidth',2,'Color','cyan');
plot(scenario5_flood_meanperformance,'LineWidth',2,'Color','magenta');
hold off
ylim([0 1.05])
title('Flood vulnerability of the infrastructure','FontSize',18);
ylabel('Network performance','Fontsize',16);
xlabel('Number of substations failed','Fontsize',16);
leg = legend('Baseline','Centralized','Distributed','Offshore wind','Import');
set(leg,'FontSize',14);
xtick =  get(gca,'XTickLabel');
set(gca,'XTickLabel',xtick,'FontSize',12)
ytick =  get(gca,'YTickLabel');
set(gca,'YTickLabel',ytick,'FontSize',12)
saveas(gcf,'output/floodperformance.fig')
saveas(gcf,'output/floodperformance.png')













