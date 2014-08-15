%% PLOT THE FLOOD RESULTS - LINE PLOTS

clear

% load the results
load('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/1000iterationsb/floodperformance.mat');
outputpath = ('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/finalplots_1000iterations/');

%load('C:\Users\Andrew\Documents\TUD\AndrewBollinger\projects\MatlabResilienceModel\ThesisVersion\output\1000iterationsb\floodperformance.mat')
%outputpath = ('C:\Users\Andrew\Documents\TUD\AndrewBollinger\projects\MatlabResilienceModel\ThesisVersion\output\finalplots_1000iterations\');

figure(1)
z = figure(1);

% plot the default case results
xvalues = 1:numberofbusestoremove;
xvalues = xvalues - 1;
subplot(1,2,1);
plot(xvalues, networkperformancematrix00,'LineWidth',1,'Color',[0.7 0.7 0.7])
hold on
plot(xvalues, mean(networkperformancematrix00'),'LineWidth',3,'Color','black')
hold off
title({'Infrastructure performance vs. Event magnitude';'(no adaptation measures)'},'FontSize',16)
xlabel('Event magnitude (# substations failed)','FontSize',15)
ylabel('Performance (fraction demand served)','FontSize',15)
xlim([0 max(xvalues)])
ylim([min(min(networkperformancematrix00)) - 0.05 1.01])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',14)
set(gca,'XTickMode','manual','YTickMode','manual')
saveas(gcf,strcat(outputpath,'floodperformance_noadaptations.fig'));
saveas(gcf,strcat(outputpath,'floodperformance_noadaptations.png'));

% plot the combined results
subplot(1,2,2);
plot(xvalues, meanperformancematrix(:,1),...
    xvalues, meanperformancematrix(:,2),...
    xvalues, meanperformancematrix(:,3),...
    xvalues, meanperformancematrix(:,4),...
    xvalues, meanperformancematrix(:,5),...
    xvalues, meanperformancematrix(:,6),'LineWidth',3);
title({'Infrastructure performance vs. Event magnitude';'(comparison of measures)'},'FontSize',16)
xlabel('Event magnitude (# substations failed)','FontSize',15)
ylabel('Performance (fraction demand served)','FontSize',15)
leg = legend('0% protected','20% protected','40% protected','60% protected','80% protected','100% protected','Location','SouthWest');
set(leg,'FontSize',12)
v = get(leg,'title');
set(v,'string','Adaptation measure','FontSize',14);
xlim([0 max(xvalues)])
ylim([min(mean(networkperformancematrix00)) - 0.01 1.01])
colormap('Lines')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',14)
set(gca,'XTickMode','manual','YTickMode','manual')
saveas(gcf,strcat(outputpath,'floodperformance_adaptationscomparison.fig'));
saveas(gcf,strcat(outputpath,'floodperformance_adaptationscomparison.png'));

set(z,'PaperUnits','inches','PaperPosition',[0 0 14 6])
%print -dpng floodresults.png;
print(z,'-dpng','/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/finalplots_1000iterations/floodresults.png');

min(min(networkperformancematrix00))
min(mean(networkperformancematrix00))
max(mean(networkperformancematrix00))
std(mean(networkperformancematrix00))

%% PLOT THE FLOOD RESULTS - SUMMARY BOX PLOT

clf('reset')
figure(2)
z = figure(2);
reset(z)

% box plot of combined results
networkperformancematrix_total = [mean(networkperformancematrix00);...
    mean(networkperformancematrix20);...
    mean(networkperformancematrix40);...
    mean(networkperformancematrix60);...
    mean(networkperformancematrix80);...
    mean(networkperformancematrix100)];
h = boxplot(networkperformancematrix_total','boxstyle','filled','symbol','k+',...
    'labels',{'0%','20%','40%','60%','80%','100%'});
set(h,'color','black')
title({'Flood resilience of the infrastructure under different degrees of flood protection';'(prioritization on the basis of substation vulnerability)'},'FontSize',16)
xlabel({'Adaptation measure';'(% vulnerable substations protected)'},'FontSize',15)
ylabel('Resilience','FontSize',15)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',14)
saveas(gcf,strcat(outputpath,'floodperformance_adaptationscomparison_box.fig'));
saveas(gcf,strcat(outputpath,'floodperformance_adaptationscomparison_box.png'));


%% UNUSED FLOOD RESULTS PLOTS 

% create a histogram and kernel density plot
hist(mean(networkperformancematrix00),20)
set(get(gca,'child'),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
hold on
ksdensity(mean(networkperformancematrix00))
set(get(gca,'child'),'Color','black','LineWidth',2);
hold off

% create a bar plot with error bars
addpath('C:\Users\Andrew\Documents\TUD\AndrewBollinger\projects\MatlabResilienceModel\ThesisVersion')

%   Symmetric Example:
y = [mean(mean(networkperformancematrix00));...
    mean(mean(networkperformancematrix20));...
    mean(mean(networkperformancematrix40));...
    mean(mean(networkperformancematrix60));...
    mean(mean(networkperformancematrix80));...
    mean(mean(networkperformancematrix100))];
errY = zeros(6,1,2);
errY(:,:,1) = [mean(mean(networkperformancematrix00)) - min(mean(networkperformancematrix00));...
    mean(mean(networkperformancematrix20)) - min(mean(networkperformancematrix20));...
    mean(mean(networkperformancematrix40)) - min(mean(networkperformancematrix40));...
    mean(mean(networkperformancematrix60)) - min(mean(networkperformancematrix60));...
    mean(mean(networkperformancematrix80)) - min(mean(networkperformancematrix80));...
    mean(mean(networkperformancematrix100)) - min(mean(networkperformancematrix100))];
errY(:,:,2) = [max(mean(networkperformancematrix00)) - mean(mean(networkperformancematrix00));...
    max(mean(networkperformancematrix20)) - mean(mean(networkperformancematrix20));...
    max(mean(networkperformancematrix40)) - mean(mean(networkperformancematrix40));...
    max(mean(networkperformancematrix60)) - mean(mean(networkperformancematrix60));...
    max(mean(networkperformancematrix80)) - mean(mean(networkperformancematrix80));...
    max(mean(networkperformancematrix100)) - mean(mean(networkperformancematrix100))];
h = barwitherr(errY, y);% Plot with errorbars
set(gca,'XTickLabel',{'0%','20%','40%','60%','80%','100%'})
ylim([min(mean(networkperformancematrix00)) - 0.1 1.01])
set(h(1),'FaceColor',[0.7 0.7 0.7]);
title('Flood resilience of the infrastructure (comparison of measures)','FontSize',17)
xlabel('Adaptation measure (% vulnerable substations protected)','FontSize',14)
ylabel('Resilience','FontSize',14)

%% PLOT THE HEAT WAVE RESULTS - LINE PLOTS

clear

% load the results
load('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/1000iterationsb/heatwaveperformance.mat')
outputpath = ('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/finalplots_1000iterations/');

%load('C:\Users\Andrew\Documents\TUD\AndrewBollinger\projects\MatlabResilienceModel\ThesisVersion\output\100iterations\heatwaveperformance_repaired.mat')
%outputpath = ('C:\Users\Andrew\Documents\TUD\AndrewBollinger\projects\MatlabResilienceModel\ThesisVersion\output\finalplots_1000iterations\');

figure(1)
z = figure(1);
axislabelfontsize = 10;
titlelabelfontsize = 12;

% plot the default case results
xincrement = capacitytosubtracteachiteration / 1000;
xvalues = xincrement:xincrement:length(meanperformance)*xincrement;
subplot(1,2,1);
plot(xvalues, networkperformancematrix00,'LineWidth',1,'Color',[0.7 0.7 0.7])
hold on
plot(xvalues, mean(networkperformancematrix00'),'LineWidth',3,'Color','black')
hold off
title({'Infrastructure performance vs. Event magnitude';'(no adaptation measures)'},'FontSize',16)
xlabel('Event magnitude (GW generation capacity disabled)','FontSize',15)
ylabel('Performance (fraction demand served)','FontSize',15)
xlim([0 max(xvalues)])
ylim([min(min(networkperformancematrix00)) - 0.01 1.01])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',14)
set(gca,'XTickMode','manual','YTickMode','manual')
saveas(gcf,strcat(outputpath,'heatwaveperformance_noadaptations.fig'));
saveas(gcf,strcat(outputpath,'heatwaveperformance_noadaptations.png'));

% plot the combined results
subplot(1,2,2);
plot(xvalues, meanperformancematrix(:,1),...
    xvalues, meanperformancematrix(:,2),...
    xvalues, meanperformancematrix(:,3),...
    xvalues, meanperformancematrix(:,4),...
    xvalues, meanperformancematrix(:,5),...
    xvalues, meanperformancematrix(:,6),'LineWidth',3);
title({'Infrastructure performance vs. Event magnitude';'(comparison of measures)'},'FontSize',16)
xlabel('Event magnitude (GW generation capacity disabled)','FontSize',15)
ylabel('Performance (fraction demand served)','FontSize',15)
leg = legend('0% reduction','5% reduction','10% reduction','15% reduction','20% reduction','25% reduction','Location','SouthWest');
set(leg,'FontSize',12)
v = get(leg,'title');
set(v,'string','Adaptation measure','FontSize',14);
xlim([0 max(xvalues)])
ylim([mean(min(networkperformancematrix00)) - 0.001 1.0003])
colormap('Lines')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',14)
set(gca,'XTickMode','manual','YTickMode','manual')
saveas(gcf,strcat(outputpath,'heatwaveperformance_adaptationscomparison.fig'));
saveas(gcf,strcat(outputpath,'heatwaveperformance_adaptationscomparison.png'));

set(z,'PaperUnits','inches','PaperPosition',[0 0 14 6])
%print -dpng heatwaveresults.png;
print(z,'-dpng','/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/finalplots_1000iterations/heatwaveresults.png');

min(min(networkperformancematrix00))
min(mean(networkperformancematrix00))
max(mean(networkperformancematrix00))
std(mean(networkperformancematrix00))

%% PLOT THE HEAT WAVE RESULTS - SUMMARY BOX PLOT

clear

% load the results
load('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/1000iterationsb/heatwaveperformance.mat')
outputpath = ('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/finalplots_1000iterations/');

%load('C:\Users\Andrew\Documents\TUD\AndrewBollinger\projects\MatlabResilienceModel\ThesisVersion\output\100iterations\heatwaveperformance_repaired.mat')
%outputpath = ('C:\Users\Andrew\Documents\TUD\AndrewBollinger\projects\MatlabResilienceModel\ThesisVersion\output\finalplots\');

% box plot of combined results
networkperformancematrix_total = [mean(networkperformancematrix00);...
    mean(networkperformancematrix05);...
    mean(networkperformancematrix10);...
    mean(networkperformancematrix15);...
    mean(networkperformancematrix20);...
    mean(networkperformancematrix25)];
h = boxplot(networkperformancematrix_total','boxstyle','filled','symbol','k+',...
    'labels',{'0%','5%','10%','15%','20%','25%'});
set(h,'color','black')
title({'Heat wave resilience of the infrastructure';'(comparison of measures)'},'FontSize',16)
xlabel({'Adaptation measure';'(% reduction in demand)'},'FontSize',15)
ylabel('Resilience','FontSize',15)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',14)
saveas(gcf,strcat(outputpath,'heatwaveperformance_adaptationscomparison_box.fig'));
saveas(gcf,strcat(outputpath,'heatwaveperformance_adaptationscomparison_box.png'));

%% PLOT THE FLOOD RESULTS -- FREQUENCY VS MAGNITUDE (CURRENT SITUATION)

clear

% load the results
load('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/100iterations_flooddefensescenarios/floodperformance.mat')
outputpath = ('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/finalplots/');

% line plot
xvalues = 1:numberofbusestoremove;
xvalues = xvalues - 1;
semilogy(xvalues,mean(eventprobabilitymatrix00a'),xvalues,mean(eventprobabilitymatrix00b'),xvalues,mean(eventprobabilitymatrix00c'),xvalues,mean(eventprobabilitymatrix00d'))
xlabel('Event frequency (# substations failed)')
ylabel('Event probability')

%% PLOT THE FLOOD RESULTS -- FRQUENCY-ADJUSTED PERFORMANCE VS MAGNITUDE (CURRENT SITUATION)

clear

% load the results
load('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/100iterations_flooddefensescenarios/floodperformance.mat')
outputpath = ('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/finalplots/');

xvalues = 1:numberofbusestoremove;
xvalues = xvalues - 1;
%xvalues(1) = [];
meanfrequencyadjustedperformancematrix_cropped = meanfrequencyadjustedperformancematrix;
meanfrequencyadjustedperformancematrix_cropped(1,:) = [];

plot(xvalues, meanfrequencyadjustedperformancematrix(:,1))

xlabel('Event magnitude (# substations failed)')
ylabel('Frequency-adjusted performance')


%% PLOT THE FLOOD RESULTS -- FRQUENCY VS MAGNITUDE (FUTURE SCENARIOS)



%% PLOT THE FLOOD RESULTS -- FREQUENCY-ADJUSTED PERFORMANCE VS. MAGNITUDE (FUTURE SCENARIOS)



%%