
if printstuff == true
    disp('Number of links remaining = ');
    disp(numrows);
end

% construct the adjacency matrix from the bus and branch data
adjacencymatrix = zeros(length(mpc.bus(:,1)),length(mpc.bus(:,1)));
for g = 1:size(mpc.branch(:,1))
    adjacencymatrix(mpc.branch(g,1), mpc.branch(g,2)) = 1;
    adjacencymatrix(mpc.branch(g,2), mpc.branch(g,1)) = 1;
end

% add ones to the diagonal of the adjacency matrix
for h = 1:length(adjacencymatrix)
    adjacencymatrix(h,h) = 1;
end

% identify the nodes belonging to the different components
m = adjacencymatrix;
[a,b,c,d] = dmperm(m);
components = zeros(length(c) - 1, numrows);
for y = 2:length(c)
    nodevector = a(c(y - 1):c(y) - 1);
    for z = 1:length(nodevector)
        components(y-1, z) = nodevector(z);
    end
end

if printstuff == true
    disp('components = ');
    disp(components);
end

% create empty matrices to contain the bus results and branch flows for the whole network
busresults2 = horzcat(mpc.bus(:,1), zeros(length(mpc.bus(:,1))));
branchflows2 = horzcat(mpc.branch(:,1:2), zeros(length(mpc.branch(:,1)),2));

% for each component
for component = 1:length(components(:,1))

    % create a vector of the nodes in the component
    nodesinthiscomponent = nonzeros(components(component,:));

    mpc2.bus = mpc.bus(nodesinthiscomponent,:);
    mpc2.gen = mpc.gen(ismember(mpc.gen(:,1), nodesinthiscomponent),:);
    mpc2.branch = mpc.branch(ismember(mpc.branch(:,1), nodesinthiscomponent),:);
    mpc2.branch = mpc2.branch(ismember(mpc2.branch(:,2), nodesinthiscomponent),:);

    if printstuff == true
        disp('bus list = ');
        disp(mpc2.bus);
        disp('gen list = ');
        disp(mpc2.gen);
        disp('branch list = ');
        disp(mpc2.branch);
    end
    
    % redispatch generators if there is not sufficient generation
    if sum(mpc2.bus(:,3)) > sum(mpc2.gen(:,2))
        generationdeficit = sum(mpc2.bus(:,3)) - sum(mpc2.gen(:,2));
        remaininggenerationavailable = sum(mpc2.gen(:,9)) - sum(mpc2.gen(:,2));
        if remaininggenerationavailable > generationdeficit
            availablegenerationvector = mpc2.gen(:,9) - mpc2.gen(:,2);
            generationfactor = generationdeficit / remaininggenerationavailable;
            supplementarygenerationvector = availablegenerationvector * generationfactor;
            mpc2.gen(:,2) = mpc2.gen(:,2) + supplementarygenerationvector;
        else 
            mpc2.gen(:,2) = mpc2.gen(:,9);
        end
    end 

    % reset the demand of the buses in case there is not enough
    % generation capacity
    totalgenerationcapacityinthiscomponent = sum(mpc2.gen(:,2));
    totalloadinthiscomponent = sum(mpc2.bus(:,3));
    if totalgenerationcapacityinthiscomponent < totalloadinthiscomponent
        mpc2.bus(:,3) = mpc2.bus(:,3) * totalgenerationcapacityinthiscomponent / totalloadinthiscomponent;
    end
    
%     if sum(mpc2.gen(:,9)) < sum(mpc2.bus(:,3))
%         disp('Demand exceeds supply')
%     end

    % initially set the bus types of all buses to 1
    mpc2.bus(:,2) = 1;

    % set the bus types of all buses with generators attached to 2
    buseswithgeneratorsattached = [];
    numbersofbuseswithattachedgenerators = mpc2.gen(:,1);
    for j = 1:length(numbersofbuseswithattachedgenerators)
        buseswithgeneratorsattached = vertcat(buseswithgeneratorsattached, find(mpc2.bus(:,1) == numbersofbuseswithattachedgenerators(j)));
    end
    mpc2.bus(buseswithgeneratorsattached,2) = 2;

    %identify the isolated buses and set their bus types = 4
    %isolated buses are buses with no attached demand, no attached 
    %generators and only one attached line
    buseswithnodemand = find(mpc2.bus(:,3) == 0); 
    buseswithnogenerators = find(mpc2.bus(:,2) == 1);
    isolatedbuses = intersect(buseswithnodemand, buseswithnogenerators);

    busesofbranches = vertcat(mpc2.branch(:,1), mpc2.branch(:,2));
    [busfrequencies, busnos] = hist(busesofbranches, unique(busesofbranches));
    buseswithonebranch = [];
    numbersofbuseswithonebranch = busnos(find(busfrequencies == 1));
    for k = 1:length(numbersofbuseswithonebranch)
        buseswithonebranch = vertcat(buseswithonebranch, find(mpc2.bus(:,1) == numbersofbuseswithonebranch(k)));
    end

    isolatedbuses = intersect(isolatedbuses, buseswithonebranch);
    mpc2.bus(isolatedbuses,2) = 4;

    %create a matrix to contain the bus results for this component
    busresults = horzcat(mpc2.bus(:,1), zeros(length(mpc2.bus(:,1)),1));

    % if there are generators and branches in this component
    if length(find(mpc2.bus(:,2) == 2)) > 0 && length(mpc2.branch(:,1)) > 0 && length(find(mpc2.bus(:,2) ~= 4)) > 1

        % if there is no slack bus, create one
        if length(find(mpc2.bus(:,2) == 3)) == 0  
            buseswithgenerators = find(mpc2.bus(:,2) == 2);
            mpc2.bus(buseswithgenerators(1),2) = 3;
        end

        % fill in the missing matrices for the power flow analysis
        mpc2.version = mpc.version;
        mpc2.baseMVA = mpc.baseMVA;
        
        % remove the unnecessary columns from the bus and branch matrices
        mpc2.bus = mpc2.bus(:,1:13);
        mpc2.branch = mpc2.branch(:,1:13);
        mpc2.gen = mpc2.gen(:,1:21);
        
        % run the power flow analysis and save the results
        results = runpf(mpc2, mpopt);
        
        % get the branch injections
        branchinjectedfromend = results.branch(:,14);
        branchinjectedtoend = results.branch(:,16);
        branchresults = horzcat(mpc2.branch(:,1:2), branchinjectedfromend, branchinjectedtoend);

        % determine the injections to/from each bus with the connected branches
        for y = 1:length(branchresults(:,1))
            for z = 1:length(busresults(:,1))

                if busresults(z,1) == branchresults(y,1) 
                    busresults(z,2) = busresults(z,2) - branchresults(y,3);
                end

                if busresults(z,1) == branchresults(y,2)
                    busresults(z,2) = busresults(z,2) - branchresults(y,4);
                end

            end
        end

        % get the generator results
        genresults = results.gen(:,1:2);

        % determine the injections to each bus from the connected generators
        for v = 1:length(genresults(:,1))
            for w = 1:length(busresults(:,1))

                if busresults(w,1) == genresults(v,1)
                    busresults(w,2) = busresults(w,2) + genresults(v,2);
                end

            end
        end
        
        % add the bus results for this component to the bus results matrix
        % for the whole network.
        for v = 1:length(busresults(:,1))
            for w = 1:length(busresults2(:,1))
                if busresults(v,1) == busresults2(w,1)
                    busresults2(w,2) = busresults(v,2);
                end
            end
        end
            
        if printstuff == true
            disp('bus results = ');
            disp(busresults2);
        end
        
        % determine the flows through each branch
        branchflows = max(abs(results.branch(:,14)), abs(results.branch(:,16)));
        branchflows = horzcat(mpc2.branch(:,1:2), branchflows, zeros(length(branchflows(:,1)),1));
        
        % add the branch results for this component to the branch results
        % matrix for the whole network
        for v = 1:length(branchflows(:,1))
            for w = 1:length(branchflows2(:,1))
                if branchflows(v,1) == branchflows2(w,1) && branchflows(v,2) == branchflows2(w,2) 
                    branchflows2(w,3) = branchflows(v,3);
                end
            end
        end
        
    end
end

% determine whether the flows through each branch is exceeding its capacity
% we use a value of 1.2 * nominal capacity as maximum capacity.
% we exclude branches with capacity <= 0 (because these are transformers, 
% which we don't care about here).
for s = 1:length(branchflows2(:,1))
    for t = 1:length(capacitylist(:,1))
        if branchflows2(s,1) == capacitylist(t,1) && branchflows2(s,2) == capacitylist(t,2) 
            branchflows2(s,4) = capacitylist(t,3);
        end
    end
end
branchesexceedingcapacity = branchflows2(branchflows2(:,3) > 1.2 * branchflows2(:,4),:);
branchesexceedingcapacity = branchesexceedingcapacity(branchesexceedingcapacity(:,4) > 0,:); %ignore the transformers
branchesexceedingcapacity = branchesexceedingcapacity(:,1:2);

% if one or more branches are exceeding their capacity, randomly redispatch generators until this is no longer the case
% or until the redispatch attempts limit has been reached.
if allowgeneratorredispatching == true && ~isempty(branchesexceedingcapacity) && redispatchattempts <= redispatchattemptslimit

    % create a random generation dispatch distribution
    generatoroutputs = random('uniform', 0, 1, length(mpc.gen(:,1)), 1);
    generatoroutputs = generatoroutputs .* mpc.gen(:,9);
    generatoroutputs = generatoroutputs * sum(mpc.bus(:,3)) / sum(generatoroutputs);
    mpc.gen(:,2) = generatoroutputs;

    redispatchattempts = redispatchattempts + 1;

    CalculatePowerFlows;

% delete the branches that are exceeding their capacity and calculate the
% power flows again.  If there are no lines left, then set all the bus
% injections to zero.
elseif allowcascadingfailures == true && ~isempty(branchesexceedingcapacity)
    
    %disp('Lines exceeding capacity')

    indicesofbranchestoremove = [];
    for q = 1:length(mpc.branch(:,1))
        for r = 1:length(branchesexceedingcapacity(:,1))
            if mpc.branch(q,1) == branchesexceedingcapacity(r,1) && mpc.branch(q,2) == branchesexceedingcapacity(r,2)
                indicesofbranchestoremove = vertcat(indicesofbranchestoremove, q);
            end
        end
    end
    mpc.branch(indicesofbranchestoremove,:) = [];

    if length(mpc.branch(:,1)) > 0
        CalculatePowerFlows;
    else
        busresults2(:,2) = zeros(length(busresults2(:,1)),1);
    end
end

