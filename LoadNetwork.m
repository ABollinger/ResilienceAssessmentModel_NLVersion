% load the network file
networklist = fileread('current_network');

% parse the file
networklist = networklist(2:end-1);
networklist = regexp(networklist,'[[','split');

% create some variables
tick = networklist(1);
buses = networklist(2);
generators = networklist(3);
links = networklist(4);
capacities = networklist(5);

% clean up the variable values
tick = regexprep(tick,' ','');
buses = regexprep(buses,']','');
generators = regexprep(generators,']','');
links = regexprep(links,']','');
capacities = regexprep(capacities,']','');

% parse the variable values
buses = regexp(buses,'[','split');
generators = regexp(generators,'[','split');
links = regexp(links,'[','split');
capacities = regexp(capacities,'[','split');

% create the generator list
genlist = zeros(size(generators{1},2), size(sscanf(generators{1}{1},'%f'),1)); 
for x = 1:size(generators{1},2)
    genlist(x,:) = sscanf(generators{1}{x},'%f')';
end

% create the bus list
buslist = zeros(size(buses{1},2), size(sscanf(buses{1}{1},'%f'),1)); 
for x = 1:size(buses{1},2)
    buslist(x,:) = sscanf(buses{1}{x},'%f')';
end

% create the branch list
branchlist = zeros(size(links{1},2), size(sscanf(links{1}{1},'%f'),1)); 
for x = 1:size(links{1},2)
    branchlist(x,:) = sscanf(links{1}{x},'%f')';
end

% create the branch capacity list
capacitylist = zeros(size(capacities{1},2), size(sscanf(capacities{1}{1},'%f'),1)); 
for x = 1:size(capacities{1},2)
    capacitylist(x,:) = sscanf(capacities{1}{x},'%f')';
end

% create the structs
mpc.bus = buslist;
mpc.gen = genlist;
mpc.branch = branchlist;
mpc.baseMVA = 100;
mpc.version = '2';