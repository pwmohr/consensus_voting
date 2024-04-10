
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear global;

global SimParams = struct(
    "voting",
    struct(
        % numChoices is the number of individual items that are available to vote upon
        "numChoices", 60,

        % participants is the number of participants being simulated
        "participants", [40:1:50],

        % votesPP is how many votes each participant can cast
        "votesPP", [3:1:15],

        % method describes how the participants distribute their votes
        % available methods:
        %
        % 'linear':
        %       All participants agree that 1 is best and that 'numChoices' is worst.
        %
        % 'random':
        %       Participants randomly allocate votes. Only 1 vote allowed per item - no vote stacking
        %
        % 'random-groupings-single-faction':
        %       The items to be voted upon are grouped (according to 'choice_groups.boundaries', below).
        %       All agree that the first group is best, next is 2nd best, etc. But within each group
        %       participants vote randomly.
        %
        % 'random-groupings-multi-faction':
        %       Like 'random-groupings-single-faction', but users are divided into factions (according
        %       to 'factions.boundaries', below), which in general do not agree on the ordering of
        %       group preferences from favorite to least favorite.
        %
        % 'random-MV'
        %       Like 'random', but allows users to stack votes.
        "method", 'random-groupings-multi-faction'
    ),
    "csr_criteria",
    struct(
        % Fractions for determining consensus, saturation, and rejection values
        "consensus_fraction", 0.8,
        "saturation_fraction", 1.0,
        "rejection_fraction", 0.1
    ),
    "choice_groups",
    struct(
        % boundaries specifies the distribution of how items being voted on are divided into
        % different preference groups
        % e.g., [0.25 0.5 0.75 1.0] defines the dividing lines to be every 25% of the total number
        % of items, so if there were 100 items, group 1 would contain 1-25, group 2 would have 26-50, etc.
        % the number of boundaries and their values can be changed, so for example if you wanted 3 groups
        % where the first one is 10% of the total, the next one is 75% and the last one is 15%, you could
        % use [0.1 0.85 1.0]
        "boundaries", [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]
    ),
    "factions",
    struct(
        % boundaries specifies the distribution of people into the various factions. The dividing lines
        % operate similarly to group boundaries, above. As with group boundaries, factions do not need to be
        % the same size as one another, and they can be added and taken away by the manner in which 'boundaries'
        % is edited.
        "boundaries", [0.25 0.5 0.75 1.0],

        % preferences determines the order of preference that each faction has for each choice group.
        % the rows represent factions 1-N, and the columns represent the order of preference for that faction,
        % with the 1st column being the favorite, and the last column the least favogrite.
        %
        % NOTE that it is important that the preferences array match the dimensions implied by the number of choice
        % groups (columns), and the number of factions (rows).
        "preferences", [1 2 3 4 5 6 7 8 9 10;
                        4 1 2 3 6 5 8 7 9 10;
                        3 4 1 2 5 6 7 8 9 10;
                        2 3 4 1 6 5 8 7 9 10],

        % vote_distribution represents what fraction of total available votes a participant will cast in each
        % of their descending list of preferences. For example, [ 0.5 0.35 0.15 ] would mean that a voter
        % would allocate 50% of their votes to their favorite choice group, 35% of their votes to their 2nd favorite
        % group, and 15% of their votes to their 3rd favorite group.
        %
        % The sum of the columns should be 1.0
        "vote_distribution", [ 0.3 0.25 0.2 0.15 0.1 ]
    ),
    "precision",
    struct(
        % epsilon and deltaExitCriteria define the conditions under which
        % the simulation ends. When the 'consensus' average value has
        % changed less than 'epsilon' for at least 'deltaExitCriteria' passes,
        % the simulation is complete
        "epsilon", 0.001,
        "deltaExitCriteria", 10
    ),
    "debug",
    struct(
        % set display_vote_monitor to true to show a bar graph of voting results during the simulation,
        % set display_vote_monitor to false to turn it off
        "display_vote_monitor", true
    )
);

%
% Allocate Memory
%
c = s = r = n = zeros(length(SimParams.voting.participants), length(SimParams.voting.votesPP));

%
% Set up bar monitor plot window
%
if SimParams.debug.display_vote_monitor == true
    ss = get(0, 'screensize');
    width = ceil(0.4 * ss(3));
    height = ceil(0.67 * width);
    x = ceil(0.55 * ss(3));
    y = ceil(0.25 * ss(4));
    figure('Position', [x y width height]);
endif

%
% Run the simulations for each # of participants and # of votes per participant
%
for i = 1:length(SimParams.voting.participants)
  for j = 1:length(SimParams.voting.votesPP)
    fprintf('Simulating %2d participants %2d votes: ', SimParams.voting.participants(i), SimParams.voting.votesPP(j));
    [c(i,j), s(i,j), r(i,j), n(i,j)] = simulate_voting(SimParams.voting.participants(i), SimParams.voting.votesPP(j));
    fprintf('\b  [c = %f, s = %e, r = %f, n = %d]\n', c(i,j), s(i,j), r(i,j), n(i,j));
  end
end

%
% Plot Results
%
width = ceil(0.8 * ss(3));
height = ceil(0.25 * width);
x = 0.5 * (ss(3) - width);
y = 0.5 * (ss(4) - height);
figure("Position", [x y width height]);
subplot(1,4,1);
axis off;
str = disp(SimParams);
underscores = strfind(str,'_');
str(underscores) = ' ';
text(0,0.5,str);

subplot(1,4,2);
[x,y] = meshgrid(SimParams.voting.participants, SimParams.voting.votesPP);
mesh(x, y, c');
title(strcat("Consensus"," ( c >= ", sprintf("%4.2f ",SimParams.csr_criteria.consensus_fraction)," )"), "FontSize", 16);
xlabel("# of Participants", "FontSize", 14);
ylabel("# of Votes per Participant", "FontSize", 14);
zlabel("Fraction", "FontSize", 14);

subplot(1,4,3);
[x,y] = meshgrid(SimParams.voting.participants, SimParams.voting.votesPP);
mesh(x, y, r');
title(strcat("Rejected"," ( r <= ", sprintf("%4.2f ",SimParams.csr_criteria.rejection_fraction)," )"), "FontSize", 16);
xlabel("# of Participants", "FontSize", 14);
ylabel("# of Votes per Participant", "FontSize", 14);
zlabel("Fraction", "FontSize", 14);

subplot(1,4,4);
[x,y] = meshgrid(SimParams.voting.participants, SimParams.voting.votesPP);
mesh(x, y, s');
title(strcat("Saturated"," ( s >= ", sprintf("%4.2f ",SimParams.csr_criteria.saturation_fraction)," )"), "FontSize", 16);
xlabel("# of Participants", "FontSize", 14);
ylabel("# of Votes per Participant", "FontSize", 14);
zlabel("Fraction", "FontSize", 14);


