
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear global;

global SimParams = struct(
    "voting",
    struct(
        % numChoices is the number of individual items that are available to vote upon
        "numChoices", 100,

        % participants is the number of participants being simulated. This can be an
        % array to run the simulation multiple times.
        % e.g., [30 30 30 30], or equivalently, [30]*ones(1,5), would run the simulation 5 times
        % each time with 30 participants.
        % e.g., [10:5:25] would run the simulation with 10, then 15, then 20, then 25 participants.
        "participants", [30]*ones(1,5),

        % votesPP is how many votes each participant can cast
        "votesPP", [1:100],

        % method describes how the participants distribute their votes
        % available methods:
        %
        % 'linear':
        %       All participants agree that 1 is best and that 'numChoices' is worst.
        %
        % 'random':
        %       Participants randomly allocate votes. Only 1 vote allowed per item - no vote stacking
        %
        % 'random-groupings-single-faction': [DEPRECATED]
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
        %
        "method", 'random-groupings-multi-faction',

        % exit-criteria describes what causes the voting routine to exit and return the results
        %
        % 'epsilon-deltaExitCriteria':
        %       Looks at how much the average consensus number changes with each
        %       iteration and once it is below epsilon for deltaExitCriteria cycles, it exits.
        %
        % 'add-votes-until-threshold':
        %       Ignores epsilon / deltaExitCriteria, and instead exits once consensus is first achieved,
        %       and reports number of votes required to achieve first consensus.
        "exit_criteria", 'add-votes-until-threshold',

        % Min votes required to consider one of the choices 'selected'. It should be a number from 0 to 1 and is a fraction of the number
        % of voting participants. This is used in the 'add-votes-until-threshold' exit criteria and is ignored in the 'epsilon-deltaExitCriteria'
        % exit criteria.
        "votes_required_for_selection", .75,

        % Min choices selected to exit. This is the number of choices that must be selected (according to 'votes-required-for-selection')
        % in order to exit the voting. It is a fraction of the number of choices, so should be a number from 0 to 1.
        "choices_selected_to_exit", .2
   ),
    "csr_criteria",
    struct(
        % Fractions for determining consensus, saturation, and rejection values - used in the 'epsilon-deltaExitCriteria' exit criteria,
        % ignored in 'add-votes-until-threshold'
        "consensus_fraction", 10/60,
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
        "boundaries", [.25 .5 .75 1.0],

        % preferences determines the order of preference that each faction has for each choice group.
        % the rows represent factions 1-N, and the columns represent the order of preference for that faction,
        % with the 1st column being the favorite, and the last column the least favorite.
        %
        % NOTE that it is important that the preferences array match the dimensions implied by the number of choice
        % groups (columns), and the number of factions (rows).
        "preferences", [1 2 3 4  5 6 7 8 9 10;
                        4 1 2 3  6 5 8 7 9 10;
                        7 8 9 10 3 4 1 2 5 6;
                        8 7 9 10 2 3 4 1 6 5],

        % vote_distribution represents what fraction of each group of choices the participant wants to vote for
        % e.g., [ .75 .5 .25 .1 ] would mean the participant would want to place votes on 75% of their favorite
        % choice group, but then would move on to start voting on their second favorite choice group. Once 40%
        % of those choices had been voted for, the participant would move on to the 3rd group and so on until they
        % are out of votes. These do not need to add up to anything (such as 1). They simply reflect how strongly
        % a participant wants to max out their votes in each category.
        %
        % If after they've hit those percentage targets, they have votes remaining, they will fill up the max number
        % of votes in their first preference choice group, then fill up the max number in their 2nd preferred choice group, etc.
        "vote_distribution", linspace(.8,0,10)
    ),
    "precision",
    struct(
        % epsilon and deltaExitCriteria define the conditions under which
        % the simulation ends for the 'epsilon-deltaExitCriteria' exit_critera.
        % When the 'consensus' average value has changed less than 'epsilon'
        % for at least 'deltaExitCriteria' passes, the simulation is complete
        "epsilon", 0.001,
        "deltaExitCriteria", 10
    ),
    "debug",
    struct(
        % set display_vote_monitor to true to show a bar graph of voting results during the simulation,
        % set display_vote_monitor to false to turn it off
        %
        % The display will be updated according to "display_update_frequency", 1 = every time, n = every nth pass.
        "display_vote_monitor", true,
        "display_update_frequency", 1
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
% Consensus-Saturation method
% Run the simulations for each # of participants and # of votes per participant
%
if isequal( SimParams.voting.exit_criteria, 'epsilon-deltaExitCriteria' )
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
else
    for i = 1:length(SimParams.voting.participants)
        [c,s,r,n(i)] = simulate_voting(SimParams.voting.participants(i), max(SimParams.voting.votesPP));
        if n(i) > 0
            fprintf("For %d participants, exit condition achieved at %d votes per participant.\n", SimParams.voting.participants(i), n(i));
        else
            fprintf("For %d participants, no exit with up to %d votes.\n", SimParams.voting.participants(i), max(SimParams.voting.votesPP));
        end
    end

    m = mean(n(:,1));
    sd = std(n(:,1));

    if m > 0 && length(SimParams.voting.participants) > 1
        fprintf("Mean number of votes to match exit criteria: %f\n", m);
        fprintf("Standard deviation: %f\n", sd);
    endif
endif

