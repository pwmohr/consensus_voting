%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = 0.001;
numChoices = 20;
participants = [4:2:20];
votesPP = [3:1:10];
deltaExitCriteria = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Voting simulation function
%
% Parameters:
%        numChoices: The number of individual items under consideration that can be voted for
%      participants: The number of people who are voting on the choices
%           votesPP: The number of votes allowed per participant
%           epsilon: The change in averages of consensus and dilution that must be achieved for loop exit
% deltaExitCriteria: The number of times in a row that change must be less than epsilon for loop exit
%
% Returns:
%          cons_avg: Average computed consensus (0,1)
%          dilu_avg: Average computed dilution (0,1)
%                 n: Number of iterations required to achieve exit condition
%
function [cons_avg, satu_avg, reje_avg, n] = simulate_voting(numChoices, participants, votesPP, epsilon, deltaExitCriteria)
    n = 0;
    cons_avg = 0;
    satu_avg = 0;
    reje_avg = 0;
    deltaCounts = 0;
    prefs = zeros(participants, numChoices);
    choiceVotes = zeros(1,numChoices);
    sw='\\\\\\\\||||||||////////--------';

    while deltaCounts < deltaExitCriteria
        % generate participant preferences
        for i = 1:participants
            prefs(i,:) = generate_prefs(i, participants, numChoices);
        end

        % Add up their votes
        choiceVotes = zeros(1,numChoices);
        for i = 1:participants
            for j = 1:votesPP
              choiceVotes(prefs(i,j)) = choiceVotes(prefs(i,j)) + 1;
            end
        end

        % Compute consensus, dilution, rejected
        cons = consensus(choiceVotes, participants);
        satu = saturation(choiceVotes, participants);
        reje = rejected(choiceVotes, participants);

        % Calculate weighted averages
        cons_avg_new = (cons_avg * n + cons) / (n + 1);
        satu_avg_new = (satu_avg * n + satu) / (n + 1);
        reje_avg_new = (reje_avg * n + reje) / (n + 1);

        delta_cons = abs(cons_avg - cons_avg_new);
        delta_satu = abs(satu_avg - satu_avg_new);

        delta = max( delta_cons, delta_satu );

        if delta < epsilon
            deltaCounts = deltaCounts + 1;
        else
            deltaCounts = 0;
        end

        cons_avg = cons_avg_new;
        satu_avg = satu_avg_new;
        reje_avg = reje_avg_new;

        % increment the loop counter
        n = n + 1;

        % print progress widget
        fprintf('\b%c', sw(mod(n,length(sw))+1));
    end
end

%
% Consensus
% What fraction of choices received more than 'frac' percentage of people voting for it?
%
function [c] = consensus(choiceVotes, participants)
    frac = 0.5;
    c = sum(choiceVotes >= frac * participants) / length(choiceVotes);
end

%
% Saturation
% What fraction of choices received the maximum number of votes?
%
function [d] = saturation(choiceVotes, participants)
    frac = 1.0;
    d = sum(choiceVotes >= frac * participants) / length(choiceVotes);
end

%
% Rejected
% What fractionm of choices received less than 'frac' percentage of people voting for it?
%
function [s] = rejected(choiceVotes, participants)
    frac = 0.1;
    s = sum(choiceVotes <= frac * participants) / length(choiceVotes);
end

%
% Generate Participant Prefs
%
% p_idx = participant index
%   p_n = number of participants
%     n = number of choices to vote for
function [p] = generate_prefs(p_idx, p_n, n)
    type = 'random-groupings-single-faction';
    switch type
        case 'random'
            p = randperm(n);
        case 'linear'
            p = 1:n;
        case 'random-groupings-single-faction'
            p = randomGSF(n);
        case 'random-groupings-multi-factions-no-overlap'
            p = randomGMFNO(p_idx, p_n, n)
        otherwise
            error('Unsupported option for generate_prefs().');
    endswitch
end

%
% Randomize preferences in some number of groups, all participants assumed
% to be in the same 'faction'.
%
function [p] = randomGSF(n)
    boundary_values = [ 0 0.25 0.5 0.75 1.0 ];

    % Determine group boundary indices
    b = round(n * [boundary_values]);
    for i = 2:length(boundary_values)
        p(b(i-1)+1:b(i)) = randperm(b(i)-b(i-1)) + b(i-1);
    endfor
end

%
% Randomize preferences in some number of groups, participants divided into
% some number of non-overlapping factions
%
function [p] = randomGMFNO(p_idx, p_n, n)
    option_boundaries = [ 0 0.25 0.5 0.75 1.0 ];
    faction_boundaries  = [ 0.25 0.5 0.75 1.0 ];
    fp                  = [ 1 2 3 4;
                            2 3 4 1;
                            3 4 1 2;
                            4 1 2 3];

    b = round(n * option_boundaries);
    f = round(p_n * faction_boundaries);

    % Determine the faction
    faction = sum(f < p_idx) + 1;

    % Determine preference groupings for that faction
    for i = 2:length(option_boundaries)
        p(b(i-1)+1:b(i)) = randperm(b(i)-b(i-1)) + b(fp(faction,i-1));
    endfor
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF SCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Allocate Memory
%
c = s = r = n = zeros(length(participants), length(votesPP));

%
% Run the simulations
%
for i = 1:length(participants)
  for j = 1:length(votesPP)
    fprintf('Simulating %2d participants %2d votes: ', participants(i), votesPP(j));
    [c(i,j), s(i,j), r(i,j), n(i,j)] = simulate_voting( numChoices, participants(i), votesPP(j), epsilon, deltaExitCriteria);
    fprintf('\b  [c = %f, s = %e, r = %f, n = %d]\n', c(i,j), s(i,j), r(i,j), n(i,j));
  end
end

%
% Plot Results
%
figure;
[x,y] = meshgrid(participants, votesPP);
mesh(x, y, c');
title("Consensus", "FontSize", 20);
xlabel("# of Participants", "FontSize", 14);
ylabel("# of Votes per Participant", "FontSize", 14);
zlabel("Fraction", "FontSize", 14);

figure;
[x,y] = meshgrid(participants, votesPP);
mesh(x, y, s');
title("Saturated", "FontSize", 20);
xlabel("# of Participants", "FontSize", 14);
ylabel("# of Votes per Participant", "FontSize", 14);
zlabel("Fraction", "FontSize", 14);

figure;
[x,y] = meshgrid(participants, votesPP);
mesh(x, y, r');
title("Rejected", "FontSize", 20);
xlabel("# of Participants", "FontSize", 14);
ylabel("# of Votes per Participant", "FontSize", 14);
zlabel("Fraction", "FontSize", 14);
