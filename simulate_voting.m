%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Voting simulation function
%
% Parameters:
%         SimParams: See above.
%               p_n: number of participants in this run
%               p_v: number of votes per participant this run
%
% Returns:
%          cons_avg: Average computed consensus (0,1)
%          dilu_avg: Average computed dilution (0,1)
%                 n: Number of iterations required to achieve exit condition
%
function [cons_avg, satu_avg, reje_avg, n] = simulate_voting(p_n, p_v)
    global SimParams;
    n = 0;
    cons_avg = 0;
    satu_avg = 0;
    reje_avg = 0;
    deltaCounts = 0;
    prefs = zeros(p_n, SimParams.voting.numChoices);
    choiceVotes = zeros(1,SimParams.voting.numChoices);
    sw='\\\\\\\\\\\\\\\\||||||||||||||||////////////////----------------';

    while deltaCounts < SimParams.precision.deltaExitCriteria
        % generate participant preferences
        for i = 1:p_n
            prefs(i,:) = generate_prefs(i, p_n, p_v);
        end

        % Add up their votes
        choiceVotes = zeros(1,SimParams.voting.numChoices);
        for i = 1:p_n
            for j = 1:p_v
              choiceVotes(prefs(i,j)) = choiceVotes(prefs(i,j)) + 1;
            end
        end

        % Compute consensus, dilution, rejected
        cons = consensus(choiceVotes, p_n);
        satu = saturation(choiceVotes, p_n);
        reje = rejected(choiceVotes, p_n);

        % Calculate weighted averages
        cons_avg_new = (cons_avg * n + cons) / (n + 1);
        satu_avg_new = (satu_avg * n + satu) / (n + 1);
        reje_avg_new = (reje_avg * n + reje) / (n + 1);

        delta_cons = abs(cons_avg - cons_avg_new);
        delta_satu = abs(satu_avg - satu_avg_new);

        delta = delta_cons; % was this: max( delta_cons, delta_satu );

        if delta < SimParams.precision.epsilon
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
function [c] = consensus(choiceVotes, p_n)
    global SimParams;
    c = sum(choiceVotes >= SimParams.csr_criteria.consensus_fraction * p_n) ...
            / SimParams.voting.numChoices;
end

%
% Saturation
% What fraction of choices received the maximum number of votes?
%
function [d] = saturation(choiceVotes, p_n)
    global SimParams;
    d = sum(choiceVotes >= SimParams.csr_criteria.saturation_fraction * p_n) ...
            / SimParams.voting.numChoices;
end

%
% Rejected
% What fractionm of choices received less than 'frac' percentage of people voting for it?
%
function [s] = rejected(choiceVotes, p_n)
    global SimParams;
    s = sum(choiceVotes <= SimParams.csr_criteria.rejection_fraction * p_n) ...
            / SimParams.voting.numChoices;
end

%
% Generate Participant Prefs
%
% p_idx = participant index
function [p] = generate_prefs(p_idx, p_n, p_v)
    global SimParams;

    switch SimParams.voting.method
        case 'linear'
            p = 1:SimParams.voting.numChoices;
        case 'random'
            p = randperm(SimParams.voting.numChoices);
        case 'random-groupings-single-faction'
            p = randomGSF(SimParams.voting.numChoices);
        case 'random-groupings-multi-faction'
            p = randomGMF(p_idx, p_n);
        case 'random-MV'
            p = randomMV(SimParams.voting.numChoices);
        otherwise
            error('Unsupported option for generate_prefs().');
    endswitch
end

%
% Randomize preferences in some number of groups, all participants assumed
% to be in the same 'faction'.
%
function [p] = randomGSF(n)
    global SimParams;

    boundary_values = [0 SimParams.choice_groups.boundaries];

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
function [p] = randomGMF(p_idx, p_n, n)
    global SimParams;

    option_boundaries = [ 0 SimParams.choice_groups.boundaries];
    faction_boundaries = SimParams.factions.boundaries;
    fp = SimParams.factions.preferences;

    b = round(n * option_boundaries);
    f = round(p_n * faction_boundaries);

    % Determine the faction
    faction = sum(f < p_idx) + 1;

    % Determine preference groupings for that faction
    for i = 2:length(option_boundaries)
        p(b(i-1)+1:b(i)) = randperm(b(i)-b(i-1)) + b(fp(faction,i-1));
    endfor
end

%
% Randomize preferences but allow multiple votes on a choice for each user
%   n = number of items to vote for
function [p] = randomMV(n)
    p = ceil(rand(1,n)*n);
end



