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
    global SimParams

    % Choose which voting method to use
    %
    % *_csr is the consensus, saturation, and rejection method using the
    %       epsilon / delta exit criteria
    %
    % *_avt is the add-votes-to-threshold method
    if( isequal( SimParams.voting.exit_criteria, 'epsilon-deltaExitCriteria' ) )
        [cons_avg, satu_avg, reje_avg, n] = simulate_voting_csr(p_n, p_v);
    else
        [cons_avg, satu_avg, reje_avg, n] = simulate_voting_avt(p_n, p_v);
    endif
end

% Simulates voting in which participants are granted more votes until enough of
% the choices are selected.
%
% To be considered 'selected', a choice must receive 'votes-required-for-selection',
% and 'enough' are selected when the fraction of choices meets or exceeds 'choices-selected-for-exit'
%
function [cons_avg, satu_avg, reje_avg, n] = simulate_voting_avt(p_n, p_v)
    global SimParams;
    n = 0;
    cons_avg = satu_avg = reje_avg = 0;

    prefs = zeros(p_n, SimParams.voting.numChoices);
    choiceVotes = zeros(1,SimParams.voting.numChoices);
    votesReqdForSelection = SimParams.voting.votes_required_for_selection * p_n;
    p_v = 0;

    while sum( choiceVotes >= votesReqdForSelection ) < ( SimParams.voting.choices_selected_to_exit * SimParams.voting.numChoices )

        % Determine number of votesPP
        % if we are at the end of the votesPP specified,
        % return n = -1 to indicate criteria not achieved
        if n+1 <= length( SimParams.voting.votesPP )
            p_v = SimParams.voting.votesPP(n+1);
        else
            n = -1;
            return;
        end

        % generate voter preferences
        for i = 1:p_n
            prefs(i,:) = generate_prefs(i, p_n, p_v);
        endfor

        % Add up their votes
        choiceVotes = zeros(1,SimParams.voting.numChoices);
        for i = 1:p_n
            for j = 1:p_v
                if prefs(i,j) != 0
                    choiceVotes(prefs(i,j)) = choiceVotes(prefs(i,j)) + 1;
                end
            endfor
        endfor

        if SimParams.debug.display_vote_monitor == true
            display_vote_monitor( choiceVotes, p_n, p_v );
        end
        n = n + 1;
    endwhile
    n = p_v;
end

%
% Simulates voting in which the exit criteria is based on epislon / delta
% where voting ends when the average levels of consensus specified in SimParams
% is no longer changing very much.
%
function [cons_avg, satu_avg, reje_avg, n] = simulate_voting_csr(p_n, p_v)
    global SimParams;
    n = 0;
    cons_avg = 0;
    satu_avg = 0;
    reje_avg = 0;
    deltaCounts = 0;
    prefs = zeros(p_n, SimParams.voting.numChoices);
    choiceVotes = zeros(1,SimParams.voting.numChoices);
    sw='\\\\\\\\\\\|||||||||||///////////-----------';

    % Rejection Labeling
    while deltaCounts < SimParams.precision.deltaExitCriteria
        % generate participant preferences
        for i = 1:p_n
            prefs(i,:) = generate_prefs(i, p_n, p_v);
        end

        % Add up their votes
        choiceVotes = zeros(1,SimParams.voting.numChoices);
        for i = 1:p_n
            for j = 1:p_v
              if prefs(i,j) != 0
                    choiceVotes(prefs(i,j)) = choiceVotes(prefs(i,j)) + 1;
              end
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

            if SimParams.debug.display_vote_monitor == true
                display_vote_monitor( choiceVotes, p_n, p_v );
            endif
        else
            deltaCounts = 0;
        end

        cons_avg = cons_avg_new;
        satu_avg = satu_avg_new;
        reje_avg = reje_avg_new;

        % increment the loop counter
        n = n + 1;

        % print progress widget
        % fprintf('\b%c', sw(mod(n,length(sw))+1));
    end
end

%
% Bar Graph Vote Monitor
%
function [] = display_vote_monitor( choiceVotes, p_n, p_v )
    global SimParams;

    bar( choiceVotes );
    textVertOffset = 0.5;
    axis([0 SimParams.voting.numChoices 0 max(SimParams.voting.participants)]);

    if isequal( SimParams.voting.exit_criteria, 'epsilon-deltaExitCriteria' )
        line([0 SimParams.voting.numChoices], [1 1]*SimParams.csr_criteria.rejection_fraction*p_n, 'color', 'red');
        text(0.95*SimParams.voting.numChoices, SimParams.csr_criteria.rejection_fraction*p_n+textVertOffset, 'Rejection', ...
             'fontsize', 12, 'color', 'black', 'horizontalAlignment', 'right');

        % Consensus Labeling
        line([0 SimParams.voting.numChoices], [1 1]*SimParams.csr_criteria.consensus_fraction*p_n, 'color', 'green');
        text(0.95*SimParams.voting.numChoices, SimParams.csr_criteria.consensus_fraction*p_n+textVertOffset, 'Consensus', ...
             'fontsize', 12, 'color', 'black', 'horizontalAlignment', 'right');

        % Saturation Labeling
        line([0 SimParams.voting.numChoices], [1 1]*SimParams.csr_criteria.saturation_fraction*p_n, 'color', 'blue');
        text(0.95*SimParams.voting.numChoices, SimParams.csr_criteria.saturation_fraction*p_n+textVertOffset, 'Saturation',
             'fontsize', 12, 'color', 'black', 'horizontalAlignment', 'right');
    end

    xlabel('Choices');
    ylabel('Total Votes');
    title(sprintf('Voting Bar Graph, %d Participants, %d Votes Each', p_n, p_v));
    drawnow;
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
            p = randomGMF(p_idx, p_n, p_v);
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
% p_idx = index for this participant
% p_n   = number of participants
% p_v   = number of votes per participant
function [p] = randomGMF(p_idx, p_n, p_v)
    global SimParams;

    option_boundaries = [ 0 SimParams.choice_groups.boundaries];
    faction_boundaries = SimParams.factions.boundaries;
    fp = SimParams.factions.preferences;

    b = round(SimParams.voting.numChoices * option_boundaries);
    f = round(p_n * faction_boundaries);

    % Determine the faction for this voter
    faction = sum(f < p_idx) + 1;

    % Determine how many votes this participant will allocate across each group
    numGroups = length(SimParams.choice_groups.boundaries);
    numGroupsToVote = length(SimParams.factions.vote_distribution);
    votesPerGroup = [round(SimParams.factions.vote_distribution * p_v) zeros(1,numGroups-numGroupsToVote)];

    % Due to rounding, sometimes the above few statements can cause a voter to cast more votes than allowed
    % or fewer votes than permitted. The code below corrects this.
    while sum(votesPerGroup) < p_v
        % If not enough votes allocated, add to favorite group
        votesPerGroup(1) = votesPerGroup(1) + 1;
    endwhile
    while sum(votesPerGroup) > p_v
        % If too many votes allocated, take away from least favorite group that has
        % at least one vote
        for i = numGroupsToVote:-1:1
            if votesPerGroup(i) > 0
                votesPerGroup(i) = votesPerGroup(i) - 1;
                break
            endif
        endfor
    endwhile

    % Depending on how many groups there are, it may be the case that more votes are allocated to a group%
    % than the number of choices available in that group. Unless stacking votes is permitted, this isn't okay,
    % so this code fixes that by going through preferences and reallocating votes in that manner.

    % First collect the excessVotes
    excessVotes = 0;
    for i = 2:length(option_boundaries)
        excessVotes = excessVotes + max(votesPerGroup(i-1) - (b(i)-b(i-1)),0);
        votesPerGroup(i-1) = votesPerGroup(i-1) - max(votesPerGroup(i-1) - (b(i)-b(i-1)),0);
    endfor

    % Then redistribute them in order of preferences
    for i = 2:length(option_boundaries)
        % Add votes, in order of group preference until they're full or we run out of excessVotes
        while votesPerGroup(fp(faction,i-1)) < (b(i)-b(i-1)) && excessVotes > 0
            votesPerGroup(fp(faction,i-1)) = votesPerGroup(fp(faction,i-1)) + 1;
            excessVotes = excessVotes - 1;
        endwhile
    endfor

    % Determine preference groupings for that faction
    p = zeros(1,SimParams.voting.numChoices);
    votesTaken = 0;
    for i = 2:length(option_boundaries)
        p(votesTaken+1:votesTaken+votesPerGroup(i-1)) = randperm(b(i)-b(i-1),votesPerGroup(i-1)) + b(fp(faction,i-1));
        votesTaken = votesTaken + votesPerGroup(i-1);
    endfor

end

%
% Randomize preferences but allow multiple votes on a choice for each user
%   n = number of items to vote for
function [p] = randomMV(n)
    p = ceil(rand(1,n)*n);
end



