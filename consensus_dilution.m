% Parameters
epsilon = 0.01;
numChoices = 100;
participants = 5;
votesPerPerson = 10;
deltaExitCriteria = 5;

[c, d, n] = simulate_voting( numChoices, participants, votesPerPerson, epsilon, deltaExitCriteria);

fprintf("Consensus, mean: %f\n", c);
fprintf("Dilution, mean: %f\n", d);
fprintf("Number of simulations: %d\n", n);


function [cons_avg, dilu_avg, n] = simulate_voting(numChoices, participants, votesPerPerson, epsilon, deltaExitCriteria)
    n = 0;
    cons_avg = 0;
    dilu_avg = 0;
    deltaCounts = 0;
    prefs = zeros(participants, numChoices);
    choiceVotes = zeros(1,numChoices);
    
    while deltaCounts < deltaExitCriteria
        % generate participant preferences
        for i = 1:participants
            prefs(i,:) = randperm(numChoices);
        end
    
        % Add up their votes
        for i = 1:participants
            for j = 1:votesPerPerson
                choiceVotes(prefs(i,j)) = choiceVotes(prefs(i,j)) + 1;
            end
        end
    
        % Compute consensus and dilution
        cons = consensus(choiceVotes, participants);
        dilu = dilution(choiceVotes, participants);
    
        % Calculate weighted averages
        cons_avg_new = (cons_avg * n + cons) / (n + 1);
        dilu_avg_new = (dilu_avg * n + dilu) / (n + 1);
    
        delta_cons = abs(cons_avg - cons_avg_new);
        delta_dilu = abs(dilu_avg - dilu_avg_new);
    
        delta = max( delta_cons, delta_dilu );
    
        if delta < epsilon
            deltaCounts = deltaCounts + 1;
        else
            deltaCounts = 0;
        end
    
        cons_avg = cons_avg_new;
        dilu_avg = dilu_avg_new;
    
        % increment the loop counter
        n = n + 1;

        % escape if more than max iterations
        if n > 100
            break
        end
    end
end

function [c] = consensus(choiceVotes, participants)
    frac = 0.5;
    c = 0;
    for i = 1:length(choiceVotes)
        if choiceVotes(i) >= (frac * participants)
            c = c + 1;
        end
    end

    c = c / length(choiceVotes);
end

function [d] = dilution(choiceVotes, participants)
    frac = 1.0;
    d = 0;
    for i = 1:length(choiceVotes)
        if choiceVotes(i) >= (frac * participants)
            d = d + 1;
        end
    end

    d = d / length(choiceVotes);
end


