% Parameters
%epsilon = 0.00001;
epsilon = 0.001;

numChoices = 100;
participants = [5:5:50];
votesPerPerson = [5:5:50];
%participants = [40:45];
%votesPerPerson = [25:30];
deltaExitCriteria = 10;

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
        choiceVotes = zeros(1,numChoices);
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
            if deltaCounts >= deltaExitCriteria
              fprintf('dc=%f, dd=%f', delta_cons, delta_dilu);
            endif
        else
            deltaCounts = 0;
        end

        cons_avg = cons_avg_new;
        dilu_avg = dilu_avg_new;

        % increment the loop counter
        n = n + 1;

        % print progress dot
        if mod(n,1/(1000*epsilon)) == 0
          fprintf('.');
        end
    end
end

function [c] = consensus(choiceVotes, participants)
    frac = 0.5;
    c = sum(choiceVotes >= frac * participants) / length(choiceVotes);
end

function [d] = dilution(choiceVotes, participants)
    frac = 1.0;
    d = sum(choiceVotes >= frac * participants) / length(choiceVotes);
end

% Allocate Memory
c = zeros(length(participants), length(votesPerPerson));
d = zeros(length(participants), length(votesPerPerson));
n = zeros(length(participants), length(votesPerPerson));

% Run the simulations
for i = 1:length(participants)
  for j = 1:length(votesPerPerson)
    fprintf('Simulating %d participants %d votes...', participants(i), votesPerPerson(j));
    [c(i,j), d(i,j), n(i,j)] = simulate_voting( numChoices, participants(i), votesPerPerson(j), epsilon, deltaExitCriteria);
    fprintf('[c = %f, d = %f, n = %d]. Done.\n', c(i,j), d(i,j), n(i,j));
  end
end

% Plot Results
figure;
[x,y] = meshgrid(participants, votesPerPerson);
mesh(x, y, c');
title("Consensus", "FontSize", 20);
xlabel("# of Participants", "FontSize", 14);
ylabel("# of Votes per Participant", "FontSize", 14);
zlabel("Fraction", "FontSize", 14);

%plot(votesPerPerson, c, votesPerPerson, d);
%fprintf("Consensus, mean: %f\n", c);
%fprintf("Dilution, mean: %f\n", d);
%fprintf("Number of simulations: %d\n", n);


