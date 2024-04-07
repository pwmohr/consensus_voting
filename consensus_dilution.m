%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = 0.01;
numChoices = 100;
participants = [5:5:50];
votesPP = [5:5:50];
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
function [cons_avg, dilu_avg, n] = simulate_voting(numChoices, participants, votesPP, epsilon, deltaExitCriteria)
    n = 0;
    cons_avg = 0;
    dilu_avg = 0;
    deltaCounts = 0;
    prefs = zeros(participants, numChoices);
    choiceVotes = zeros(1,numChoices);
    sw='\\\\\\\\||||||||////////--------';

while deltaCounts < deltaExitCriteria
        % generate participant preferences
        for i = 1:participants
            prefs(i,:) = randperm(numChoices);
        end

        % Add up their votes
        choiceVotes = zeros(1,numChoices);
        for i = 1:participants
            for j = 1:votesPP
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

        % print progress widget
        fprintf('\b%c', sw(mod(n,length(sw))+1));
    end
end

%
% Consensus
%
function [c] = consensus(choiceVotes, participants)
    frac = 0.5;
    c = sum(choiceVotes >= frac * participants) / length(choiceVotes);
end

%
% Dilution
%
function [d] = dilution(choiceVotes, participants)
    frac = 1.0;
    d = sum(choiceVotes >= frac * participants) / length(choiceVotes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF SCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Allocate Memory
%
c = zeros(length(participants), length(votesPP));
d = zeros(length(participants), length(votesPP));
n = zeros(length(participants), length(votesPP));

%
% Run the simulations
%
for i = 1:length(participants)
  for j = 1:length(votesPP)
    fprintf('Simulating %2d participants %2d votes: ', participants(i), votesPP(j));
    [c(i,j), d(i,j), n(i,j)] = simulate_voting( numChoices, participants(i), votesPP(j), epsilon, deltaExitCriteria);
    fprintf('\b  [c = %f, d = %e, n = %d]\n', c(i,j), d(i,j), n(i,j));
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
mesh(x, y, d');
title("Dilution", "FontSize", 20);
xlabel("# of Participants", "FontSize", 14);
ylabel("# of Votes per Participant", "FontSize", 14);
zlabel("Fraction", "FontSize", 14);
