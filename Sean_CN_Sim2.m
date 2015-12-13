% Sean CN Second Graph Langevin Gaussian Leaping

X0 = [2;
      4;
      0;
      2;
      0;
      0]; %Initial State For X

Z0 = zeros(10,1); %Initial State for Z

S = [1, -1, 0, 0, 0, 0, 0, 0, -2, 2;
     0, 0, 0, 0, -1, 1, -1, 1, 1, -1;
     0, 0, 1, -1, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, -1, 1, 0, 0, 0, 0;
     0, 0, 0, 0, 1, -1, -1, 1, 0, 0;
     0, 0, 0, 0, 0, 0, 1, -1, 0, 0]; %The Stoichiometry Matrix for the System

c1 = 0.043;
c2 = 0.0007;
c3 = 0.0715;
c4 = 0.0039;
c5 = 0.02;
c6 = 0.4791;
c7 = 0.002;
c8 = 0.8765E-11;
c9 = 0.083;
c10 = 0.5; %These are the constants for the transitions between states

C = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10]; %Placed in vetor for function input

M=10;
tmax = 600;
numbertrials = 1000;

X1 = zeros(numbertrials,1);
X2 = zeros(numbertrials,1);
X3 = zeros(numbertrials,1);

for i = 1:numbertrials
    [X,Z,T] = PoissonApprox(X0, tmax, S, M, C);
    X1(i) = X(1);
    X2(i) = X(2);
    X3(i) = X(3);
end

large_x1 = max(X1); % Find largest number of occurrences that occurred 
large_x2 = max(X2);
large_x3 = max(X3);

figure()
subplot(1,3,1);
h1 = histogram(X1,'Normalization', 'probability')
title('Probability of Monomers')
xlabel('Number of Monomers')
ylabel('Probability at 10 min')

subplot(1,3,2);
h2 = histogram(X2, 'Normalization','probability')
title('Probability of Dimers')
xlabel('Number of Dimers')

subplot(1,3,3);
h3 = histogram(X3,'Normalization','probability')
title('Probability of mRNA')
xlabel('Number of mRNA')

