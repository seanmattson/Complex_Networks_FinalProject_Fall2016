function [X, Z, T] = GillespieDimerization(X0, tmax, S, M, C)

DA_Identity = eye(M); %Identity Matrix for DA

%Initializing Values
t = 0;
Z = zeros(M,1);

while t<tmax
    a1 = C(1)*(Z(3) - Z(4));
    a2 = C(2)*(2 + Z(1) - Z(2) - 2*Z(9) + 2*Z(10));
    a3 = C(3)*(Z(5) - Z(6) - Z(7) + Z(8));
    a4 = C(4)*(Z(3) - Z(4));
    a5 = C(5)*(4 - Z(5) + Z(6) - Z(7) + Z(8) + Z(9) - Z(10))*(2 - Z(5) + Z(6));
    a6 = C(6)*(Z(5) - Z(6) - Z(7) + Z(8));
    a7 = C(7)*(4 - Z(5) + Z(6) - Z(7) + Z(8) + Z(9) - Z(10))*(Z(5) - Z(6) - Z(7) + Z(8));
    a8 = C(8)*(Z(7) - Z(8));
    a9 = (C(9)*(2 + Z(1) - Z(2) - 2*Z(9) + 2*Z(10))*(1 + Z(1) - Z(2) - 2*Z(9) + 2*Z(10)))/2;
    a10 = C(10)*(4 - Z(5) + Z(6) - Z(7) + Z(8) + Z(9) - Z(10));
    
    A = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10]; %Vector of Propensities
    A_sum = sum(A); %Sum of different Propensities
    
    mean_exp = 1/A_sum;
    tau = exprnd(mean_exp); %Draws value from the exponential distribution
    
    pdf = A./A_sum; % Creates discrete pdf from the propensities
    
    m_rand = rand; % Choose value randomly (0,1)
    m = find(cumsum(pdf)>m_rand,1); % Randomly chooses index value from pdf, based on probability to occur
    
    Z = Z + DA_Identity(:,m); % Updates DA
    t = t + tau; %Updates Time
end

X = X0 + S*Z;
Z = Z;
T = t;
end
