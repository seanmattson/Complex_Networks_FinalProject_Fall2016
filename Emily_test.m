function Emily_test
c1 = 0.043; % [s^-1] RNA -> RNA + M
c2 = 0.0007; % [s^-1] M -> 0
c3 = 0.0715; % [s^-1] DNA*D -> RNA + DNA*D
c4 = 0.0039; % [s^-1] RNA -> 0
c5 = 0.02; % [s^-1] DNA + D -> DNA*D
c6 = 0.4791; % [s^-1] DNA*D -> DNA + D
c7 = 0.002; % [s^-1] DNA*D + D -> DNA*2D
c8 = 0.8765*(10^-11); % [s^-1] DNA*2D -> DNA*D + D
c9 = 0.083; % [s^-1] M + M -> D
c10 = 0.5; % [s^-1] D -> M + M
    
%     Z = [Z1; Z2; Z3; Z4; Z5; Z6; Z7; Z8; Z9; Z10];
%     
%     X1 = 2 + Z1 - Z2 - 2*Z9 + 2*Z10; % M (Protein monomer)
%     X2 = 4 - Z5 + Z6 - Z7 + Z8 + Z9 - Z10; % D (Transcription factor dimer)
%     X3 = Z3 - Z4; % mRNA
%     X4 = 2 - Z5 + Z6; % DNA template free of dimers
%     X5 =  Z5 - Z6 - Z7 + Z8; % DNA template bound at R1
%     X6 = Z7 - Z8; % DNA template bound at R1 and R2
    
X = [2; %M
     4; %D
     0; %RNA
     2; %DNA
     0; %DNA*D
     0];%DNA*2D

%State Matrix
D = [1 -1  0  0  0  0  0  0 -2  2;
     0  0  0  0 -1  1 -1  1  1 -1;
     0  0  1 -1  0  0  0  0  0  0;
     0  0  0  0 -1  1  0  0  0  0;
     0  0  0  0  1 -1 -1  1  0  0;
     0  0  0  0  0  0  1 -1  0  0];
 
tmax = 600; %max time
Nalloc = 10000; %memory allocation
Nchs = 6; %number of species

iter = 1;
T = zeros(1,Nalloc); %empty time array
N = zeros(Nchs,Nalloc); %empty species array
t=0;

tic
while t <= tmax
    
    a1 = c1*X(3);
    a2 = c2*X(1);
    a3 = c3*X(5);
    a4 = c4*X(3);
    a5 = c5*X(2)*X(4);
    a6 = c6*X(5);
    a7 = c7*X(2)*X(5);
    a8 = c8*X(6);
    a9 = c9*X(1)*(X(1)-1)/2;
    a10 = c10*X(2);
    
    a = [a1; a2; a3; a4; a5; a6; a7; a8; a9; a10];
    a_sum = sum(a);
    
    tau=-1/a_sum*log(rand);
    u = find(cumsum(a)>a_sum*rand,1);
    X = X + D(:,u);
    
    T(iter) = t;
    N(:,iter) = X;
    t=t+ tau;
    iter = iter + 1;
    if iter>=Nalloc %adjust memory
        T = [T,zeros(1,Nalloc)];
        N = [N,zeros(Nchs,Nalloc)];
        Nalloc = 2*Nalloc;
    end
end
toc

%remove extra zeros
T = T(1:iter-1);
N = N(:,1:iter-1);   
t_m = T/60;

%figure 2 plots
figure
subplot(1,3,1)
plot(t_m,N(1,:))
axis([0 tmax/60 0 60])
% figure
subplot(1,3,2)
plot(t_m,N(2,:))
axis([0 tmax/60 0 80])
% figure
subplot(1,3,3)
plot(t_m,N(3,:))
axis([0 tmax/60 0 10])

%figure 3a plots
figure
subplot(1,3,1)
histogram(N(1,:),'Normalization','probability')
xlabel('Number of Monomers')
ylabel('Probability')
Y = get(gca,'YLim');
axis([0 45 0 max(.1,Y(2))])
% figure
subplot(1,3,2)
histogram(N(2,:),'Normalization','probability')
xlabel('Number of Dimers')
% ylabel('Probability')
Y = get(gca,'YLim');
axis([0 45 0 max(.06,Y(2))])
% figure
subplot(1,3,3)
histogram(N(3,:),'Normalization','probability')
xlabel('Number of mRNAs')
% ylabel('Probability')
Y = get(gca,'YLim');
axis([0 45 0 max(.25,Y(2))])
return