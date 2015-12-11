function Emily_test
c = [0.043; % [s^-1] RNA -> RNA + M
	 0.0007; % [s^-1] M -> 0
	 0.0715; % [s^-1] DNA*D -> RNA + DNA*D
	 0.0039; % [s^-1] RNA -> 0
	 0.02; % [s^-1] DNA + D -> DNA*D
	 0.4791; % [s^-1] DNA*D -> DNA + D
	 0.002; % [s^-1] DNA*D + D -> DNA*2D
	 0.8765*(10^-11); % [s^-1] DNA*2D -> DNA*D + D
	 0.083; % [s^-1] M + M -> D
	 0.5]; % [s^-1] D -> M + M
    
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
Z = [0  0  0  0  0  0  0  0  0  0];
 
tmax = 600; %max time
%tmax = 2100; %If we are just recreating graph at 10 minutes, why go to 35 minutes? Since Markov mathematically same
mem = 10000; %memory allocation... what is math behind this?
spec = 6; %number of species
MC = 1000;
mem_all = mem*MC;
T_all = zeros(1,mem_all);
N_all = zeros(spec,mem_all);
iter = 1;

tic %Never knew about this command, love it
for i = 1:MC
    [T,N] = EE_Gillespie(tmax,c,X,mem,spec,Z,D);
    s = size(T,2);
    if iter+s>=mem_all
        [T_all,N_all,mem_all]=incr_mem(mem_all, T_all, N_all,spec);
    end
    T_all(1,iter:iter+s-1)=T;
    N_all(:,iter:iter+s-1)=N;
    iter = iter+s;
%     fprintf('Doubling memory for data. (Reallocation.)  t(%d)=%f', m, t(m))
end
toc
T_all = T_all(1:iter-1);
N_all = N_all(:,1:iter-1);   
t_m = T/60;

%figure 2 plots
figure('name','2 (gray) Exact Solution');
xlabel('Time (min)')
subplot(1,3,1)
plot(t_m,N(1,:))
axis([0 tmax/60 0 60])
ylabel('Number of Monomers')
% figure
subplot(1,3,2)
plot(t_m,N(2,:))
axis([0 tmax/60 0 80])
ylabel('Number of Dimers')
% figure
subplot(1,3,3)
plot(t_m,N(3,:))
axis([0 tmax/60 0 10])
ylabel('Number of mRNAs')

%figure 3a plots
figure('name','3a) Exact Solution');
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
axis([0 85 0 max(.06,Y(2))])
% figure
subplot(1,3,3)
histogram(N(3,:),'Normalization','probability')
xlabel('Number of mRNAs')
% ylabel('Probability')
Y = get(gca,'YLim');
axis([0 20 0 max(.25,Y(2))])

mem_all = mem*MC;
T_all = zeros(1,mem_all);
N_all = zeros(spec,mem_all);
iter = 1;

tic
for i = 1:MC
    [T,N] = EE_SQEA(tmax,c,X,mem,spec,Z);
    s = size(T,2);
    if iter+s>=mem_all
        [T_all,N_all,mem_all]=incr_mem(mem_all, T_all, N_all,spec);
    end
    T_all(1,iter:iter+s-1)=T;
    N_all(:,iter:iter+s-1)=N;
    iter = iter+s;
end
toc
T_all = T_all(1:iter-1);
N_all = N_all(:,1:iter-1); 


t_m = T/60;

%figure 2 plots
figure('name','7 (gray) SQEA');
xlabel('Time (min)')
subplot(1,3,1)
plot(t_m,N(1,:))
axis([0 tmax/60 0 60])
ylabel('Number of Monomers')
% figure
subplot(1,3,2)
plot(t_m,N(2,:))
axis([0 tmax/60 0 80])
ylabel('Number of Dimers')
% figure
subplot(1,3,3)
plot(t_m,N(3,:))
axis([0 tmax/60 0 10])
ylabel('Number of mRNAs')

%figure 3a plots
figure('name','3e) SQEA');
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
axis([0 85 0 max(.06,Y(2))])
% figure
subplot(1,3,3)
histogram(N(3,:),'Normalization','probability')
xlabel('Number of mRNAs')
% ylabel('Probability')
Y = get(gca,'YLim');
axis([0 20 0 max(.25,Y(2))])

return

function [T,N] = EE_Gillespie(tmax,c,X,mem,spec,Z,D) %Where do you update Z?
iter = 1;
T = zeros(1,mem); %empty time array
N = zeros(spec,mem); %empty species array
t=0;
% MC_num = 1;

while t <= tmax
    
    a1 = c(1)*(Z(3) - Z(4));
    a2 = c(2)*(2 + Z(1) - Z(2) - 2*Z(9) + 2*Z(10));
    a3 = c(3)*(Z(5) - Z(6) - Z(7) + Z(8));
    a4 = c(4)*(Z(3) -Z(4));
    a5 = c(5)*(4 - Z(5) + Z(6) - Z(7) + Z(8) + Z(9) - Z(10))*(2 - Z(5) + Z(6));
    a6 = c(6)*(Z(5) - Z(6) - Z(7) + Z(8));
    a7 = c(7)*(4 - Z(5) + Z(6) - Z(7) + Z(8) + Z(9) - Z(10))*(Z(5) - Z(6) - Z(7) + Z(8));
    a8 = c(8)*(Z(7) - Z(8));
    a9 = c(9)*(2 + Z(1) - Z(2) - 2*Z(9) + 2*Z(10))*(1 + Z(1) - Z(2) - 2*Z(9) + 2*Z(10))/2;
    a10 = c(10)*(4 - Z(5) + Z(6) - Z(7) + Z(8) + Z(9) - Z(10));
    
    a = [a1; a2; a3; a4; a5; a6; a7; a8; a9; a10];
    a_sum = sum(a);
    
    tau = exprnd(1/a_sum);
%    tau=1/a_sum*log(1/rand); %make sure rand is from exponential distribution, with exponential rate of the sum
%     tau=-1/a_sum*log(rand);
%     u = zeros(1,MC_num);
%     for i = 1:MC_num
%         u(i) = find(cumsum(a)>a_sum*rand,1);
%     end
%     u = mode(u);
    u = find(cumsum(a)>a_sum*rand,1); %pull from discrete distribution... Why don't you keep it normalized and do cumsum(a/a_sum)>rand?
    X = X + D(:,u);
%     X(1) = 2 + Z(1) - Z(2) - 2*Z(9) + 2*Z(10);
%     X(2) = 4 - Z(5) + Z(6) - Z(7) + Z(8) + Z(9) - Z(10);
%     X(3) = Z(3) -Z(4);
%     X(4) = 2 - Z(5) + Z(6);
%     X(5) = Z(5) - Z(6) - Z(7) + Z(8);
%     X(6) = Z(7) - Z(8);
  
    T(iter) = t;
    N(:,iter) = X;
    t=t+ tau;
    iter = iter + 1;
    if iter>=mem %adjust memory
%         T = [T,zeros(1,mem)];
%         N = [N,zeros(spec,mem)];
%         mem = 2*mem;
        [T,N,mem]=incr_mem(mem, T, N,spec);
    end
end
%remove extra zeros
T = T(1:iter-1);
N = N(:,1:iter-1); 
return

function [T,N,mem]=incr_mem(mem, T, N,spec)
    T = [T,zeros(1,mem)];
    N = [N,zeros(spec,mem)];
    mem = 2*mem;
return

function [T,N] = EE_SQEA(tmax,c,X,mem,spec,Z) %try doing fast first
iter = 1;
T = zeros(1,mem); %empty time array
N = zeros(spec,mem); %empty species array
t=0;
MC_num = 1;

while t <= tmax
    
    %slow reaction approximations
    A = 2 + Z(1) - Z(2) + c(10)/(2*c(9)) - (1/2);
    B = (1/4)*(2 + Z(1) - Z(2))*(2 + Z(1) - Z(2) - 1) ...
        - (c(10)/(2*c(9)))*(4 - Z(5) + Z(6) - Z(7) + Z(8));
    mu = .5*(A - sqrt((A*A)-4*B));
    
    %fast reaction approximations
    a1 = c(1)*(Z(3) - Z(4));
    a2 = c(2)*(2 + Z(1) - Z(2) - 2*mu);
    a3 = c(3)*(Z(5) - Z(6) - Z(7) + Z(8));
    a4 = c(4)*(Z(3) -Z(4));
    a5 = c(5)*(4 - Z(5) + Z(6) - Z(7) + Z(8) + mu)*(2 - Z(5) + Z(6));
    a6 = c(6)*(Z(5) - Z(6) - Z(7) + Z(8));
    a7 = c(7)*(4 - Z(5) + Z(6) - Z(7) + Z(8) + mu)*(Z(5) ...
        - Z(6) - Z(7) + Z(8));
    a8 = c(8)*(Z(7) - Z(8));
    
    a = [a1; a2; a3; a4; a5; a6; a7; a8];
    a_sum = sum(a);
    
    tau=1/a_sum*log(1/rand);
%     tau=-1/a_sum*log(rand);
    u = zeros(1,MC_num);
    for i = 1:MC_num
        u(i) = find(cumsum(a)>a_sum*rand,1);
    end
    u = mode(u);
    Z(u) = Z(u)+1;
    X(1) = 2 + Z(1) - Z(2) - 2*mu;
    X(2) = 4 - Z(5) + Z(6) - Z(7) + Z(8) + mu;
    X(3) = Z(3) -Z(4);
    X(4) = 2 - Z(5) + Z(6);
    X(5) = Z(5) - Z(6) - Z(7) + Z(8);
    X(6) = Z(7) - Z(8);
  
    T(iter) = t;
    N(:,iter) = X;
    t=t+ tau;
    iter = iter + 1;
    if iter>=mem %adjust memory
%         T = [T,zeros(1,mem)];
%         N = [N,zeros(spec,mem)];
%         mem = 2*mem;
        [T,N,mem]=incr_mem(mem, T, N,spec);
    end
end
%remove extra zeros
T = T(1:iter-1);
N = N(:,1:iter-1); 
return

% function [T,N] = EE_Poisson(tmax,c,X,mem,spec,D)
% iter = 1;
% T = zeros(1,mem); %empty time array
% N = zeros(spec,mem); %empty species array
% t=0;
% 
% while t <= tmax
%     
%     a1 = c(1)*X(3);
%     a2 = c(2)*X(1);
%     a3 = c(3)*X(5);
%     a4 = c(4)*X(3);
%     a5 = c(5)*X(2)*X(4);
%     a6 = c(6)*X(5);
%     a7 = c(7)*X(2)*X(5);
%     a8 = c(8)*X(6);
%     a9 = c(9)*X(1)*(X(1)-1)/2;
%     a10 = c(10)*X(2);
%     
%     a = [a1; a2; a3; a4; a5; a6; a7; a8; a9; a10];
%     a_sum = sum(a);
%     
% %     tau=-1/a_sum*log(rand);
%     u = find(cumsum(a)>a_sum*rand,1);
%     X = X + D(:,u);
%     
%     T(iter) = t;
%     N(:,iter) = X;
%     t=t+ .05;
%     iter = iter + 1;
%     if iter>=mem %adjust memory
%         T = [T,zeros(1,mem)];
%         N = [N,zeros(spec,mem)];
%         mem = 2*mem;
%     end
% end
% %remove extra zeros
% T = T(1:iter-1);
% N = N(:,1:iter-1); 
% return

