clear
close all
clc

profile on

sequenceLength = 1e2;
symbols = [-1 1];
M = length(symbols);
channelCoef = [1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2];
memory = length(channelCoef)-1;
noStates = M^memory;

%% Generate states
% initialize states
states = zeros(noStates,memory);
for i=1:memory
    temp = 1;
    for j=1:(M^(i-1))
        for k=1:M
            for m=1:(M^(memory-i))
                states(temp,i) = symbols(k);
                temp = temp + 1;
            end
        end
    end
end

%% Check possible transitions
nextState = zeros(1:memory-1);
currState = zeros(1:memory-1);
possibleTransition = zeros(noStates,noStates);

for j = 1:noStates      % future state loop
    for tmp1 = 2:memory
        nextState(tmp1-1) = states(j,tmp1);
    end
    for m = 1:noStates   % present state loop
        for tmp2 = 1:memory-1
            currState(tmp2) = states(m,tmp2);
        end
        %for loop over all possible paths from/to a state
        if (nextState==currState)   % Verifies trellis connection between state m and state j
            possibleTransition(j,m) = 1;
        end
    end
end    

bitSequence = randi([0 1], sequenceLength, 1);
TxSequence = bitToAntipodal(bitSequence, 1);

ISISequence = conv(channelCoef, TxSequence);
RxSequence = addAWGN(ISISequence, 6);

trellis = zeros(noStates,length(RxSequence)+1);
tmpx = zeros(noStates,1);
for i = 1:length(RxSequence)
    tic
    a = RxSequence(i);
    b = channelCoef(1);
    c = trellis(:,i);
    parfor j = 1:noStates      % future state loop
        temp1 = -2*states(j,1)*a;
        temp3 = states(j,1)^2*b;
        flag_first=1;

        v = zeros(noStates,1);

        for m = 1:noStates   % present state loop
            %for loop over all possible paths from/to a state
            if (possibleTransition(j,m))   % Verifies trellis connection between state m and state j
                temp2 = 0;
                for p = 1:memory
                    temp2 = states(m,p)*channelCoef(p+1) + temp2;
                end
                temp2 = states(j,1)*temp2;
                branch_metric = temp1+temp2+temp3;

                % test2 stores survivor sequences at time i   
                if flag_first==1
%                     trellis(j,i+1) = trellis(m,i)+branch_metric;
                    v(j) = c(m)+branch_metric;
                    test2(j,i)=m;
%                 elseif (trellis(j,i+1)>(trellis(m,i)+branch_metric))
                elseif (v(j)>(c(m)+branch_metric))
%                     trellis(j,i+1)=trellis(m,i)+branch_metric; % updates
%                     the best path merging into state j. Flag=1 no update needed, flag=2 update needed\
                    v(j) = c(m)+branch_metric;
                    test2(j,i)=m;
                end

                flag_first=0;
            end %if transition is possible
        end %second state loop
        tmpx(j) = v(j);
    end %first state loop
    trellis(:,i+1) = tmpx;
    toc
end  
    
profile viewer 

function [antipodalMap] = bitToAntipodal(sequence, energy)
    antipodalMap = zeros(length(sequence),1); 
    for i = 1:length(sequence)
        if sequence(i) == 1
            antipodalMap(i) = energy;   %map bit 1 to 1
        else
            antipodalMap(i) = -energy;  %map bit 0 to -1
        end
    end
end

function [noisySignal] = addAWGN(signal, SNRDB)
    N0=10.^-(SNRDB/10);     %PSD of the noise (linear)
    z = randn(length(signal),1);    %array of normal distributed noise
    noisySignal = signal+sqrt(N0/2)*z;      %adding the noise to the signal
end