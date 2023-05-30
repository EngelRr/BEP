function [symbols] = viterbi(RxSequence,symbols,channelCoef)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M = length(symbols);
memory = length(channelCoef)-1;
noStates = M^memory;

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

trellis = zeros(noStates,length(RxSequence)+1);
test2 = zeros(noStates,length(RxSequence));

for i = 1:length(RxSequence)
    %for loop over number of states
    for j = 1:noStates      % future state loop
        temp1 = -2*states(j,1)*RxSequence(i);
        temp3 = states(j,1)^2*channelCoef(1);
        flag_first=1;
        
        for m = 1:noStates   % present state loop
            %for loop over all possible paths from/to a state
            if (states(j,memory)==states(m,1))   % Verifies trellis connection between state m and state j
                temp2 = 0;
                for p = 1:memory
                    temp2 = states(m,p)*channelCoef(p+1) + temp2;
                end
                temp2 = states(j,1)*temp2;
                branch_metric = temp1+temp2+temp3;
                   
                % test2 stores survivor sequences at time i   
                if flag_first==1
                    trellis(j,i+1) = trellis(m,i)+branch_metric;
                    test2(j,i)=m;
                elseif (trellis(j,i+1)>(trellis(m,i)+branch_metric))
                    trellis(j,i+1)=trellis(m,i)+branch_metric; % updates the best path merging into state j. Flag=1 no update needed, flag=2 update needed
                    test2(j,i)=m;
                end

                flag_first=0;
            end
            
        end
    end
end

k = zeros(1,length(test2)+1);
test = length(test2)+1;
for i = 1:test
    if (i == 1)
        [~,k(test)] = min(trellis(:,length(trellis)));
    else
        k(test-i+1) = test2(k(test-i+2),length(test2)+1-i+1);
    end
end

% decoding the Rx
decodedSyms = zeros(length(k)+1,1);
for i = 1:length(k)
    if (i == 1)
        for j = 1:memory
            decodedSyms(j) = states(k(i), memory-j+1);
        end
    else
        decodedSyms(i+memory-1) = states(k(i),1);
    end
end
symbols = decodedSyms;
end