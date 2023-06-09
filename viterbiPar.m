
function [symbols] = viterbiPar(RxSequence,states,channelCoef,possibleTransition)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

memory = length(channelCoef)-1;
noStates = size(states,1);

trellis = zeros(noStates,length(RxSequence)+1);
test2 = zeros(noStates,length(RxSequence));

% f = parallel.pool.Constant(states);
% disp(f.Value);

for i = 1:length(RxSequence)
    %for loop over number of states
    a = RxSequence(i);
    b = channelCoef(1);
    d = channelCoef(2:end);
    g = trellis(:,i);
    tmpx =[];
    parfor j = 1:noStates      % future state loop
        temp1 = -2*states(j,1)*a;
        temp3 = states(j,1)^2*b;
        flag_first=1;
        
        c = zeros(noStates,1);
        
        for m = 1:noStates   % present state loop
            %for loop over all possible paths from/to a state
            if (possibleTransition(j,m))   % Verifies trellis connection between state m and state j
                %temp2 = sum(f.Value(m,:).*d);
                temp2 = sum(states(m,:).*d);
%                 temp2 = 0;
%                 for p = 1:memory
%                     temp2 = states(m,p)*channelCoef(p+1) + temp2;
%                 end
                temp2 = states(j,1)*temp2;
                branch_metric = temp1+temp2+temp3;
                   
                % test2 stores survivor sequences at time i   
                if flag_first==1
                    c(j) = g(m)+branch_metric;
                    test2(j,i)=m;
                elseif (c(j)>(g(m)+branch_metric))
                    c(j)=g(m)+branch_metric; % updates the best path merging into state j. Flag=1 no update needed, flag=2 update needed
                    test2(j,i)=m;
                end

                flag_first=0;
            end 
        end
        tmpx(j) = c(j);
    end
    trellis(:,i+1) = tmpx;
end

k = zeros(1,length(RxSequence)+1);
test = length(RxSequence)+1;
for i = 1:test
    if (i == 1)
        [~,k(test)] = min(trellis(:,end));
    else
        k(test-i+1) = test2(k(test-i+2),size(test2,2)+1-i+1);
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