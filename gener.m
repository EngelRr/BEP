sequenceLength = 1e3;    % Sequence length per iteration 
symbols = [-1 1];
M = length(symbols);
channelCoef = [1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2];
memory = length(channelCoef)-1;
noStates = M^memory;

tic
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
toc
tic

nextState = zeros(1:memory-1);
currState = zeros(1:memory-1);
possibleTransition = zeros(noStates,noStates);

parfor j = 1:noStates      % future state loop
    v = zeros(1:memory-1);
    for tmp1 = 2:memory
        v(tmp1-1) = states(j,tmp1);
    end
    nextState = v;
    for m = 1:noStates   % present state loop
        w = zeros(1:memory-1);
        for tmp2 = 1:memory-1
            w(tmp2) = states(m,tmp2);
        end
        currState = w;
        %for loop over all possible paths from/to a state
        if (nextState==currState)   % Verifies trellis connection between state m and state j
            possibleTransition(j,m) = 1;
        end
    end
end
toc