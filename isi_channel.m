clear
close all
clc
%


%Example 7.13
channelMemory = [1, 0.4, -0.2];     %given [S_0, S_(+-1), S_(+-2)]
states = [1,1;-1,1;1,-1;-1,-1];     %given [(a_(l-1), a_(l-2)]
%states = [  1,1,1 ; -1,1,1 ; 1,-1,1 ; -1,-1,1; 1,1,-1; -1,1,-1; 1,-1,-1; -1,-1,-1];
prevState = zeros(length(states),1);             %initialize

ISISequence = [1, -1.2,  0.5, -1.5, -0.2, 1, 0.8, 0.9]; %recieved sequence [Z_0 ... Z_7]
seqLen = length(ISISequence);   %length of recieved sequence

U(:,1) = initialUValue(ISISequence, states, channelMemory);     %step 1

for i = length(channelMemory):seqLen
    V = getV(states, ISISequence, i, channelMemory);    %step 2
    %the matrix U stores all the state metrics with states in vertical and
    %step in the horizontal 
    %the matrix prevState stores all the previous states of each state with
    %the vertical being the state and the value being the previous state of
    %that node in the trellis
    [U(:,i-length(channelMemory)+2), prevState(:,i-length(channelMemory)+1)] = getUNext(U(:,i-length(channelMemory)+1), V); %step 3
end
%seqLen = seqLen - length(channelMemory)+2;
k = shortestPath(seqLen, prevState, U); %step %5
decodedBits = isiToAntipodal(states, k);

%% functions 
function [k] = shortestPath(seqLen, prevState, U)
    k = zeros(1,seqLen-1);  %initialize
    [a,b] = min(U(:,seqLen-1)); %get the lowest end state of the trellis
    for i = 1:seqLen-1
        if (i == 1)
            k(seqLen-1) = b;    %2 for example in book
        else
            k(seqLen-i) = prevState(k(seqLen-i+1),seqLen-i);    %trace back the path taken through the trellis 
        end
    end
end

function [decodedBits] = isiToAntipodal(states, k)
    %state consist of sigma_k = (a_k, a_(k-1)) where a_k referes to a +1 or
    %-1 
    decodedBits = zeros(1,length(k)-1);
    for i = 1:length(k)
        if (i == 1)
            decodedBits(i) = states(k(i), 2);
            decodedBits(i+1) = states(k(i), 1);
        else
            decodedBits(i+1) = states(k(i), 1);
        end
    end
end

%intialUValue is used in the first step of the viterbi algorithm to
%determine the values of the states in the first step.
function [U] = initialUValue(ISISequence, states, channelMemory)
    temp1 = 0;
    temp2 = 0;
    for i = 1:(length(channelMemory)-1)     %sum from 0 to L-1 (L =0,1,2)-> matlab 1 to L-1 (L=1,2,3)
        temp1 = ISISequence(i).*states(:, length(channelMemory)-i) + temp1;   %first term Z_i * a_i 
        for j = 1:(length(channelMemory)-1)   %sum from 1 to L-1 (L=1,2,3) used for second term in 7.119
            temp2 = states(:,i).*states(:,j).*channelMemory(abs(i-j)+1) + temp2; %second term a_i*a_j*s(i-j)
        end
    end
    U = -2*temp1+temp2;     %7.119 -2*first term + second term
end

function [V] = getV(states, ISISequence, recievedSymbol, channelMemory)
    V = zeros(4:4);     %initialize
    
    %nextstate refers to sigma_(k+1)
    %currstate refers to sigma_(k) = (a_(k-1), a(k-2))
    %state consists of 1 or -1 and another 1 or -1 in case of example 7.13
    for currstate = 1:length(states)   %for each possible state (in this case 4)
        %all terms are multiplied by Z_k where k is the variable
        %recievedSymbol this will be filtered out later 
        temp1 = -2*states(currstate,1)*ISISequence(recievedSymbol); %first term in 7.120 -2*a_k*Z_k
        temp3 = states(currstate,1)^2*channelMemory(1); %third term (a_k)^2*s_0
        for prevstate = 1:length(states)
            temp2 = 0;
            %the if statements ensure only possible state transitions are
            %taken into account
            if (states(currstate, 2) == states(prevstate,1))
                %second term in 7.120 2*a_k*(sum_(m=k-L)^(k-L)a_m*s_(k-m)) k = L = 3
                for i = 1:(length(channelMemory)-1)
                    temp2 = temp2 + states(prevstate, i)*channelMemory(i+1);
                end
                temp2 = 2*states(currstate,1)*temp2;
                V(prevstate, currstate) = temp1+temp2+temp3; %all possible branches of the trellis
            else
                V(prevstate, currstate) = NaN;  %impossible branches are entered NaN in matrix
            end
        end
    end
end

function [Unext, prevState] = getUNext(U, V)
    %sum the values at each state with their possible trellis branches
    test = U+V;
    %pick the next state with the lowest weight
    %Unext stores the weight of a node
    %prevState stores the previous state of a node
    [Unext, prevState] = min(test);
end