clear
close all
clc

tic

SNRDB = -4:2:20; %SNR in dB
SNR=10.^(SNRDB/10); %linear SNR (needed for qfunction)
energy = 1;
sequenceLength = 10000;
iterations = 10;
avBER = 0;
avBERv = 0;
channelMemory = [1, 0.4, -0.2];
states = [1,1;-1,1;1,-1;-1,-1];

bitSequence = randomBitSequence(sequenceLength);
antipodalSequence = bitToAntipodal(bitSequence, energy);
ISISequence = conv(channelMemory, antipodalSequence);
%ISISequence = [1, -1.2,  0.5, -1.5, -0.2, 1, 0.8, 0.9];
seqLen = length(ISISequence);

for j = 1:iterations %start simulation and run iterations number of time
    j
    bitSequence = randomBitSequence(sequenceLength);   %generate random array of 1 and 0
    antipodalSequence = bitToAntipodal(bitSequence, energy); 
    ISISequence = conv(channelMemory, antipodalSequence);
    
    for i = 1:length(SNRDB)
        m = addAWGN(ISISequence, SNRDB(i));
        mhatv = viterbi(m, states, channelMemory);
        decodedBits = antipodalToBit(mhatv, energy);
        
        mhat = estimateSignal(m, energy);       %estimate signal from recieved values
        bhat = antipodalToBit(mhat, energy);    %convert estimates into bits sequence

        BER(i) = calcultateBitErrorRate(bitSequence, bhat);
        BERv(i) = calcultateBitErrorRate(bitSequence, decodedBits);
        
%         bit = transpose(bitSequence)
%         vit = transpose(decodedBits)
%         MLD = transpose(bhat)
    end
    
    avBER = avBER + BER;
    avBERv = avBERv + BERv; %sum the BER of all iterations
end


avBER = avBER / iterations;
avBERv = avBERv / iterations;

%% Plot
figure
semilogy(SNRDB,avBER,'-o');
axis([-4 20 10^-6 10^0])
hold on
semilogy(SNRDB,avBERv,'-x')
grid on
legend('ML', 'viterbi');
xlabel('SNR [dB]')
ylabel('BER')

toc

%% Reciever viterbi
function [decodedAntipodal] = viterbi(m, states, channelMemory)
    U(:,1) = initialUValue(m, states, channelMemory);
    prevState = zeros(4,1);

    for i = length(channelMemory):length(m)
        V = getV(states, m, i, channelMemory);
        [U(:,i-1), prevState(:,i-2)] = getUNext(states, U(:,i-2), V);
    end

    k = shortestPath(length(m), prevState, U);
    decodedAntipodal = transpose(isiToAntipodal(states, k));
end

%% functions 
function [bhat] = antipodalToBit(sequence, energy)
    bhat = zeros(length(sequence),1);
    for i = 1:length(sequence)
        if (sequence(i) == energy)
            bhat(i) = 1;            %if signal is 1 bit is 1
        else
            bhat(i) = 0;            %if signal is -1 bit is 0
        end
    end
end

function [k] = shortestPath(seqLen, prevState, U)
    k = zeros(1,seqLen-1);
    [a,b] = min(U(:,seqLen-1));
    for i = 1:seqLen-1
        if (i == 1)
            k(seqLen-1) = b;%2 for example in book
        else
            k(seqLen-i) = prevState(k(seqLen-i+1),seqLen-i);
        end
    end
end

function [decodedBits] = isiToAntipodal(states, k)
    decodedBits = zeros(1,length(k)+1);
    for i = 1:length(k)
        if (i == 1)
            decodedBits(i) = states(k(i), 2);
            decodedBits(i+1) = states(k(i), 1);
        else
            decodedBits(i+1) = states(k(i), 1);
        end
    end
end

function [sequence] = randomBitSequence(length)
   sequence = randi([0 1], length, 1);      %generate random 1's and 0's
end

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
    N0=sqrt(10.^-(SNRDB/10)/2);     %PSD of the noise (linear)
    z = randn(length(signal),1);    %array of normal distributed noise
    noisySignal = signal+N0*z;      %adding the noise to the signal
end

function [U] = initialUValue(ISISequence, states, channelMemory)
    temp1 = 0;
    temp2 = 0;
    for i = 1:(length(channelMemory)-1)
        temp1 = ISISequence(i).*states(:, 2-i+1) + temp1;
        for j = 1:2
            temp2 = states(:,i).*states(:,j).*channelMemory(abs(i-j)+1) + temp2;
        end
    end
    U = -2*temp1+temp2;
end

function [V] = getV(states, ISISequence, recievedSymbol, channelMemory)
    V = zeros(4:4);    
    
    for currstate = 1:length(states)  

        temp1 = -2*states(currstate,1)*ISISequence(recievedSymbol);
        temp3 = states(currstate,1)^2*channelMemory(1); 
        for prevstate = 1:length(states)
            temp2 = 0;
            if (states(currstate, 2) == states(prevstate,1))
           
                for i = 1:(length(channelMemory)-1)
                    temp2 = temp2 + states(prevstate, i)*channelMemory(i+1);
                end
                temp2 = 2*states(currstate,1)*temp2;
                V(prevstate, currstate) = temp1+temp2+temp3; 
            else
                V(prevstate, currstate) = NaN; 
            end
        end
    end
end

function [Unext, prevState] = getUNext(states, U, V)
    test = U+V;
    [Unext, prevState] = min(test);
end

function [BER] = calcultateBitErrorRate(inputSequence, outputSequence)
    errorCount = 0;
    for i = 1:length(inputSequence)
        if inputSequence(i) ~= outputSequence(i)
            errorCount = errorCount + 1;        %sum the total ammount of errors
        end
    end
    
    BER = errorCount / length(inputSequence);   %number of errors/total number of bits = error rate
end

function [mhat] = estimateSignal(signal, energy)
    mhat = zeros(length(signal),1);
    for i = 1:length(signal)
        temp1 = abs(signal(i) - energy);        %distance to 1
        temp2 = abs(signal(i) + energy);        %distance to -1
        if temp1 < temp2
            mhat(i) = energy;                   %if distance to 1 is smaller mhat = 1
        else
            mhat(i) = -energy;                  %if distance to -1 is smaller mhat = -1
        end
    end
end