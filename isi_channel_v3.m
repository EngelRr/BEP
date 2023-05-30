clear
close all
clc
tic
rng(1);
profile on

%% Parameters
sequenceLength = 1e3;
errorCount = 0;
memory = 2;
symbols = [-1 1];
M = length(symbols);
SNRDB = 4:1:8; %SNR in dB
%SNR=10.^(SNRDB/10); %linear SNR 
energy = 1;
noStates = M^memory;


%% Generate states
%initialize states
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

% Get random sequence
bitSequence = randi([0 1], sequenceLength, 1);
TxSequence = bitToAntipodal(bitSequence, energy);
channelCoef = [1,0.2,0.1];
ISISequence = conv(channelCoef, TxSequence);

BER = zeros(1,length(SNRDB));
BERref = zeros(1,length(SNRDB));
BERSxS = zeros(1,length(SNRDB));



parfor SNRDBvalue = 1:length(SNRDB)
    disp(SNRDB(SNRDBvalue));
    %errorCount and loopCount used to get more data for plot
    errorCount = 0;
    loopCount = 0;
    while (errorCount < 1e3)
        loopCount = loopCount+1;
        
        %initializing variables that are used later to compare if a
        %connection between states is possible
        nextState = zeros(1:memory-1);
        currState = zeros(1:memory-1);
        option = zeros(M,2);

        %creating the RxSequence
        RxSequence = addAWGN(ISISequence, SNRDB(SNRDBvalue));
        %% Viterbi
        %initialize trellis
        trellis = zeros(noStates,length(RxSequence));
        path = zeros(noStates,length(RxSequence)-1);  

        %for loop over length of sequence
        for i = 1:length(RxSequence)
            %for loop over number of states
            for j = 1:noStates
                %x used for indexing the option array later on
                x = 1;
                %calculate path metrics according to path metric =
                %2*a_k*Z_k + a_k*sum(a_m*s_(k-m)) + a_k^2*s_0 = temp1 +
                %temp2 + temp3
                temp1 = -2*states(j,1)*RxSequence(i);
                temp3 = states(j,1)^2*channelCoef(1);
                
                %selecting the overlap of the next state with the previous
                %state, used to check for possible transitions. It is done
                %in this way to improve speed as compared to selecting
                %indexes from array. Same is done with currState later on.
                for test2 = 2:memory
                    nextState(test2-1) = states(j,test2);
                end

                for m = 1:noStates
                    %for loop over all possible paths from/to a state
                    %see nextState for explanation
                    for test3 = 1:memory-1
                        currState(test3) = states(m,test3);
                    end
                    %check for possible state transitions
                    if(states(j,memory)==states(m,1))
                        temp2 = 0;
                        for p = 1:memory
                            temp2 = states(m,p)*channelCoef(p+1) + temp2;
                        end
                        temp2 = states(j,1)*temp2;
                        total = temp1+temp2+temp3;
                        %all possible options of state transitions and
                        %the previous state are stored in option
                        option(x,1) = trellis(m,i)+total;
                        option(x,2) = m;
                        x = x+1;
                    end
                end
                
                %determine the transition which leads to the lowest overall
                %metric. The value is stored in trellis, the previous state
                %is stored in path.
                minimum = option(1,1);
                minimumPos = 1;
                for optieLen = 1:length(option)  
                    if (option(optieLen,1) <= minimum)
                        minimum = option(optieLen,1); 
                        minimumPos = optieLen;
                    end
                end
                trellis(j,i+1) = option(minimumPos,1);
                path(j,i) = option(minimumPos,2);

    %             [~,b] = min(optie);
    %             trellis(j,i+1) = optie(b(1),1);
    %             path(j,i) = optie(b(1),2);
            end
        end
    % end of trellis processing 

    % Final sequence selection 
        k = zeros(1,length(path)+1);
        test = length(path)+1;
        for i = 1:test
            if (i == 1)
                [~,k(test)] = min(trellis(:,length(trellis)));
            else
                k(test-i+1) = path(k(test-i+2),length(path)+1-i+1);
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

        %aligning bits
        recievedBits = antipodalToBit(decodedSyms);
        alignedBitSeq = recievedBits(memory+1:length(recievedBits)-memory);

        % Comparison with ISI-free AWGN performance 
        reference = addAWGN(TxSequence, SNRDB(SNRDBvalue));
        referenceRecieved = estimateSignal(reference, energy);
        referenceDecoded = antipodalToBit(referenceRecieved);

        %Comparison with SxS 
        SxSRecieved = estimateSignal(RxSequence, energy);
        SxSDecoded = antipodalToBit(SxSRecieved);

        %Calculating the bit error rate
        BER(SNRDBvalue) = sum(bitSequence ~= alignedBitSeq)/length(bitSequence) +BER(SNRDBvalue);
        BERref(SNRDBvalue) = sum(bitSequence ~= referenceDecoded)/length(bitSequence) + BERref(SNRDBvalue);
        BERSxS(SNRDBvalue) = sum(bitSequence ~= SxSDecoded(1:length(bitSequence)))/length(bitSequence) + BERSxS(SNRDBvalue);

        errorCount = errorCount + sum(bitSequence ~= alignedBitSeq);
    end
    
    %getting the average BER for every loop
    BER(SNRDBvalue) = BER(SNRDBvalue)/loopCount;
    BERref(SNRDBvalue) = BERref(SNRDBvalue)/loopCount;
    BERSxS(SNRDBvalue) = BERSxS(SNRDBvalue)/loopCount;
end

profile viewer
toc
%% Plot
figure
semilogy(SNRDB,BER,'-o');
hold on;
semilogy(SNRDB,BERref,'-+');
semilogy(SNRDB,BERSxS,'-x')
PE2=qfunc(sqrt(2*10.^(SNRDB/10)));     %theoretical BER qfunc(sqrt(2Eb/N0)) Eb/N0 is linear SNR
semilogy(SNRDB,PE2,'-');
axis([-4 16 10^-6 10^0])
grid on
legend('viterbi', 'reference','SxS');
xlabel('SNR [dB]')
ylabel('BER')

%% Function
function [bhat] = antipodalToBit(sequence)
    bhat = zeros(length(sequence),1);
    for i = 1:length(sequence)
        if (sequence(i) == 1)
            bhat(i) = 1;            %if signal is 1 bit is 1
        else
            bhat(i) = 0;            %if signal is -1 bit is 0
        end
    end
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