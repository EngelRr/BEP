clear
close all
clc
tic

%rng('default')

%% Parameters
sequenceLength = 1e5;
memory = 2;
symbols = [-1 1];
M = length(symbols);
iterations = 10;
SNRDB = -4:2:16; %SNR in dB
SNR=10.^(SNRDB/10); %linear SNR 
avBER = 0;
avBERref = 0;
avBERSxS = 0;
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


   
    
%% Get random sequence
bitSequence = randi([0 1], sequenceLength, 1);
TxSequence = bitToAntipodal(bitSequence, energy);
%channelCoef = [1,normrnd(0,.3,1,memory)];
channelCoef = [1,0.4,-0.2];
ISISequence = conv(channelCoef, TxSequence);

BER = zeros(1,length(SNRDB));
BERref = zeros(1,length(SNRDB));
BERSxS = zeros(1,length(SNRDB));

for asdf = 1:length(SNRDB)
    disp(asdf);
    RxSequence = addAWGN(ISISequence, SNRDB(asdf));
    %ISISequence = [1, -1.2,  0.5, -1.5, -0.2, 1, 0.8, 0.9];
    %% Viterbi
    %initialize trellis
    trellis = zeros(noStates,length(RxSequence));
    %trellis(:,1) = [2.4;-3.2;5.6;3.2];
    path = zeros(noStates,length(RxSequence)-1);
    pathMetrics = zeros(1,M);

    %for loop over length of sequence
    for i = 1:length(RxSequence)
        %for loop over number of states
        for j = 1:noStates
            x = 1; 
            for m = 1:noStates
                %for loop over all possible paths from/to a state
                if(states(m,1:size(states,2)-1) == states(j,2:size(states,2)))
                temp1 = -2*states(j,1)*RxSequence(i);
                temp3 = states(j,1)^2;
                temp2 = 0;
                    for k = 1:M
                        temp2 = states(m,k)*channelCoef(k+1) + temp2;
                    end
                temp2 = 2*states(j,1)*temp2;
                total = temp1+temp2+temp3;
                optie(x,:) = [trellis(m,i)+total,m];
                x = x+1;
                end
            end
            [~,b] = min(optie);
            trellis(j,i+1) = optie(b(1),1);
            path(j,i) = optie(b(1),2);
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
    alignedBitSeq = recievedBits(memory+1:length(recievedBits)-2);


    % Comparison with ISI-free AWGN performance 
    reference = addAWGN(TxSequence, SNRDB(asdf));
    referenceRecieved = estimateSignal(reference, energy);
    referenceDecoded = antipodalToBit(referenceRecieved);
    
    %Comparison with SxS 
    SxSRecieved = estimateSignal(RxSequence, energy);
    SxSDecoded = antipodalToBit(SxSRecieved);

    BER(asdf) = sum(bitSequence + alignedBitSeq)/length(bitSequence);
    BERref(asdf) = sum(bitSequence + referenceDecoded)/length(bitSequence);
    BERSxS(asdf) = sum(bitSequence + SxSDecoded(1:length(bitSequence)))/length(bitSequence);

%     BER(asdf) = sum(bitSequence ~= alignedBitSeq)/length(bitSequence);
%     BERref(asdf) = sum(bitSequence ~= referenceDecoded)/length(bitSequence);
%     BERSxS(asdf) = sum(bitSequence ~= SxSDecoded(1:length(bitSequence)))/length(bitSequence);
end

%% Plot
figure
semilogy(SNRDB,BER,'-o');
hold on;
semilogy(SNRDB,BERref,'-+');
semilogy(SNRDB,BERSxS,'-x')
PE2=qfunc(sqrt(2*SNR));     %theoretical BER qfunc(sqrt(2Eb/N0)) Eb/N0 is linear SNR
semilogy(SNRDB,PE2,'-');
axis([-4 20 10^-6 10^0])
grid on
legend('viterbi', 'reference','SxS');
xlabel('SNR [dB]')
ylabel('BER')

toc
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