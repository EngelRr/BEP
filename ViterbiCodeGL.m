clear
close all
clc
tic
profile on
%rng('default')

%% Parameters
sequenceLength = 1e3;    % Sequence length per iteration 
symbols = [-1 1];
M = length(symbols);
channelCoef = [1,0.2,0.1];
memory = length(channelCoef)-1;
noStates = M^memory;


%iterations = 10;
SNRDB = 4:1:8; %SNR in dB
SNR=10.^(SNRDB/10); %linear SNR 
avBER = 0;
avBERref = 0;
avBERSxS = 0;
energy = 1;

TargetErrorCount=1e3; 

%% Variable Initialisation
BER = zeros(1,length(SNRDB));
BERref = zeros(1,length(SNRDB));
BERSxS = zeros(1,length(SNRDB));

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


   
for asdf = 1:length(SNRDB)
    BitErrorCount=0;
    BitErrorCountRef=0;
    BitErrorCountSxS=0;
    it=0;
    
    disp(['SNR=' num2str(SNRDB(asdf)) 'dB']);

    while BitErrorCount<TargetErrorCount
        %% Get random sequence
        bitSequence = randi([0 1], sequenceLength, 1);
        TxSequence = bitToAntipodal(bitSequence, energy);
        %channelCoef = [1,normrnd(0,.3,1,memory)];

        ISISequence = conv(channelCoef, TxSequence);
        RxSequence = addAWGN(ISISequence, SNRDB(asdf));
        %RxSequence=ISISequence;
        %ISISequence = [1, -1.2,  0.5, -1.5, -0.2, 1, 0.8, 0.9];
        %% Viterbi
        %initialize trellis
        trellis = zeros(noStates,length(RxSequence)+1);
        %trellis(:,1) = [2.4;-3.2;5.6;3.2];
        path = zeros(noStates,length(RxSequence));
        pathMetrics = zeros(1,M);

        %for loop over length of sequence
        for i = 1:length(RxSequence)
            %for loop over number of states
            BestState=zeros(1,noStates);

            for j = 1:noStates      % future state loop

                %BestState=1;        % State index indicating the originary state for the best path into state j at time i
                flag_first=1;
                for m = 1:noStates   % present state loop
                    %for loop over all possible paths from/to a state
                    if (states(j,memory)==states(m,1))   % Verifies trellis connection between state m and state j
                        a_vec=[states(j,1),states(m,1),states(m,memory)];   % contains input symbols at time a_k, a_k-1, ....a_k-M
                        temp1 = -2*a_vec(1)*RxSequence(i);
                        temp3 = channelCoef(1);
                        temp2 = 0;
                        for k = 1:memory  % This should be the memory length
                            temp2 = a_vec(k+1)*channelCoef(k+1) + temp2;
                        end
                        temp2 = temp2*a_vec(1);
                        branch_metric = temp1+temp2+temp3;
                        if flag_first==1
                            trellis(j,i+1) = trellis(m,i)+branch_metric;
                            BestState(j)=m;
                        elseif (trellis(j,i+1)>(trellis(m,i)+branch_metric))
                            trellis(j,i+1)=trellis(m,i)+branch_metric; % updates the best path merging into state j. Flag=1 no update needed, flag=2 update needed
                            BestState(j)=m;
                        end

                        flag_first=0;
                    end
                end
            end

            % Stores survivor sequences at time i
            if i==1
                path(:,i)=states(:,1);
            else
                path(:,1:i)=[path(BestState,1:i-1),states(:,1)];    % Updates surviving path in state j at trellis section i
            end
        end
        % end of trellis processing

        % Final sequence selection
        [MinMetric,SurvPathState]=min(trellis(:,i));

        % decoding the Rx
        %decodedSyms = zeros(length(k)+1,1);
        decodedSyms = path(SurvPathState,:);

        %aligning bits
        recievedBits = antipodalToBit(decodedSyms(1:end-memory));
        alignedBitSeq = recievedBits;

        % Comparison with ISI-free AWGN performance
        reference = addAWGN(TxSequence, SNRDB(asdf));
        referenceRecieved = estimateSignal(reference, energy);
        referenceDecoded = antipodalToBit(referenceRecieved);

        %Comparison with SxS
        SxSRecieved = estimateSignal(RxSequence, energy);
        SxSDecoded = antipodalToBit(SxSRecieved((1:end-memory)));
        BitErrorCount=sum(bitSequence~=alignedBitSeq)+BitErrorCount;
        BitErrorCountRef=sum(bitSequence~=referenceDecoded)+BitErrorCountRef;
        BitErrorCountSxS=sum(bitSequence~=SxSDecoded)+BitErrorCountSxS;
        it=it+1;
    end
        BER(asdf) =  BitErrorCount/(it*length(bitSequence));
        BERref(asdf) = BitErrorCountRef/(it*length(bitSequence));
        BERSxS(asdf) = BitErrorCountSxS/(it*length(bitSequence));
end
profile viewer
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
legend('Viterbi', 'ISI-free','SxS');
xlabel('SNR [dB]')
ylabel('BER')

%toc
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
    N0=10.^-(SNRDB/10);     %PSD of the noise (linear)
    z = randn(length(signal),1);    %array of normal distributed noise
    noisySignal = signal+sqrt(N0/2)*z;      %adding the noise to the signal
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