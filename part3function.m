function [ber_bench1, bler_bench1, ber_bench2, bler_bench2, ber_bench3, bler_bench3] = part3function(z1, z2, z3)

load("OFDM_PILOT.mat"); % Load all necessary Matlab files
load("CODE.mat");
load("INTRLVR.mat");
load("MESG.mat");
load("ofdm_map.mat");
load("zw.mat");
load("Z_all.mat");

K = 2048;  % # of subcarriers
L = 199;   % # of taps - 1 for L1 = 0
KN = 112;   % # of null subcarriers
KP = 512;   % # of pilot subcarriers
Kd = 1420;  % # of data subcarriers
W = 21;    % # of OFDM symbols
% z1 = bb_rece_data_172648_1474; % First Hydrophone Signal
% z2 = bb_rece_data_172648_1475; % Second Hydrophone Signal
% z3 = bb_rece_data_172648_1476; % Third Hydrophone Signal
%z1 = Z_all;
%zw = Z_w(1:2048, 1);
dp = zeros(KP, 1); % pilot data
pilotind = zeros(KP, 1); % Index of pilot subcarriers
zp1 = zeros(KP, W); % pilot signal 1
zp2 = zeros(KP, W); % pilot signal 2
zp3 = zeros(KP, W); % pilot signal 3
zk1 = zeros(Kd, W); % Data signal 1
zk2 = zeros(Kd, W); % Data signal 2
zk3 = zeros(Kd, W); % Data signal 3
ncount = 1; % Counters
pcount = 1;
dcount = 1;
Hk1 = zeros(Kd, W); % Data only channel 1
Hk2 = zeros(Kd, W); % Data only channel 1
Hk3 = zeros(Kd, W); % Data only channel 1
zn1 = zeros(KN, W); % Null signal 1
zn2 = zeros(KN, W); % Null signal 2
zn3 = zeros(KN, W); % Null signal 3
Hw1 = zeros(K, W); % Whole signal channel 1
Hw2 = zeros(K, W); % Whole signal channel 2
Hw3 = zeros(K, W); % Whole signal channel 3

% Collect pilot and null subcarriers/indices
for k = 1:length(ofdm_map)% Null Index 
    if ofdm_map(k, 1) == 0
        zn1(ncount, :) = z1(k, 1:W);
        zn2(ncount, :) = z2(k, 1:W);
        zn3(ncount, :) = z3(k, 1:W);
        ncount = ncount + 1;
    elseif ofdm_map(k, 1) == 1 % Pilot Index 
        zp1(pcount, :) = z1(k, 1:W);
        zp2(pcount, :) = z2(k, 1:W);
        zp3(pcount, :) = z3(k, 1:W);
        dp(pcount, 1) = OFDM_PILOT(k, 1);
        pilotind(pcount, 1) = k;
        pcount = pcount + 1;
    end 
end

Ds = diag(dp); % Diagonalize pilot data
hls1 = zeros(L+1, W); % Time domain channel 1
hls2 = zeros(L+1, W); % Time domain channel 2
hls3 = zeros(L+1, W); % Time domain channel 3
l = 0:L;
V = exp(-1j*2*pi*(pilotind*l)/K); % estimating time channel
kv = transpose(0:K-1);
Vw = exp(-1j*2*pi*(kv*l)/K); % time to frequency matrix

for w = 1:W % Create Frequency ch
    % +annel from time channel 
    hls1(:, w) = (1/KP)*V'*Ds'*zp1(:, w);
    hls2(:, w) = (1/KP)*V'*Ds'*zp2(:, w);
    hls3(:, w) = (1/KP)*V'*Ds'*zp3(:, w);
    Hw1(:, w) = Vw*hls1(:, w);
    Hw2(:, w) = Vw*hls2(:, w);
    Hw3(:, w) = Vw*hls3(:, w);
end

for k = 1:1:length(ofdm_map) % Get data subcarriers/indices
    if ofdm_map(k, 1) == 2
        zk1(dcount, :) = z1(k, 1:W);
        zk2(dcount, :) = z2(k, 1:W);
        zk3(dcount, :) = z3(k, 1:W);
        Hk1(dcount, :) = Hw1(k, :);
        Hk2(dcount, :) = Hw2(k, :);
        Hk3(dcount, :) = Hw3(k, :);
        dcount = dcount + 1;
    end
end
Task 1b
twosigsquare1 = 1/KN * sum(abs(zn1).^2); % Noise Varience Estimation 1

x1 = 1/sqrt(2)+1i*1/sqrt(2);% QPSK Symbols
x2 = -1/sqrt(2)+1i*1/sqrt(2);
x3 = 1/sqrt(2)-1i*1/sqrt(2);
x4 = -1/sqrt(2)-1i*1/sqrt(2);

    normA = sqrt(conj(zk1 - Hk1*x1).*(zk1 - Hk1*x1));% LLR1 = log((A+C)/(B+D))
    normB = sqrt(conj(zk1 - Hk1*x2).*(zk1 - Hk1*x2));% LLR2 = log((A+B)/(C+D))
    normC = sqrt(conj(zk1 - Hk1*x3).*(zk1 - Hk1*x3));% Lines 101 - 117 accomplish this
    normD = sqrt(conj(zk1 - Hk1*x4).*(zk1 - Hk1*x4));
   
    A = -normA.^2 ./ twosigsquare1;
    B = -normB.^2 ./ twosigsquare1;
    C = -normC.^2 ./ twosigsquare1;
    D = -normD.^2 ./ twosigsquare1;

Lb1_zk = (max(A, C) + log(1+exp(-abs(C-A)))) - (max(B, D) + log(1+exp(-abs(D-B))));
Lb2_zk = (max(A, B) + log(1+exp(-abs(B-A)))) - (max(C, D) + log(1+exp(-abs(D-C))));

LR = zeros(W, 2*Kd); %  Combing both bits post LLR
for w = 1:W
    for n = 1:1:length(Lb1_zk)
        LR(w, 2*n-1) = Lb1_zk(n, w);
        LR(w, 2*n) = Lb2_zk(n, w);
    end
end
Task 1c
Name = './5G_LDPC_M10_N20_Z142_Q2_nonVer.txt'; % One tap equalizer
LPDC_Decoder = open('5G_LDPC_M10_N20_Z142_Q2_nonVer.txt');
[address, LPDC_INFOLEN] = ldpc_mex_initial_CAPI([1420, 2840, 2], Name);
est_code = zeros(length(LR), W);
LR_in_de = zeros(length(LR), W);
APP_code = zeros(length(LR), W);

ber = zeros (1, W);
bec = zeros (1, W);
bler = 0;
for w = 1:W
    LR_in_de(INTRLVR, w) = LR(w, :);
    APP_code(:, w) = ldpcDecoder_CAPI(address, LR_in_de(:, w));
    est_code(:, w) = (APP_code(:, w) < 0);
    bec(1, w) = sum(abs(est_code(:, w) - CODE(:, w)));
    ber(1, w) = bec(1, w) / length(CODE);
    if ber(1, w) ~= 0 && w>1 
        bler = bler + 1;
    end
end
bec_bench1 = sum(bec(1, 2:end));
ber_bench1 = sum(ber(1, 2:end))/(W-1);
bler_bench1 = bler/(W-1);
Task 2
%%% 2 Hydrophones %%%

twosigsquare2 = 1/KN * sum(abs(zn2).^2); % Noise Varience Estimation 2
normHk2 = zeros(Kd, W);
for w = 1:W
    for kd = 1:Kd
        normHk2(kd, w) = sqrt(Hk1(kd, w)'*Hk1(kd, w) + Hk2(kd, w)'*Hk2(kd, w));
    end
end
hermHk1 = Hk1';
hermHk2 = Hk2';
qk2 = zeros(Kd, W);
for w = 1:W
    for kd = 1:Kd
        qk2(kd, w) = [hermHk1(w, kd), hermHk2(w, kd)]*[zk1(kd, w); zk2(kd, w)];
    end
end
qk2 = qk2./normHk2;

sumSigSquared2 = twosigsquare1.*(conj(Hk1).*Hk1)+twosigsquare2.*(conj(Hk2).*Hk2);
twosigsquare2comb = 1./(normHk2.^2).*sumSigSquared2;

absA2 = abs(qk2 - normHk2.*x1).^2;
absB2 = abs(qk2 - normHk2.*x2).^2;
absC2 = abs(qk2 - normHk2.*x3).^2;
absD2 = abs(qk2 - normHk2.*x4).^2;

A2 = -absA2 ./ twosigsquare2comb;
B2 = -absB2 ./ twosigsquare2comb;
C2 = -absC2 ./ twosigsquare2comb;
D2 = -absD2 ./ twosigsquare2comb;

Lb1_zk2 = (max(A2, C2) + log(1+exp(-abs(C2-A2)))) - (max(B2, D2) + log(1+exp(-abs(D2-B2))));
Lb2_zk2 = (max(A2, B2) + log(1+exp(-abs(B2-A2)))) - (max(C2, D2) + log(1+exp(-abs(D2-C2))));

LR2 = zeros(W, 2*Kd); %  Combing both bits post LLR
for w = 1:W
    for n = 1:1:length(Lb1_zk2)
        LR2(w, 2*n-1) = Lb1_zk2(n, w);
        LR2(w, 2*n) = Lb2_zk2(n, w);
    end
end

est_code2 = zeros(length(LR2), W);
LR_in_de2 = zeros(length(LR2), W);
APP_code2 = zeros(length(LR2), W);

ber2 = zeros (1, W);
bec2 = zeros (1, W);
bler2 = 0;
for w = 1:W
    LR_in_de2(INTRLVR, w) = LR2(w, :);
    APP_code2(:, w) = ldpcDecoder_CAPI(address, LR_in_de2(:, w));
    est_code2(:, w) = (APP_code2(:, w) < 0);
    bec2(1, w) = sum(abs(est_code2(:, w) - CODE(:, w)));
    ber2(1, w) = bec2(1, w) / length(CODE);
    if ber2(1, w) ~= 0 && w>1 
        bler2 = bler2 + 1;
    end
end
bec_bench2 = sum(bec2(1, 2:end));
ber_bench2 = sum(ber2(1, 2:end))/(W-1);
bler_bench2 = bler2/(W-1);

%%% 3 Hydrophones %%%

twosigsquare3 = 1/KN * sum(abs(zn3).^2); % Noise Varience Estimation 3
normHk3 = zeros(Kd, W);
for w = 1:W
    for kd = 1:Kd
        normHk3(kd, w) = sqrt(Hk1(kd, w)'*Hk1(kd, w)+Hk2(kd, w)'*Hk2(kd, w)+Hk3(kd, w)'*Hk3(kd, w));
    end
end
hermHk3 = Hk3';
qk3 = zeros(Kd, W);
for w = 1:W
    for kd = 1:Kd
        qk3(kd, w) = [hermHk1(w, kd), hermHk2(w, kd), hermHk3(w, kd)]*[zk1(kd, w); zk2(kd, w); zk3(kd, w)];
    end
end
qk3 = qk3./normHk3;

sumSigSquared3 = twosigsquare1.*(conj(Hk1).*Hk1)+twosigsquare2.*(conj(Hk2).*Hk2)+twosigsquare3.*(conj(Hk3).*Hk3);
twosigsquare3comb = 1./(normHk3.^2).*sumSigSquared3;

absA3 = abs(qk3 - normHk3.*x1).^2;
absB3 = abs(qk3 - normHk3.*x2).^2;
absC3 = abs(qk3 - normHk3.*x3).^2;
absD3 = abs(qk3 - normHk3.*x4).^2;

A3 = -absA3 ./ twosigsquare3comb;
B3 = -absB3 ./ twosigsquare3comb;
C3 = -absC3 ./ twosigsquare3comb;
D3 = -absD3 ./ twosigsquare3comb;

Lb1_zk3 = (max(A3, C3) + log(1+exp(-abs(C3-A3)))) - (max(B3, D3) + log(1+exp(-abs(D3-B3))));
Lb2_zk3 = (max(A3, B3) + log(1+exp(-abs(B3-A3)))) - (max(C3, D3) + log(1+exp(-abs(D3-C3))));

LR3 = zeros(W, 2*Kd); %  Combing both bits post LLR
for w = 1:W
    for n = 1:1:length(Lb1_zk3)
        LR3(w, 2*n-1) = Lb1_zk3(n, w);
        LR3(w, 2*n) = Lb2_zk3(n, w);
    end
end

est_code3 = zeros(length(LR2), W);
LR_in_de3 = zeros(length(LR2), W);
APP_code3 = zeros(length(LR2), W);

ber3 = zeros (1, W);
bec3 = zeros (1, W);
bler3 = 0;
for w = 1:W
    LR_in_de3(INTRLVR, w) = LR3(w, :);
    APP_code3(:, w) = ldpcDecoder_CAPI(address, LR_in_de3(:, w));
    est_code3(:, w) = (APP_code3(:, w) < 0);
    bec3(1, w) = sum(abs(est_code3(:, w) - CODE(:, w)));
    ber3(1, w) = bec3(1, w) / length(CODE);
    if ber3(1, w) ~= 0 && w>1 
        bler3 = bler3 + 1;
    end
end
bec_bench3 = sum(bec3(1, 2:end));
ber_bench3 = sum(ber3(1, 2:end))/(W-1);
bler_bench3 = bler3/(W-1);

end