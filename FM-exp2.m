clc;
clear all;
close all;

Fs_pcm = 50e3;
f_msg = 1e4;
L = 8;
A = 1;
Fs_sdr = 1e6;
T = 0.05;

sps = Fs_sdr / Fs_pcm;
n = (0:1/Fs_pcm:T-1/Fs_pcm);

msg = A * sin(2 * pi * f_msg * n);

delta = 2*A/L;
q_index = floor((msg + A)/ delta);
q_index(q_index >= L) = L-1;
q_signal = -A + delta*(q_index + 0.5);

pcm_bits = de2bi(q_index, log2(L), "left-msb");
pcm_stream = pcm_bits;
pcm_stream = pcm_stream(:);

tx_bits = 2 * pcm_stream;
tx_upsampled = upsample(tx_bits, sps);
tx_bb = complex(tx_upsampled, zeros(size(tx_upsampled)));

tx = sdrtx("Pluto");
tx.CenterFrequency = 2.45e9;
tx.BasebandSampleRate = Fs_sdr;
tx.Gain = 0;

rx = sdrrx("Pluto");
rx.CenterFrequency = tx.CenterFrequency;
rx.BasebandSampleRate = Fs_sdr;
rx.SamplesPerFrame = length(tx_bb);
rx.GainSource = "AGC Slow Attack";
rx.OutputDataType = "double";

tx(tx_bb);
rx_bb = rx();

rx_sym = downsample(real(rx_bb), sps);
rx_sym = -rx_sym;
rx_bits = rx_sym > 0;

rx_bits = rx_bits(1:length(pcm_stream));
rx_bits = reshape(rx_bits, [], log2(L));
rx_index = bi2de(rx_bits, 'left-msb');

rx_msg = -A + delta * (rx_index + 0.5);

[b, a] = butter(4, (2*f_msg)/(Fs_pcm/2));
rx_filtered = filter(b, a, rx_msg);
Nplot = round(1e-3*Fs_pcm);
N = min([Nplot, numel(n), numel(q_signal), numel(rx_msg), numel(rx_filtered)]);
t1 = n(1:N);
figure
subplot(4,1,1);
plot(t1*1e3 , msg(1:Nplot), 'c' , 'LineWidth' , 1.2);
title('original analog signal');
xlabel('time(ms)');
ylabel('Amplitude');
grid on

subplot(4,1,2);
plot(t1*1e3, q_signal(1:Nplot) , 'm' , 'LineWidth', 1.2);
title('quantized signal');
xlabel('time(ms)');
ylabel('Amplitude');
grid on

subplot(4,1,3);
plot(t1*1e3, rx_msg(1:N) , 'y' , 'LineWidth', 1.2);
title('reconstructed signal (without filtering)');
xlabel('time(ms)');
ylabel('Amplitude');
grid on

subplot(4,1,4);
plot(t1*1e3, rx_filtered(1:N) , 'b' , 'LineWidth', 1.2);
title('reconstructed signal (after LPF)');
xlabel('time(ms)');
ylabel('Amplitude');
grid on
