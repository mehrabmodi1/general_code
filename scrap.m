clear all
close all

%generating kernel
kernel = zeros(1,10);
kernel(3) = 5;
kernel(4:8) = [4:-1:0];

%input signal
signal = rand(1, 10000);

%convolved output
conv_sig = conv(signal, kernel);
conv_sig = conv_sig + rand(1, 10009).*0.01;

deconv_kernel = deconv(conv_sig, signal);

plot(kernel, 'r.')
hold on
plot(kernel, 'r')
plot(deconv_kernel, 'b')

% figure(2)
% plot(conv_sig)