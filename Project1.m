clear;
% 1. Load speech signal
Fs = 4000;
[mSpeech, Fs] = audioread("MaleSpeech-16-4-mono-20secs.wav");
% sound(mSpeech,Fs)
% Consider the speech signal in 1.5s
t = 0:1/Fs:1.5;
plot(t, mSpeech(1:length(t)), 'LineWidth', 2);
hold on

% 2. Quantize the sample signal
L = 16; % the number of quantization levels
V_p = 0.5625; % the peak voltage of signal
% Determine the single quantile interval ?-wide
q = 2 * V_p / (L - 1); % Use the exact equation
quantized_pos = q/2:q:V_p;
quantized_neg = -q/2:-q:-V_p;
q = sort([quantized_neg, quantized_pos]);
s_q_2 = quan_uni(mSpeech(1:length(t)), q); % Uniform quantization

% Plot the sample signal and the quantization signal
plot(t, s_q_2, 'ro', 'MarkerSize', 6, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

% 3. Calculate the average quantization noise power,
% the average power of the sample signal and SNR
e_uni = mSpeech(1:length(t)) - s_q_2; % error between sample signal and quantized signal

pow_noise_uni = 0;
pow_sig = 0;
    for i = 1:length(t)
        pow_noise_uni = pow_noise_uni + e_uni(i)^2;
        pow_sig = pow_sig + mSpeech(i)^2;
    end
    pow_noise_uni = pow_noise_uni / length(t);
    pow_sig = pow_sig / length(t);
    SNR = pow_sig / pow_noise_uni;
% --------compression-------------

% 5. Compress the sample signal ‘mSpeech’
mu = 255; %mu-Law
% A = 87.6; %use the standard value A-Law
y_max = V_p;
x_max = V_p;
% Replace the compress equation for u-law and A-law
% with x is the 'mSpeech' signal
s_c_5 = ulaw(x_max, y_max, mSpeech(1:length(t)), mu); %mu-Law
% s_c_5 = Alaw(x_max, y_max, mSpeech(1:length(t)), A); %A-Law
% Plot the compress signal;
plot(t, s_c_5, 'r--');

% 6. Quantize the compress signal and plot the quantized signal
s_q_6 = quan_uni(s_c_5, q);
plot(t, s_q_6, 'b^', 'MarkerSize', 6, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');

% 7. Expand the quantized signal
s_e_7 = inv_ulaw(x_max, y_max, s_q_6, mu); %mu-Law
% s_e_7 = inv_Alaw(x_max, y_max, s_q_6, A); %A-Law

plot(t, s_e_7, 'g*', 'MarkerSize', 6, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');

legend('Sample signal', 'Uniform quantized values', 'Compress signal', ...
    'Compress quantized values', 'Nonuniform quantized values');

% 9. Calculate the average quantization noise power,
% the average power of the analog signal and SNR
e_com = mSpeech(1:length(t)) - s_e_7;
pow_noise_com = 0;

for i = 1:length(t)
    pow_noise_com = pow_noise_com + e_com(i)^2;
end
pow_noise_com = pow_noise_com / length(t);
SNR_a_com = pow_sig / pow_noise_com;

%%Function
function x = quan_uni(a, combined_quantized)
    x = zeros(size(a));
    for i = 1:length(a)
        for t = 1:length(combined_quantized)-1
            if (a(i) >= combined_quantized(t) && a(i) <= combined_quantized(t + 1))
                ave = (combined_quantized(t) + combined_quantized(t + 1))/2;
                if (a(i) < ave)
                    x(i) = combined_quantized(t);
                else 
                    x(i) = combined_quantized(t + 1);
                end
            end
        end
    end
end

%Function compand mu-law
function y = ulaw(xmax, ymax, a, mu)
    y = zeros(size(a));   
    for i = 1:length(a)
        x = a(i);
        if x >= 0
            sign_x = 1;
        else
            sign_x = -1;
        end
        y(i) = ymax * (1/log(1+mu))*log(1+mu*abs(x)/xmax)*sign_x;
    end
end

%Function expand mu-law
function x = inv_ulaw(xmax, ymax, y, mu)
    x = zeros(size(y));
    for i = 1:length(y)
        n = y(i); 
        if n >= 0
            sign_x = 1;
        else
            sign_x = -1;
        end
        x(i) = sign_x * ((exp(abs(n) * ...
        (log(1 + mu)/ymax)) - 1) * (xmax/mu));
    end
end

%Function compand A-law
function y = Alaw(xmax, ymax, a, A)
    y = zeros(size(a));   
    for i = 1:length(a)
        x = a(i);
		temp = abs(a(i))/xmax;
        if x >= 0
            sign_x = 1;
        else
            sign_x = -1;
        end
		if (0 < temp && temp <= 1/A)
			y(i) = ymax * sign_x * (A * abs(x) / xmax) / (1 + log(A));
		elseif (1/A < temp && temp < 1)
			y(i) = ymax * sign_x * (1 + log(A * abs(x) / xmax)) / (1 + log(A));
        end
    end
end

%Function expand A-law
function x = inv_Alaw(xmax, ymax, y, A)
    x = zeros(size(y));
    for i = 1:length(y)
        n = y(i);
        if n >= 0
            sign_x = 1;
        else
            sign_x = -1;
        end
        x(i) = sign_x * abs(n) * (1 + log(A)) * xmax / ( ymax * A);
        temp = abs(x(i))/xmax;
        if (1/A < temp && temp < 1)
            x(i) = sign_x * exp(abs(n)*(1+log(A))/(ymax)-1) * (xmax / A);     
        end
    end
end


