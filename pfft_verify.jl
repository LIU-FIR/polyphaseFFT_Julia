using PyPlot, DSP, FFTW
include("utils.jl")
##
L = 200000; # 原始序列长度.the length of the original time sequence
N = 256; # pFFT-length
R = 4; # Taps in each polyphase branch. 
M = N * R; # pFFT block size
n = 0:L-1;   # sampling index (time)
k = -N/2:N/2-1; # FFT-bin index (normalized frequency) 
srate = 1024e6; # sampling rate.

oversampled_rate = 1

# K = div(L-N*R,slide_len)
# @show cnt_reload, K

ham_win = hamming(M); # hamming window function
responsetype = Lowpass(1 / N);
designmethod = FIRWindow(rect(M)); # use DSP.jl design sinc_win function
sinc_win = digitalfilter(responsetype, designmethod); # sinc window function
sinc_win_norm = sinc_win ./ maximum(sinc_win) # 归一化sinc window function
sys_win = sinc_win_norm .* ham_win # composed system window function

fs_pfb_cs = srate / N # PFB critically-sampled rate
fs_pfb_os = srate / N * oversampled_rate # PFB over-sampled rate

# generate a sequence of 3 tones, 
# v0 is at the multiple-integer frequency
# v1 and v2 are shifted from v0.
# v2 is configured such to aliase to v3
# (v0+v1+v2) has the same effect with (v0+v1+v3), except v2 is filtered in PFB.
kν = 16
ν0 = srate * kν / N
ν1 = ν0 + fs_pfb_cs / 4
ν2 = ν0 + 3fs_pfb_cs / 4
ν3 = ν0 - fs_pfb_cs / 4

tL = 0:L-1
x_ν0 = exp.(2im * pi * ν0 / srate .* tL)
x_ν1 = 2 * exp.(2im * pi * ν1 / srate .* tL)
x_ν2 = 3 * exp.(2im * pi * ν2 / srate .* tL)
x_ν3 = 3 * exp.(2im * pi * ν3 / srate .* tL)

x = x_ν0 + x_ν1 + x_ν2
# x = x_ν0 + x_ν1 + x_ν3
# x = x_ν0 + x_ν1
# x = x_ν0 


X_pfft_fra = pfft(x, sys_win, N, R, oversampled_rate)
# 
fig1, ax1 = subplots()
ax1.plot(abs.(X_pfft_fra[:, 1]), "b.-")

N_zft = 512
kz = -N_zft/2:N_zft/2-1

Xkν_PFB_fra = X_pfft_fra[kν+1, 1:512]

Yν_zft = fft(Xkν_PFB_fra)

fig2, ax2 = subplots()
ax2.plot(kz, fftshift((abs.(Yν_zft))), "k-")