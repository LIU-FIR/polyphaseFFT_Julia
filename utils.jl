using FFTW

function pfft(x, win, N, R, OSR)
# Input:
#   x: L-length time squence, with multiple M-pts(N*R) long-FFT
#   win: M-pts window function，applying on M-pts(N*R) long-FFT
#   N: poly-phase FFT point
#   R: poly-phase branches
#   OSR: oversampled ratio，ranging from N/(N-1) to N
# Output：
#   X_PFB_frames: (N*frame_num) 2-d matrix，cols:polyphase-FFT，rows: frames in time   
    
    L = length(x)
    M = N * R # M-pts long-FFT
    win_regs = reshape(win, N, R)
    
    slide_len = Int(N/OSR)# compute slide-length，limited in 1~N-1
    frame_num = div(L - M, slide_len) + 1 #compute loaded-frames
    X_PFB_frames = zeros(ComplexF64,N,frame_num)

    for cnt = 1:frame_num
        # For each sliding, compute the head and tail of current pfft
        head_ptr = 1+(cnt-1)*slide_len
        tail_ptr = head_ptr + M-1
        #% load the current pfft into register.
        # xk_reg = x[head_ptr:tail_ptr]
        xk_pfb_regs = reshape(x[head_ptr:tail_ptr],N,R)
        X_PFB_frames[:,cnt] = fft(sum(win_regs .* xk_pfb_regs,dims=2))

    end
    return X_PFB_frames
end