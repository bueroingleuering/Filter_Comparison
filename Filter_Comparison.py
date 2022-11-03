"""@author: robin leuering"""
import numpy as np
import matplotlib.pyplot as plt
"""Plot Function"""
def compare_plot(t,u,y,F,U,Y,title):
    plt.figure(); plt.rcParams['figure.figsize'] = [16, 12]
    plt.plot(t,u,'k',label="u_org")
    plt.plot(t,y,'b',label="y_filter")
    plt.xlim([t[0],t[-1]]), plt.xlabel('t [s]',fontsize=16)
    plt.title('Signal time domain '+title,fontsize=20)
    plt.legend(fontsize=16), plt.grid(True)
    plt.savefig(title+"time.svg")
    
    plt.figure(); plt.rcParams['figure.figsize'] = [16, 12]
    plt.stem(F,U,'k', markerfmt='ko',label="U_org")
    plt.stem(F,Y,'b', markerfmt='bo',label="Y_filter")
    plt.xlabel('f [Hz]',fontsize=16)
    plt.title('Signal freq. domain '+title,fontsize=20)
    plt.legend(fontsize=16), plt.grid(True)
 
plt.close('all')
"""Signal Parameter """
A   = 100                                           # Amplitude
f   = 50                                            # Fundamental frequency
w   = 2*f*np.pi
N   = 128                                           # Number of samples
K   = int(N/2+1)                                    # Wave number
dt  = 1/(N*f)                                       # Time difference
t   = np.arange(0,1/f,dt)                           # Time array
F   = np.arange(0,K)*f                              # Frequency array
u   = A*np.sin(w*t+0*np.pi)                         # Fundamental sine wave
u   = u + (A/5)*np.sin(20*w*t)                      # 20th harmonic
u   = u + (A/10)*np.sin(50*w*t)                     # 50th harmonic
y   = 0*u
"""Filter Parameter """
fc  = 10*f                                          # Cut off frequency
tau = 1/(2*np.pi*fc)                                # Time constant
[a1,a2]  = [dt/(dt+tau),tau/(dt+tau)]               # Filter coefficients
"""Filter Algorithm"""
for i in range(1,N):
    y[i]   = a1*u[i]+a2*y[i-1]
"""Filter Gain Compensation"""
y = y*np.sqrt(1+pow(w*tau,2))                       # Fdt damping compensation
"""FFT Filter Matrix"""
H = np.eye(N)
[H[20][20],H[N-20][N-20]] = [0,0]                   # Set U[20*f] to zero
[H[50][50],H[N-50][N-50]] = [0,0]                   # Set U[50*f] to zero
"""FFT"""
Y = 2*np.fft.fft(y/N)
U = 2*np.fft.fft(u/N)
Z = np.matmul(H,U)
"""inverse FFT"""
z = (N*np.fft.ifft((Z))/2).real                     # Amplitude gain
"""Plot"""
compare_plot(t, u, y, F, np.abs(U[0:K]),np.abs(Y[0:K]), "Low Pass")    
compare_plot(t, u, z, F, np.abs(U[0:K]),np.abs(Z[0:K]), "FFT/iFFT") 