### **Plastic corrector algorithm** ###
#### Author: Abhishek Palchoudhary ####
############### 25/07/24 ##############

import numpy as np
import math
import matplotlib
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy import optimize
import time as tt
font = {'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

class load_func:

    def __init__(self, amplitude, period):
        self.amplitude = amplitude
        self.period = period

    def eval(self, t):
        a = self.amplitude
        p = self.period
        return ((4*a/p)*np.abs(((t-p/4)%(p))-p/2)-a)

def loading_function(n_cycles, amplitude, num_time_steps, end_time):

    period = end_time / n_cycles  # Calculate period based on number of cycles
    time = np.linspace(0, end_time, num_time_steps)
    load = load_func(amplitude, period).eval(time)

    peaks, _ = find_peaks(load, height=0)
    valleys, _ = find_peaks(-load, height=0)
    peaksvalleys = np.concatenate([peaks, valleys])
    peaksvalleys = np.sort(peaksvalleys)

    plt.figure()
    plt.plot(time, load, 'k')
    plt.ylabel(r'$f(t)$')
    plt.xlabel('Time (s)')
    plt.grid()
    plt.show()

    return time, load, peaksvalleys

def Newton_vectorized(f, x0, tol=1e-3, max_iter=1000, eps=1.0e-5, verbose=False):
    converged = True
    x = x0.copy()  # Make a copy of x0 to avoid modifying the original array
    f_cur = f(x0)
    f_cur_ini = f_cur
    i = 0

    while not np.all(f_cur < 1e-10):
        i += 1

        # Compute the gradient
        f_pert = f(x + eps)
        f_prime = (f_pert - f_cur) / eps

        # Update only the elements of x where f_cur is positive (i.e. update points outside yield surface)
        positive_indices = f_cur > 0
        x[positive_indices] -= f_cur[positive_indices] / f_prime[positive_indices]

        # Recompute the function values for the updated x
        f_cur = f(x)

        if verbose:
            print("Newton iteration", i, "/ f_cur:", f_cur, "/ f_ini:", f_cur_ini)

        if i >= max_iter + 1 or np.any(np.isnan(f_cur)):
            print("Newton not converged")
            x = np.full_like(x0, np.nan)
            converged = False
            break

    return x, converged

class visco_plast:

  def initialise(self, n_points, e_p_0, s_0, e_0, E, nu, sig_y, b, Q, C, D, x_0, sig_vm_e, p_reinitialize):

    self.e_p_0 = e_p_0          #Passing initialized e_p (explicitly specified as array of zeros)
    self.e_p_n = e_p_0
    self.s_0 = s_0
    self.e_0 = e_0
    self.x_o = x_0
    self.x_n = x_0

    self.E = E
    self.nu = nu
    self.sig_y = sig_y

    self.b = b
    self.Q = Q
    self.C = C
    self.D = D

    self.mu = self.E/(2*(1+self.nu))
    self.sig_vm_e = sig_vm_e

    self.p = p_reinitialize

  def residual(self,e_p_n,e_p_0,delta_t,g_t,e_p_o,p_o,s_0,x_o): # returns residual calculated for the first step, non-incremental Neuber

    s_n = self.s_new(e_p_n, e_p_0, g_t, s_0)

    p_n = self.p_new(e_p_n,e_p_o,p_o)

    R = self.Q * (1-np.exp(-self.b*p_n))

    x_n = self.x_new(e_p_n,e_p_o,p_n,p_o,x_o)

    sig_VM = np.abs(s_n - x_n/(2*self.mu))*self.sig_vm_e

    f = sig_VM - self.sig_y - R

    return f

  def x_new(self,e_p_n,e_p_o,p_n,p_o,x_o):
    return (x_o+2.0/3.0*self.C*(e_p_n - e_p_o)) / (1.0 + self.D*(p_n - p_o))

  def p_new(self,e_p_n,e_p_o,p_o):
    return p_o + np.abs( e_p_n - e_p_o ) * (1/(3*self.mu)) * self.sig_vm_e

  def s_new(self,e_p_n,e_p_0,g_t,s_0):
    discriminant = (e_p_n-e_p_0)**2  +  4*g_t
    if self.point%2==0:
      s_new =( ( -(e_p_n-e_p_0)  +  np.sqrt(discriminant) ) / ( 2 ) ) + s_0
    elif self.point%2==1:
      s_new =( ( -(e_p_n-e_p_0)  -  np.sqrt(discriminant) ) / ( 2 ) ) + s_0
    return s_new

  def e_new(self,e_p_n,e_p_0,s_n,s_0,e_0):
    return (s_n-s_0) + (e_p_n-e_p_0) + e_0

  def integrate_backward_euler(self,n_points,time,load,peaksvalleys):

    self.x_n_save = np.zeros((len(time),n_points))
    self.e_p_save = np.zeros((len(time),n_points))
    self.s_n_save = np.zeros((len(time),n_points))
    self.e_n_save = np.zeros((len(time),n_points))
    self.p_save = np.zeros((len(time),n_points))

    self.sig_VM_save = np.zeros((len(time),n_points))

    self.e_p_save[0,:] = self.e_p_0
    self.p_save[0,:] = self.p
    self.s_n_save[0,:] = self.s_0
    self.e_n_save[0,:] = self.e_0
    self.newtontime=np.zeros((len(time),n_points))

    e_p_o = np.zeros(n_points)

    f_to = load[0]
    time_o = time[0]
    self.point=0

    ############ Time integration #############
    for i in range(1,len(time)):

      f_t = load[i]

      if i in peaksvalleys:   ######### Change of peaks for load reversal ###########
        f_to = load[peaksvalleys[self.point]]
        self.e_p_0=self.e_p_n
        self.s_0=self.s_n
        self.e_0=self.e_n
        self.point+=1

      delta_t = time[i] - time[i-1]

      g_t = pow((f_t-f_to),2)

      p_o = self.p
      x_o = self.x_n
      e_p_o=self.e_p_n

      start = tt.process_time()

      fun = lambda e_p_n: self.residual(e_p_n,self.e_p_0,delta_t,g_t,e_p_o,p_o,self.s_0,x_o)

      if np.any(np.array(fun(self.e_p_n))>0):
        self.e_p_n, _ = Newton_vectorized(fun, e_p_o)

      test = fun(self.e_p_n)
      if np.any(np.array(test)>1.0e-5):
        print('Time instant', i)
        print("test",test)

      self.newtontime[i,:] = tt.process_time() - start

      self.s_n=self.s_new(self.e_p_n, self.e_p_0, g_t, self.s_0)
      self.e_n=self.e_new(self.e_p_n,self.e_p_0,self.s_n,self.s_0, self.e_0)
      self.p = self.p_new(self.e_p_n,e_p_o,p_o)
      self.x_n = self.x_new(self.e_p_n,e_p_o,self.p,p_o,x_o)

      self.e_p_save[i,:] = self.e_p_n   #Saving at time t=1,2,...len(time), for n_points
      self.s_n_save[i,:] = self.s_n
      self.e_n_save[i,:] = self.e_n
      self.p_save[i,:] = self.p
      self.x_n_save[i,:] = self.x_n

      #Von Mises stress
      self.sig_VM_save[i,:] = np.abs(self.s_n - self.x_n/(2*self.mu))*self.sig_vm_e

  def plot(self,time):
    plt.plot(time,self.p_save[:,:])
    plt.grid()
    plt.ylabel(r'Cumulative plastic strain p')
    plt.xlabel('Time')
    plt.show()

    plt.plot(self.e_n_save[:,:],self.s_n_save[:,:])
    plt.grid()
    plt.ylabel('s')
    plt.xlabel('e')
    plt.show()
    return

  def extract_qtys(self,time):
    print('Time taken for time integration (seconds):',np.sum(self.newtontime, axis=0)[0])
    return self.p_save[:,:], np.squeeze(self.e_p_save[:,:]), self.s_n_save[:,:], self.e_n_save[:,:]
