#3.2.1 program
from scipy.integrate import RK45, solve_ivp
import numpy as np 
import matplotlib.pyplot as plt

def plot_trajectory(statvector):
    
    p0 = np.array([192200, 192200, 0])#km
    pd = np.array([-153760, 0, 0]) #desired position
    plt.plot(statvector[:,0], statvector[:,1])
    plt.plot([p0[0]], [p0[1]],marker='o',markersize=20, markerfacecolor='green')
    plt.plot([pd[0]], [pd[1]],marker='o',markersize=20, markerfacecolor='red')

    plt.show()

def rk4_step( f, t, y, h ):
	'''
	Calculate one RK4 step
	'''
	k1 = f( t, y )
	k2 = f( t + 0.5 * h, y + 0.5 * k1 * h )
	k3 = f( t + 0.5 * h, y + 0.5 * k2 * h )
	k4 = f( t +       h, y +       k3 * h )

	return y + h / 6.0 * ( k1 + 2 * k2 + 2 * k3 + k4 )


def dynamics(t, state):
    '''
    circular restricted 3 body problem 
    given in EQ 2.23
    '''
    rho = state[:3]
    v = state[3:]
    m_moon = 7.34767e22 #kg
    m_earth = 5.97219e24 #kg
    mu= m_moon/(m_earth + m_moon)
#    lstar = #distance between two primaries r12
    d= np.sqrt((rho[0]+mu)**2 + rho[1]**2 + rho[2]**2)#r13/lstar
    r= np.sqrt((rho[0]-1+mu)**2 + rho[1]**2 + rho[2]**2)#r23/lstar
    return np.array([v[0],
                     v[1],
                     v[2],
        2*v[1] + rho[0] - (1-mu)*(rho[0]+mu)/d**3 - mu*(rho[0]-1+mu)/r**3,
        -2*v[0] + rho[1] - (1-mu)*rho[1]/d**3 - mu*rho[1]/r**3,
        -rho[2]*(1-mu)/d**3 - mu*rho[2]/r**3
        ])

def propogate(state, t, dt, tf):
    '''
    propogate dynamics given initial state at inital time
    to a final time
    '''
#    lstar = 238900 * 1.609
#    G = 6.674e-20
#    m_moon = 7.34767e22 #kg
#    m_earth = 5.97219e24 #kg
#    mstar = (m_earth + m_moon)
#    tstar = np.sqrt(lstar**3 /(G*mstar) )
#    h = 5/tstar
#    x = solve_ivp(dynamics, (t,t+ h), state ).y
#    print('ivp_solution',x)
#    return x[:]
    #while t<tf:
#    solution= RK45(dynamics, t, state,tf ) 
#    while True:
#        solution.step()
#        if solution.status == 'finished':
#            break
#        t_values = solution.t
#        y_values=solution.y
#    y_values=np.array(y_values)
#    print('y value shape',y_values.shape)
#    return y_values 
    states = np.zeros((int(tf//dt),6))
    states[0] = np.copy(state) 
    i=1
    t+=dt 
    while t<tf:
    #for i in range(int(t),int(tf),int(dt)):
        if i>=(int(tf//dt)):
            break
        states[i] = rk4_step(dynamics,t,states[i-1],dt )
        t+=dt
        i+=1
#    print(states)
    return states 


if __name__ == "__main__":
    lstar = 238900 * 1.609
    G = 6.674e-20
    m_moon = 7.34767e22 #kg
    m_earth = 5.97219e24 #kg
    mstar = (m_earth + m_moon)
    tstar = np.sqrt(lstar**3 /(G*mstar) )
    #initial conditions given in 3.31 and 3.32
    p0 = np.array([192200, 192200, 0]) / lstar #km
    v0=np.array([-.5123, .1025, 0]) * tstar/ lstar #km
    tf = 4.3425 * 24 * 60 * 60/tstar #4.3425 days in seconds
    pd = np.array([-153760, 0, 0])/lstar #desired position
    
    dt =10/lstar#.5 

    #tolerance for answer
    epsilon = 3.844e-3 # km
    #constraint vector F=p1-p1d = 0
    X = np.array([v0]) #design vector
    
    state_guess = np.array([p0[0], p0[1], p0[2], v0[0], v0[1], v0[2]])
    prop = propogate(state_guess, 0, dt, tf)
    plot_trajectory(prop*lstar)
    print(prop)
    residual = np.linalg.norm(prop[-1][:3] - pd)* lstar 
    

    print(f'p1:: {prop[-1][:3]}, pd: {pd}')
    while residual > epsilon:
        print('residual', residual)
        rho1 = state_guess[:3]
        v_guess = state_guess[3:]
        print(f'v_guess: {v_guess}')
        
        #compute jacobian using finite difference
        h=.001* tstar / lstar
        b1 = ((propogate(state_guess + np.array([0,0,0,h, 0,0]),0,dt,tf )[-1] - propogate(state_guess - np.array([0,0,0,h, 0,0]),0,dt,tf ))[-1] / (2*h))[:3]
        b2 = ((propogate(state_guess + np.array([0,0,0,0, h,0]),0,dt,tf )[-1] - propogate(state_guess - np.array([0,0,0,0, h,0]),0,dt,tf ))[-1] / (2*h))[:3]
        b3 = ((propogate(state_guess + np.array([0,0,0,0, 0,h]),0,dt,tf )[-1] - propogate(state_guess - np.array([0,0,0,0, 0,h]),0,dt,tf ))[-1] / (2*h))[:3]

        B = np.vstack((
            b1,
            b2,
            b3
            )).T
        print('b:', B) 
        DX = np.matmul(-np.linalg.inv(B),rho1-pd)# rho1-pd*0)
        state_guess[3:] += DX*.07 
        print(f'DX: {DX}')
        p1 = propogate(state_guess, 0, dt, tf)
        
        plot_trajectory(p1*lstar)
        
        print(f'p1:: {p1[-1][:3]}, pd: {pd}')
        
        residual = np.linalg.norm( p1[-1][:3]- pd) * lstar 
#        print('iteration',residual )