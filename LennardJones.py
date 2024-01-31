import numpy as np 
from itertools import product

class LennardJones:

    def __init__(self, N, L, radius, duration, nsteps, v0):

        self.N =  N #number of particles 
        self.L = L #box length in units of sigma 
        self.radius = radius #radius of particles in units of sigma
        self.duration = duration #duration of simulation in units of tau
        self.nsteps = nsteps #number of steps 
        self.dt = duration/nsteps #small time step 
        self.v0 = v0 #initial velocity of particles in units of sigma/tau
        self.rc = 2.5 #cutoff radius in units of sigma 

        #initialise the positions, put the particles in a grid so they don't overlap
        grid_size = int(np.ceil(np.sqrt(N)))
        spacing = L/ grid_size 
        x = np.linspace(radius + spacing/2, L - radius - spacing/2, grid_size) 
        pos = list(product(x, x))
        self.positions = np.array(pos[:N])  #initial positions of particles

        #initialise the velocities, start with a fixed velocity magnitude and random directions 
        theta = np.random.uniform(0, 2*np.pi, size=N)
        vx,vy = v0*np.cos(theta), v0*np.sin(theta)
        self.velocities = np.stack((vx,vy), axis=1) #initial velocities of particles


    def pair_vector(self,r1,r2):
        """Calculate the vector r2 - r1 between two particles 
        with periodic boundary conditions."""
        r = r2 - r1
        r -= np.rint(r / self.L) * self.L
        return r 
    

    def force(self,r):
        """Calculate the force between two particles separated by a vector r."""

        dist = np.linalg.norm(r)
        if dist < self.rc:
            return (-48/dist**14 + 24/dist**8) * r #force is the gradient of the potential
        else:
            return np.zeros(2) # no force if particles are too far apart
        

    def get_forces(self):
        """Calculate the forces on all particles.
        To do this, we loop over all pairs of particles and sum the forces between them."""

        forces = np.zeros_like(self.positions)
        for i in range(self.N):
            for j in range(i+1,self.N):
                r = self.pair_vector(self.positions[i],self.positions[j]) #pair vector between particles i and j
                forces[i] += self.force(r) #Force felt by particle i due to particle j
                forces[j] -= self.force(r) 

        return forces
    

    def temperature(self, velocities):
        """Calculate the temperature of the system."""
        return 0.5*np.sum(velocities**2)/self.N
    

    def step(self):
        """Perform one time step of the simulation.""" 
        
        forces = self.get_forces() #calculate the forces on all particles at time t
        self.positions[:] = self.positions + self.velocities*self.dt + 0.5*forces*self.dt**2 #update the positions
        next_forces = self.get_forces() #calculate the forces on all particles at time t+dt
        self.velocities[:] = self.velocities + 0.5*(forces + next_forces)*self.dt #update the velocities


    def animate(self, T0=None, tau = 2): 
        """Evolve the system over nsteps integrations and return the positions and velocities at each time step"""

        all_positions = np.zeros((self.nsteps, self.N, 2))
        all_velocities = np.zeros_like(all_positions)
        for t in range(self.nsteps):
            all_positions[t] = self.positions
            all_velocities[t] = self.velocities
            self.step()

            if T0: 
                T = self.temperature(self.velocities)
                self.velocities *= np.sqrt(1+(T0/T-1)*self.dt/tau) #rescale the velocities to adjust the temperature


        return all_positions, all_velocities




    