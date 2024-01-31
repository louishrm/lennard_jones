# Interacting particles and the Lennard-Jones potential in Python

# Particles interacting in a Lennard Jones potential

In this notebook, we'll simulate particles interacting in a Lennard-Jones potential. The steps are going to be very similar to those we followed for the ideal gas simulation, so I really recommend
you go see [this video](https://www.youtube.com/watch?v=Gn5QoDCDGgo&t=12s) first. 

The only difference is that we'll use **periodic boundary conditions** and we have to change our method for integrating the equations of motion since there are now **forces** between the particles. 

## Assumptions

* The particles are represented by circles that all have the same radius and mass. 
* The particles start at random positions and with random velocities. 
* The particles can only interact with each other through the Lennard-Jones potential. 
* We use periodic boundary conditions: if a particle leaves the box, it appears on the other side. 


## Algorithm 
The idea behind the algorithm to run the simulation is simple: 

* Fix the number of particles, the size of the container and length of simulation.
* Initialize the positions and velocities at time $t=0$. 
* At each time step $t$, find the forces acting on all the particles. 
* Numerically integrate the equations of motion to evolve the positions and velocities from $t \to t+ \Delta t$.



## 1 The Lennard-Jones potential

The Lennard-Jones (LJ) potential between two particles is given by: 

$$U_{LJ}(\vec r) = 4 \epsilon \left[ \left( \frac{\sigma}{r}\right)^{12} - \left( \frac{\sigma}{r}\right)^{6} \right],$$
where $r = |\vec r|$ is the distance between the two particles. We will use reduced units $\epsilon \to 1$, $\sigma \to 1$ , i.e work in units of $\epsilon$ and $\sigma$. Since $U_{LJ}$ rapidly approaches $0$ as $r$ increases, we usually introduce a cutoff $r_c$ such that $U_{LJ}(r > r_c) = 0.$ For this simulation, we'll take $r_c = 2.5$ units of $\sigma$.
The term in $1/r^{12}$ captures short range **repulsive interactions** between the particles, whilst the term in $1/r^6$ captures long range attractive forces between the particles. There is a stable equilibrium corresponding to the minimum of $U_{LJ}$. You can read more about the LJ potential [here](https://lamma.engineering.unt.edu/sites/default/files/class5_handout_mtse_5010_2018.pdf).

## 2 Initial conditions
Let's initialize the variables we need for this simulation and set $\vec r(0)$ and $\vec v(0)$. This will pretty much be exactly the same than for the ideal gas video. 

## 3. Resulting force 

The LJ potential induces a force that is proportional to the gradient of $U_{LJ}$: 

$$\vec F( \vec r)  = - \frac{\text{d} U_{LJ}}{\text{d} r} \vec u_r = \left(\frac{24}{r^7} -\frac{48}{r^{13}} \right) \frac{\vec r}{r} = \left(\frac{24}{r^8} -\frac{48}{r^{14}} \right) \vec r$$


This is the force that appears in the equations of motion for the particles. In this section, we're going to first write a function that gets the vector between two particles using periodic boundary conditions. Since we're simulating an infinite system, we need to choose the vector $\vec r = \vec r_2 - \vec r_1$ such that its norm is the smallest between all replicas of the box. 

## 4. Equations of motion 

The equation of motion for particle $i$ at time $t$ is given by Newton's second law: 

$$m \ddot{ \vec r}_i = \sum_{j} \vec F_{j \to i}$$

Where the sum runs over all particles $j \neq i$ and $\vec F_{j \to i}$ is the force exerted on particle $i$ by particle $j.$ At every step of our simulation, we need to know what is the net force felt by all particles. Let's write a function which takes in the positions of our particles at a certain time and stores the force felt by each particle. To make things simple, we'll set $m = 1$.

To summarize, we need to implement the following to get the force felt by all the particles: 
* Loop through all particles $i$.
* Loop through all particles $j \neq i$.
* Calculate the force on $i$ due to all the other particles. 

We can be a little bit clever here: since we know that by Newton's third law, $\vec F_{j \to i} = - \vec F_{i \to j},$ we can do the second for loop for $j >i$, making the algorithm slightly more efficient. This step is actually the most expensive part of the simulation since we need to do $\mathcal{O}(N^2)$ operations at each step!

## 5. Numerical integration
Now we need to find a way to perform a discrete integration of the equations of motion, that is, evolve the particles' positions and velocities from time $t \to t + \Delta t$. We'll use a popular molecular dynamics algorithm called the **velocity Verlet algorithm**. You can read more about it [here](https://tonypaxton.org/Notes/MD.pdf).

To see how this works, we can start by Taylor expanding $\vec r_i(t+ \Delta t)$ about $\Delta t$:
$$\vec r_i(t+ \Delta t) = \vec r_i(t) + \vec v_i(t) \Delta t + \frac {\Delta t^2}{2m} \sum_j \vec F_{j \to i} +\mathcal{O}(\Delta t^3) $$

This gives us the positions of the particles at time $t + \Delta t$ with an error of order $\Delta t^3$. To update the velocities, we use: 

$$\vec v_i(t + \Delta t) = \vec v_i(t) + \frac{\Delta t}2 \left( \vec a_i(t + \Delta t) + \vec a_i(t) \right)+ \mathcal{O}(\Delta t^2),$$

where $\vec a_i(t + \Delta t)$ is the acceleration vector of particle $i$ at time $t+ \Delta t$. This gives us the velocities of the particles at the next time step with an error of order $\Delta t^2$. 
In short, we need to do the following at each time step:
* Use the current positions and velocities $\vec r(t),$ $\vec v(t)$ to get the forces acting on the particles at time $t$.
* Implement the equation for $\vec r(t + \Delta t)$.
* Use the new positions to get the forces acting on the particles at time $t+ \Delta t$.
* Implement the equation for $\vec v(t+ \Delta t).$

## 6. Crystallization of particles 

With a very small modification, we can use this code do visualize how a LJ gas turns into a crystal. At low temperatures, the particles will be close to the stable equilibrium of the LJ potential and oscillate about this position. Recall that the temperature is linked to the average kinetic energy of the particles by $\langle K \rangle = k_BT$.
To turn the gas into a crystal we need to gradually lower its temperature by reducing the velocity of the particles. Physically, this corresponds to using a thermostat to regulate the temperature. In this notebook, we use the [**Berendsen thermostat**](https://www2.mpip-mainz.mpg.de/~andrienk/journal_club/thermostats.pdf). The idea is to fix some target temperature $T_0$ and rescale the velocities at each time step such that:

$$\frac{dT}{dt} = \frac{T_0-T}{\tau}$$
where $\tau$ is som characteristic decay time for the temperature that we fix. This amounts to rescaling the velocities of the particles at each time step using: 

$$v(t+\Delta t) = \lambda v(t),$$
where 

$$\lambda = \sqrt{1+ \left( \frac{T_0}{\tau}-1 \right) \frac{\Delta t}{\tau} }.$$
