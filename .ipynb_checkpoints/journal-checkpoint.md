# Weekly progress journal

## Instructions

In this journal you will document your progress of the project, making use of the weekly milestones.

Every week you should 

1. write down **on Wednesday** a short plan (bullet list is sufficient) of how you want to 
   reach the weekly milestones. Think about how to distribute work in the group, 
   what pieces of code functionality need to be implemented. 
2. write about your progress **before** the Tuesday in the next week with respect to the milestones.
   Substantiate your progress with links to code, pictures or test results. Reflect on the
   relation to your original plan.

We will give feedback on your progress on Tuesday before the following lecture. Consult the 
[grading scheme](https://computationalphysics.quantumtinkerer.tudelft.nl/proj1-moldyn-grading/) 
for details how the journal enters your grade.

Note that the file format of the journal is *markdown*. This is a flexible and easy method of 
converting text to HTML. 
Documentation of the syntax of markdown can be found 
[here](https://docs.gitlab.com/ee/user/markdown.html#gfm-extends-standard-markdown). 
You will find how to include [links](https://docs.gitlab.com/ee/user/markdown.html#links) and 
[images](https://docs.gitlab.com/ee/user/markdown.html#images) particularly
useful.

## Week 1
(due before 18 February)

- We will store each particles velocity and postion in a 1xN dimensional array. For the initial values we will use a gaussian distribution. This will be done by Kadhim

- first differentiate the lennard jones potential analytically, then we use that to calculate the force. For this we need the distance for the nearest neighbor. We will examine the elements of the position vector and subtract each element to get the shortest distance. This will be done by Achmed

-  for this task a for loop will be used to calculate the time evolution for both position and velocity. This will be done by Kadhim

- We set up a box of lenght L. Once a particle crosses the boundary we will subtract the length of the box to that particle. This will be done by Achmed

- The energy function will consist of the kinetic energy and potential(lennard jones potential). This will be done by both of us

#Progress:

- At first, the packages numpy and matplotlib.pyplot were imported. The next step is to define all functions beforehand to be able to recall them later.

- input parameters such as box lentgh, number of particles, time steps and constants (SI units) are then defined.

- Position and velocity arrays were defined for a 2D problem. For that np.random is used which is based on a uniform distribution.

- In order to keep track of the system, 2 atoms were used instead of 100. 

- The analytical expression for the potential energy (Lennard Jones) is used to determine the force. This is done by a matrix calculation instead of using a for loop.

- The position matrix is filled with 1's on the diagonal to avoid division by zero. This is set again to zero after the division.

- The modules is used to implement the periodic boundry conditions. If x is longer then the side lentgh of the box, the the "rest" (modules) will be taken as numerical value. This again to avoid using for loop.

- For the evolution in time The Euler method is used. A for loop is used for the arrays to evolve in time. 

- At the end the function of the total energy is recalled and all numerical data can be printed. 

- low values for the energies were to expect since SI units are used. 




## Week 2
(due before 25 February)

- First we are going to generalize the code by defining matrices instead of column vectors for each direction.

- An analytical epxression will be derives for the kinetic energy.

- A new expression will be derived for the position, momentum and the time in order to change the molucialr dynamics representation of the problem.

- Our is to use the modules method the impliment the minimal image convention.

- A new dimension will be added to the code in order to be a 3D problem. Depending op bullet point 3, we will discuss a suitable size for the box. We may also change the distrubtion of the particles in order to be close to the boundary of the box.

- The new position coordinates that evolves in time will be stored in a array to keep track of the system.

## Week 3
(due before 3 March)


## Week 4
(due before 10 March)


## Week 5
(due before 17 March)
