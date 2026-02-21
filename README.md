# 6 Degree of Freedom Spacecraft Simulation
This project simulates the behavior of a spacecraft in low earth orbit (LEO). This sim, with its basic rigid body dynamics and astrodynamics, may act as a testbed for more advanced GNC topics. I've tried to organize this project in a way that will make code modifications easy. 

## Software dependencies:
Runs on python, needs: numpy, scipy, matplotlib, and plotly. \
Go to this directory and type in the terminal:\
`python3 -m pip install -r requirements.txt`


## How I made this:
### Astrodynamics of Low Earth Orbit
Assuming a perfectly circular orbit at an altitude of 1500km.
Spacecraft position and velocity vectors are expressed in an Earth Centered Inertial Frame
### The Mission (desired body orientation)
The desired spacecraft attitude will be to point directly towards Earth
### Underlying Math (Quaternions, DCMs)
Quaternions are used to express the Directional Cosine Matrix (DCM) from the inertial to spacecraft body frame. Quaternions are a great choice to use because they do not have singularities. 

$$

\boldsymbol{{\beta}} =
\begin{bmatrix}
{\beta}_0 \\
{\beta}_1 \\
{\beta}_2 \\
{\beta}_3
\end{bmatrix}
$$

$$
\mathrm{DCM} =
\begin{bmatrix}
{\beta}_0^2 + {\beta}_1^2 - {\beta}_2^2 - {\beta}_3^2 & 2({\beta}_1 {\beta}_2 + {\beta}_0 {\beta}_3) & 2({\beta}_1 {\beta}_3 - {\beta}_0 {\beta}_2) \\
2({\beta}_1 {\beta}_2 - {\beta}_0 {\beta}_3) & {\beta}_0^2 - {\beta}_1^2 + {\beta}_2^2 - {\beta}_3^2 & 2({\beta}_2 {\beta}_3 + {\beta}_0 {\beta}_1) \\
2({\beta}_1 {\beta}_3 + {\beta}_0 {\beta}_2) & 2({\beta}_2 {\beta}_3 - {\beta}_0 {\beta}_1) & {\beta}_0^2 - {\beta}_1^2 - {\beta}_2^2 + {\beta}_3^2
\end{bmatrix}
$$
### Underlying Kinematics and Dynamics 
#### Kinematic Differential Equation for Quaternions:
$$
\boldsymbol{\dot{{\beta}}} = \frac{1}{2}[B(\boldsymbol{{\beta}})]^{b}\boldsymbol{\omega} 
$$

where 
$$
[B(\boldsymbol{{\beta}})] = 
\begin{bmatrix}
-{\beta}_1 & -{\beta}_2 & -{\beta}_3 \\
{\beta}_0 & -{\beta}_3 & {\beta}_2 \\
{\beta}_3 & {\beta}_0 & -{\beta}_1 \\
-{\beta}_2 & {\beta}_1 & {\beta}_0
\end{bmatrix}
$$

and $\boldsymbol{\omega}$ is the angular velocity of the body frame with respect to the inertial frame, written in body coordinates 
$$
^{b}\boldsymbol{\omega}  = 
\begin{bmatrix}
{\omega}_1 \\
{\omega}_2 \\
{\omega}_3 
\end{bmatrix}

$$

#### Euler's rotational equation of motion:
$$
^{b}\boldsymbol{\dot{\omega}} = (-[I]^{-1})^{b}\boldsymbol{\omega} \times ([I]^{b}\boldsymbol{\omega}) + [I]^{-1}\textbf{M} 
$$

where [I] is the moment of inertia matrix of the spacecraft body fixed frame centered at the center of mass.
Ideally, the body fixed frame will be the principal axis frame so that $[I]$ is a diagonal matrix. \
\
$\textbf{M}$ is a vector representing the applied moment. The applied moment can come from the environment or onboard controllers. 




### Spacecraft
For simplicity, let the spacecraft be a solid box with dimensions 20cm x 20cm x 30cm (a.k.a a 12U cubesat). Assuming that each 10cm cube in a cubesat is 1kg, $\mathrm{mass} = 12kg$

![alt text](<Screenshot 2026-02-20 at 6.27.38 PM.png>)

The above image describes the body fixed frame, centered at the center of mass, which is the geometric center for simplicity's sake. 

The X,Y, and Z axes are chosen to be normal to each face of the cubesat, as that is the principal axis frame for a cuboid. 

$[I]$ for a cuboid is given as 
$$
[I] = 
\mathrm{mass} \begin{bmatrix}
\frac{1}{12}(b^2 + c^2) & 0 & 0 \\
0 & \frac{1}{12}(a^2 + c^2) & 0 \\
0 & 0 & \frac{1}{12}(a^2 + b^2)
\end{bmatrix}
$$

For our cubesat, 
$$
[I]_{cubesat} = 
 \begin{bmatrix}
.13 & 0 & 0 \\
0 & .13 & 0 \\
0 & 0 & .08
\end{bmatrix} \mathrm{kg*m^{2}}
$$
### Attitude Determination
### Sensors
### Actuators
### Future additions


