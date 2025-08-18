---
title: Monte Carlo Increased-Radius Floating Random Walk Solution For
  Potential Problems
---

**Abstract**

In this paper, a new Monte Carlo walk method is introduced. The
increased radius floating random walk combines of the two classical
Monte Carlo methods and derived from fixed-radius floating walk method.
In this paper, the method is used to solve typical Laplace's equations
in rectangular region. Also, this method is easily applied to Poisson
equations. Lower walk number and hence lower computation time are
obtained from new method compared with the fixed random walk, floating
random walk and fixed-radius random walk methods. Analyzes were
performed on an average computer and the solution time was reduced by
80%. The results are also compared with Finite Element Method. Increased
radius walk method's results are good agreement with other methods.

***Keywords:*** "Floating random walk, Laplace equations, Monte carlo
method, Numerical methods, Potential problems."

# Introduction

In many engineering disciplines, calculations of potentials or fields
are needed for many problems such as in heat conduction, fluid dynamics
and electromagnetics. The mathematical basis and techniques needed are
well known, but are often difficult to apply in the complicated
three-dimensional geometries. Nowadays finite element method or finite
difference methods are widely used for determining the potential and
field solutions. It is difficult to use finite-element or finite
difference methods to calculate electric potentials or fields near a
high-voltage transmission tower, or stress-control fitting, because of
the large number of mesh points or nodes needed to give an adequate
representation of the geometry. The Monte Carlo methods give a
convenient and flexible means of tackling these and similar problems in
electrical power engineering, and may also be used for potential theory
calculations in other branches of engineering. The goal of the method is
to solve Dirichlet's problem: to find a potential (and its derivatives)
which satisfies Laplace's equation within a given region and takes
specified values on its boundary \[1\]. There are a lot of studies in
the literature about Monte Carlo method \[2-5\]. The Monte Carlo methods
are widely used in electromagnetics problems, areas of applications
include electrostatics, waveguide analysis and antennas. It is easily
applied to Laplace's and Poisson's equations in both two and
three-dimensional cases \[6\], and time-dependent problems \[7\]. There
are many types of Monte Carlo method. Some and significant ones are as
follows:

> 1\. Fixed random walk
>
> 2\. Floating random walk
>
> 3\. Exodus method

The first two methods are the most popular. In fixed random walk, the
step size is fixed and steps of the random walks are constrained to lie
parallel to the coordinate axes. It takes a long time to compute the
potential at a given point. On the other hand, in floating random walk,
the step size is not fixed. Hence the computation is more rapid than in
the fixed random walk. Yu and friends \[8\], in their work, two
techniques are developed to accelerate the floating random walk method
for the electrostatic computation involving rectilinear shapes. Garcia
and Sadiku \[9\] have introduced a new walk procedure in addition to
these existing three types, named fixed-radius floating random walk. The
aim of this method is to decrease the runtime and to eliminate the need
to check for the shortest distance to the border. If the radius of the
circle in floating random walk is fixed, then the programming the random
walk becomes easier. Thus it is more effective and it reduces the
runtime considerably. Although Monte Carlo simulation is a slow and
costly technique, it does have some advantages that are as follows:

> 1\. Don't require input data
>
> 2\. Conceptually easier to understand
>
> 3\. Easy programming

The major disadvantage of the Monte Carlo methods is that they permit
calculating the potential only one point at a time. For whole field
calculations, other numerical techniques such as finite element and
finite difference methods are preferred. For whole field computation,
several techniques have been proposed such as shrinking boundary method
\[10\] and inscribed figure method \[11\].

In this paper, a new walk method, namely increased-radius floating
random walk, is introduced. This method is based on the fixed-radius
floating random walk. In this method, the radiuses of the circles are
increased with Δr at every step. Thus, the arriving time to a border is
dramatically reduced.

# Monte Carlo Method

The mathematical basis of the Monte Carlo method is easily found in the
literature \[12\]. Briefly, two important concepts serve as a
mathematical basis:

1.  The mean value theorem of potential

2.  Green's function of the first kind

Also, generating random numbers and application to the Laplace equation
in rectangular and axisymmetric regions are given in \[13\]. Let us
suppose that the fixed random walk is to be applied to solve Laplace's
equation:

> ![](media/image1.wmf) in region R (1)

subject to the Dirichlet boundary condition,

> ![](media/image2.wmf) on boundary Γ (2)

We now consider cases involving rectangular and axisymmetric solution
regions.

1.  **Rectangular Solution Region**

The region R is divided into a finite difference square mesh and Eq. (1)
is replaced by its finite difference equivalent as \[14\],

> ![](media/image3.wmf) (3)
>
> ![](media/image4.wmf) (4)

To calculate the potential at (*x~n~*,*y~n~*) a particle is asked to
begin a walk at that point. The particle proceeds to wander from node to
node in the grid until it reaches the boundary. It wanders through the
mesh according to the probabilities in Eq. (4) until it reaches the
boundary where it is absorbed and the prescribed potential *V~p~*(1) is
recorded. By sending out *N* particles from (*x~n~*,*y~n~*) and
recording the potential at the end of each walk, we obtain:

> ![](media/image5.wmf) (5)

2.  **Axisymmetric solution region**

For *V=V(ρ, z)*, the finite difference equivalent of Eq. (1) for
![](media/image6.wmf) is:

> ![](media/image7.wmf) (6)

where ![](media/image8.wmf) and the random walk probabilities are given
by:

> ![](media/image9.wmf) (7)

For *ρ* = 0, the finite difference equivalent of Eq. (1) is:

> ![](media/image10.wmf) (8)

so that

> ![](media/image11.wmf) (9)

The random-walking particle is instructed to begin walking at
(*ρ~n~*,*z~n~*). It wanders through the mesh according to the
probabilities in Eq. (7) and Eq. (9) until it reaches the boundary where
it is absorbed and the prescribed potential *V~p~*(1) is recorded. By
sending out *N* particles from (*ρ~n~*,*z~n~*) and recording the
potential at the end of each walk, we obtain the potential at
(*ρ~n~*,*z~n~*) as:

> ![](media/image12.wmf) (10)

# Increased-Radius Floating Random Walk Solution Algorithms

At the beginning, let's describe the solution with fixed-radius random
walk. In order to obtain the potential at a given point P(xn,yn) with
fixed-radius random walk, a circle is placed on point P whose radius is
r0. To perform a step, a random point is selected on this circle. On
this new point, a new circle is placed with the same radius. The steps
are repeated until it reaches to a boundary. Walk is completed when the
boundary is reached.

**3.1. Algorithm for Rectangular Region**

For increased-radius walk algorithm which has been developed in this
paper, same procedure can be followed up. Here is a simplified algorithm
for solving Laplace's equation at a particular point
$P\left( x_{n},y_{n} \right)$:

Step 1. Determine the "radius $r_{0}$" of the walk and the tolerance T:
These values are determined empirically. Record the initial radius r0 as
r for reaching to the border within the tolerance.

Step 2. Perform the step: From the point $P\left( x_{n},y_{n} \right)$,
generate a circle of radius $r_{0}$ and randomly select a point on that
circle. The point is determined by generating a random number in the
interval \[0,1\] and multiplying by 2π to obtain the angle $\theta$. The
following equations are used to determine the new point:

> ![](media/image13.wmf) (11)

![](media/image14.wmf) (12)

Step 3. Increase the radius of circle for next step: Increase the radius
of circle by Δr:

![](media/image15.wmf) (13)

Step 4. Determine if a walk has been performed: From new position, check
to see if the step has reached the border within tolerance T. If T ≥ r,
constraint is satisfied then no step will go outside the border. If T \<
r, suppose that border is reached and record the potential (Vp(i)) at
that particular border and go to Step 2.

Step 5. Compute the potential for point $P\left( x_{n},y_{n} \right)$:
After a sufficient number of walks (N), the potential is determined
using the following formula:

![](media/image16.wmf) (14)

An illustrative description of a typical walk is shown in Fig. 1. The
solution region is a rectangular and the potential values at the
boundaries are V1, V2, V3 and V4 respectively. As shown in Fig. 1, after
four steps, the border has been reached. At fourth step, two possible
points which are can be determined. If walk has reached to the point P1,
then the potential to be recorded is V1, otherwise V2.

![](media/image17.png){width="5.573611111111111in"
height="3.5652777777777778in"}

**Fig. 1. Increased-radius floating random walk for an arbitrary
rectangular region.**

**3.2. Algorithm for Axisymmetric Region**

The algorithm for axisymmetric region is quite similar to the
rectangular case. Because of line of symmetry, two considerations must
be taken into account for axisymmetric regions \[10\].

> 1\. If a step lies on the left side of the line of symmetry, the point
> is reflecting back to the right side.
>
> 2\. If a step lies on the line of symmetry, the next step must have
> angle ($\phi$) as such that $- \pi\  < \ \phi < \ \pi.$

**3.3. Algorithm for Three-Dimensional Rectangular Region**

Similar to the original algorithm, the modified version of the algorithm
for a three-dimensional rectangular solution region is the same as for
the two-dimensional case. Steps in three-dimensional case are performed
on spheres instead of circles. The Radius of spheres on the next step is
increased by $\Delta R$. The next point on a sphere can be determined by
the following equations:

> ![](media/image18.wmf) (15)

In the Eq. (15), $R_{0}$ is the initial radius of sphere. The radius of
sphere on the next step is $R_{0} = R_{0} + \Delta R$.

# Illustrative Examples

**4.1. Comparison with Fixed-Radius Floating Random Walk: Example I**

The new walk method described in this paper is compared with
fixed-radius floating random walk methods. Methods are compared in terms
of speed and accuracy. The methods were implemented on a Core2Duo 3.16
GHz PC. The source codes were written in MATLAB. Example is provided to
illustrate the validity for the increased radius floating random walk.
The exact solution for the problem can be found in \[13\]. As shown in
Fig. 2, the aim of the example is to solve Laplace's equation:

![](media/image19.wmf) (16)

subjected to the following Dirichlet boundary conditions: V(x,0)=0,
V(0,y)=0, V(x,1)=100, V(1,y)=0

The fixed-radius floating random walk and increased-radius floating
random walk methods were carried out with 500, 1000, 1500 and 2000
walks. The initial radius of circle $r_{0}$ for increased radius
floating random walks was chosen as 0.01 and this value is constant for
fixed-radius floating random walk. The increment of radius of circle in
increased-radius floating random walk was chosen as 0.001. Tolerance
value (T) was chosen as 0.005. Every walks in both methods was repeated
10 times. Three typical points were chosen in solution region and the
results are given in Table 1. The results of exact solution for these
three points are 43.20, 25 and 6.79 V, respectively.

![](media/image20.png){width="3.5215277777777776in"
height="3.1909722222222223in"}

**Fig. 2. Solution region for example I.**

Table 1 presents the solutions for three points. In Table 1, shows the
average step size and shows the average runtime values of 10 repeated
walk. Random numbers for each walk are different for both methods. The
and values are depending on random numbers.

**Table 1. Results for example 1 (r0 = 0.01, Δr = 0.001).**

  -------------- --------- ---------------- ------------------------ ------------------------------- -------------------- ------------------------ -------------------------------
  **Point        ***N***   **Fixed-radius                                                            **Increased-radius                            
  (x,y)**                  walk**                                                                    walk**                                        

                           **Pot. (V)**     ![](media/image21.wmf)   ![](media/image22.wmf)**(s)**   **Pot. (V)**         ![](media/image21.wmf)   ![](media/image22.wmf)**(s)**

  (0.25,0.75)    500       42.16            136                      2.72                            42.67                6                        0.12

                 1000      42.92            67                       5.36                            43.28                7                        0.27

                 1500      43.69            50                       8.37                            43.16                8                        0.31

                 2000      43.50            350                      10.80                           43.09                9                        0.42

  (0.5,0.5)      500       24.88            155                      4.50                            24.74                7                        0.13

                 1000      25.22            90                       9.00                            24.86                10                       0.27

                 1500      25.10            148                      13.05                           25.16                9                        0.40

                 2000      25.35            400                      18.00                           25.05                11                       0.53

  (0.75,0.25)    500       6.76             52                       2.82                            6.52                 22                       0.51

                 1000      6.85             210                      5.60                            6.76                 26                       1.03

                 1500      6.80             394                      8.27                            6.75                 42                       1.55

                 2000      6.78             208                      10.90                           6.90                 18                       2.05
  -------------- --------- ---------------- ------------------------ ------------------------------- -------------------- ------------------------ -------------------------------

Fig. 3, 4 and 5 shows the average run time comparison for three points
respectively.

![](media/image23.png){width="4.182638888888889in"
height="2.5215277777777776in"}

**Fig. 3. Average run time comparison for 0.25 mm, 0.75 mm point.**

![](media/image24.png){width="4.190972222222222in"
height="2.5131944444444443in"}

**Fig. 4. Average run time comparison for 0.5 mm, 0.5 mm point.**

![](media/image25.png){width="4.147916666666666in"
height="2.3652777777777776in"}

**Fig. 5. Average run time comparison for 0.75 mm, 0.25 mm point.**

Table 1 and Fig. 3, 4 and 5 show that the average step sizes of
increased-radius floating random walk are less than that of fixed-radius
floating random walk. Hence the runtime are reduced. The runtime is
important when the solution points in the solution region are more than
one point. In this situation, more than one walk must be started
simultaneously. The results are in good agreement with the exact and
fixed-radius floating random walk solutions.

**4.2. Comparison with Finite Element Method: Example II**

In addition, the method proposed here was compared with the finite
element method. The solution region is also in rectangular coordinates
and shown in Fig. 6. The problem is to solve Laplace's equation with the
following Dirichlet boundary conditions:

![](media/image26.wmf)

![](media/image27.png){width="3.252083333333333in"
height="2.9479166666666665in"}

**Fig. 6. Solution region for example II.**

Finite element solution was performed with FEMM software \[15\]. The
mesh created in FEMM model is with 3201 nodes. The results for four
typical points are shown in Table 2. The initial radius for
increased-radius floating random walk was chosen as $r_{0}$ = 0.001 and
the increment was taken as Δr = 0.0001. The tolerance was set to T =
0.005. Each walk is repeated 10 times and the average values of the
steps and time are also shown in Table 2.

As can be seen from Table 2, the results are very close compared with
finite element method. The initial radius of circle and the increment
value can be chosen as small as possible.

Table 3 shows the potential values at three points for various $r_{0}$
values. The increment of radius was chosen as Δr = $0.1r_{0}$. The
number of walks is set to 2000. In Table 3, $\delta$ is the error
estimate, which is obtained by repeating each calculation 10 times
\[16\]. For the increased-radius floating random walk, it was noted that
2000 walks were sufficient for the solution to converge. It seems from
the results that the 0.01 value of initial radius of the circle and the
increment value are also sufficient.

**Table 2. Results for example 2 (r0 = 0.001, Δr = 0.0001).**

  ------------- --------- ---------------- -------------------- ------------------------ ------------------------
  **Point       ***N***   **FEM**          **Increased-radius                            
  (x,y)**                                  walk**                                        

                          **Potential      **Potential (V)**    ![](media/image21.wmf)   ![](media/image22.wmf)
                          (V)**                                                          **(s)**

  (3, 1)        500       8.13             8.00                 45                       1.26

                1000                       8.12                 109                      2.40

                1500                       8.14                 58                       3.70

                2000                       8.15                 74                       4.94

  (3, 3)        500       19.71            19.57                93                       1.33

                1000                       19.60                71                       2.70

                1500                       19.64                80                       4.03

                2000                       19.67                66                       5.40

  (1, 3)        500       11.57            11.53                103                      1.23

                1000                       11.56                77                       2.47

                1500                       11.63                103                      3.70

                2000                       11.52                121                      4.90

  (2.5, 3.5)    500       23.93            23.63                58                       1.16

                1000                       23.80                92                       2.34

                1500                       23.75                88                       3.51

                2000                       23.85                59                       4.66
  ------------- --------- ---------------- -------------------- ------------------------ ------------------------

**Table 3. Potential values and error estimates for example 2 (N = 2000,
Δr = 0.1r0).**

  -------------- ------------- ------------ -------------------------------
  **Point        **FEM (V)**   ***r~0~***   **Increased-radius walk
  (x,y)**                                   (*V±****δ***)**

  (3, 1)         8.13          0.01         8.1510±0.029760

                               0.001        8.1518±0.083900

                               0.0001       8.0840±0.080030

  (3, 3)         19.71         0.01         19.6440±0.05710

                               0.001        19.6740±0.07056

                               0.0001       19.6920±0.05370

  (1, 3)         11.57         0.01         11.6310±0.06174

                               0.001        11.6890±0.09238

                               0.0001       11.5350±0.13117
  -------------- ------------- ------------ -------------------------------

# Conclusions

A new walk method for Monte Carlo method has been implemented in this
paper. The method has been applied to Laplace's solution in rectangular
region in two dimensions. It can be easily applied to other coordinate
systems and three dimensional cases. The results obtained from the new
walk method are in good agreement with those obtained using finite
element solution, fixed-radius floating random walk and exact solution
methods[. Analyzes were performed on an average computer and the
solution time was reduced by 80%.]{.mark}

The advantages of the new method can be summarized as follows:

> 1\. Runtimes are reduced considerably.
>
> 2\. The user can choose the initial radius of circle very small.
>
> 3\. The implementation of the algorithm will be very suitable for
> multi-point calculations performed at the same time.
>
> 4\. It will be very useful when starting one more particle
> simultaneously at a given point.

**References**

\[1\] J. H. Pickles, "Monte carlo field calculations," Proceedings of
the Institution of Electrical Engineers, vol. 124, no. 12, pp. 1271-
1276, 1997.

\[2\] W. Yu, K. Zhai, H. Zhuang, and J. Chen, "Accelerated floating
random walk algorithm for the electrostatic computation with 3-D
rectilinear-shaped conductors," Simulation Modelling Practice and
Theory, vol. 34, pp. 20--36, 2013.

\[3\] K. Chatterjee, "A Parallelized 3D floating random-walk algorithm
for the solution of the nonlinear poisson-boltzmann equation," Progress
in Electromagnetics Research, PIER, vol. 57, pp. 237--252, 2006.

\[4\] R. Sotner, A. Lahiri, A. Kartci, N. Herencsar, J. Jerabek, and K.
Vrba, "Design of novel precise quadrature oscillators employing ECCIIs
with electronic control," Advances in Electrical and Computer
Engineering, vol. 13, no. 2, pp. 65-72, 2013.

\[5\] O. A. Mousavi, M. S. F. Astaneh, and G. B. Gharehpetian,
"Improving power system risk evaluation method using monte carlo
simulation and gaussian mixture method," Advances in Electrical and
Computer Engineering, vol. 9, no. 2, pp. 38-44, 2009.

\[6\] M. N. O. Sadiku, "Monte carlo methods in an introductory
electromagnetic course," IEEE Transactions Education, vol. 33, pp.
47-50, 1990.

\[7\] M. N. O Sadiku, C. M. Akujuobi, and S. M. Musa, "Monte carlo
analysis of time-dependent problems," in Southeastcon' 06 Proceedings of
the IEEE, 2006, pp. 7-10.

\[8\] W. Yu, K. Zhai, H. Zhuang, and J. Chen, "Accelerated floating
random walk algorithm for the electrostatic computation with 3-D
rectilinear-shaped conductors," Simulation Modelling Practice and
Theory, vol. 34, pp. 20--36, 2013.

\[9\] R. C. Garcia, and M. N. O. Sadiku, "Monte carlo fixed-radius
floating random walk solution for potential problems," In Southeastcon
'96. Bringing Together Education, Science and Technology Proceedings of
the IEEE, 1996, pp. 88-91.

\[10\] G. E. Zinsmeister, and J. A. Sawyerr, "A method for improving the
efficiency of monte carlo calculation of heat conduction problems,"
Journal of Heat Transfer, vol. 96, pp. 246--248, 1974.

\[11\] G. E. Zinsmeister, and S. S. Pan, "A modification of the monte
carlo method," International Journal for Numerical Methods Engineering,
vol. 10, pp. 1057--1064, 1976.

\[12\] M. N. O. Sadiku, Monte Carlo Methods for Electromagnetics, CRC
Press, 1st Ed., 2009.

\[13\] M. N. O. Sadiku, Numerical Techniques in Electromagnetics, CRC
Press, 2nd Ed., 2001.

\[14\] M. N. O. Sadiku, and R. C. Garcia, "Whole field computation using
monte carlo method," International Journal of Numerical Modelling," vol.
10, pp. 303--312, 1997.

\[15\] D. Meeker, FEMM (Finite Element Method Magnetics), www.femm.info

\[16\] B. Bowerman, and R. T. O'Connell, Applied Statistics, The
McGraw-Hill Companies, Inc., 1997.
