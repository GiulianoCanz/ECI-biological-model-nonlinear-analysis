# ECI biological nonlinear model

 The code in this repo is inspired by the matlab analysis performed in the paper by Mingyang Lu, Bin Huang, Samir M. Hanash and Eshel Ben-Jacob, available here https://www.pnas.org/doi/full/10.1073/pnas.1416745111. The experimental parameters used in the analysis are shown in the appendix of the paper. 

 The code was tested with Matlab R2021a.
 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ## Introduction

  Matlab analysis of the biological non-linear model known as ECI model (Exosome exchange-based Cancer–Immunity interplay model), which takes into account the role of the exosomes in the interplay between cancer cells (C) and immune cells (in particular two kind of immune cells: Dendritic (DC or D) and Killer (K) ).
  
  The state variables are the densities of DC, C and K, whose time dynamics are specified by the three deterministic non-linear differential equations that describe the generatrix functions of the state flow:
  ![2022-11-30 (2)](https://user-images.githubusercontent.com/48331066/204782450-202ef4ed-e0dd-4483-971d-a1a0c6187e91.png)

  Where 'g' is the basal cell proliferation ratio, 'k' is the basal cell apoptosis ratio, ' H(*,*) ' is the Shifted Hill function that represents the intercellular communication (it’s a parameterized version of the Hill function). In this function a coefficient 'lambda' is involved, it’s the fold change:
 – If lambda > 1 → activated communication;
 – If lambda < 1 → inhibited communication;
 – If lambda = 1 → no effect, in this case H=1.
 It can be said that λ<sup>+</sup> is the term referred to the activation, instead λ<sup>−</sup> is referred to the inhibition. It’s important to specify that the positive terms of the equations describe the proliferation, instead the negative terms are the apoptosis contributions although the presence of an activating Shifted Hill functions (λ+).
 Furthermore, also the effect caused by the immune recognition 'ρ' is included in the model, since it modulates the basal values of λ<sup>+</sup><sub>DC</sub> and λ<sup>+</sup><sub>KC</sub>. during the tumour growth, ρ gradually increases from 0 to 1, then it can eventually decrease. Since the dynamics of ρ is much slower comparing to the interplay between C, K and D cells, it can be neglected, and set static to a constant value. So a 'quasi-static approximation' is made along all the analysis.

 It's possible to briefly describe the ECI model in this way:

 • C are self-activating among themselves, through the exchange of exosomes (reprogramming of neighboring cells);

 • DC are self-activating;

 • DC activate K's proliferation (chemokine-based mechanism);

 • K induce the maturation of other DC;

 • K induce the apoptosis of C;

 • C inhibit K, slowing their proliferation and increasing their
   apoptosis;

 • DC kill C directly with exosomes (exosomes with tumor sup-
   pressor factors);

 • C induce DC maturation (antigen mechanism). This happens if
   the concentration of C is below a certain threshold, otherwise C inhibit the differentiation of DC through exosomes.

 However two cases are analysed in this work:

 - In case I the exosomes action is considered in the interaction between C and DC, and also between DC and C;
 - In case II the exosomes action included in the interactions of the case I is not considered.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## The code in ECI_model.m
ECI_model.m is the main script that the user has to run. Here the workflow followed in the script:
1. Line 23 : firstly the proper value of the immune recognition 'rho' has to be set. In the paper 3 values of rho have been chosen: 0.2 (low) , 1 (normal), 2 (enhanced) ;
2. Line 29 : Then the user has to choose the case to analyze. It is needed to uncomment just the case of interest, according to the considerations written in the       introduction of this readme file ;
3. Line 51 : choose the time interval to provide as input to the model ;
4. Line 52 : set the initial conditions, i.e. the initial densities (thousands of cells / microliter) of the cells (respectively C, DC and K) ;
5. Line 53 : the ode45 method is used to solve the system of differential equations ;
6. Line 56 : the three dimensional phase space is plotted taking into account DC density as x axis, C density as y axis and K density as z axis;
7. Line 60 : time trajectories are drawn. 
8. Line 84 : The phase plane and the vector field of C cells vs DC cells is drawn using the quiver() function. The plot is obtained by computing the density value of K cells for which dK/dt = 0 for each pair of C and DC ;
9. Line 99 : the two nullclines are computed and plotted on the phase plane, by using the contour() function. The first nullcline is the curve where dD = 0 & dK = 0, and the second is the one for dC = 0 & dK = 0  ;
10. Line 121 : the (pseudo) fixed points are computed and drawn on the phase plane, by exploiting the functions getContourLineCoordinates() and intersections() ;
11. Line 133 : the stability of the fixed points is analyzed. In the 'eigvi' matrix we find the eigenvalues of the jacobian matrix linearized around the equilibrium points ;
12. Line 146 : the bifurcation diagrams are computed for C, DC, an K cells, with respect to a variation of parameter rho (the immune recognition). In this diagrams, blue dots represent stable fixed points (both stable nodes and stable foci), and red ones represents unstable fixed points (both unstable nodes and foci) ;

Some output images related to rho = 2 and case 1 :
![3d_phasespace_rho2](https://user-images.githubusercontent.com/48331066/204828044-1e52869f-7cf5-4404-b326-2a9cad720a4a.png)
![Kcells_rho2](https://user-images.githubusercontent.com/48331066/204828119-ec373e68-9df3-43d0-bd3f-77f8e33b1058.png)
![phase_plane_rho2](https://user-images.githubusercontent.com/48331066/204828194-91f560ed-1cb1-4bca-b48b-82a55dedcfae.png)
![bifurcations_rho2](https://user-images.githubusercontent.com/48331066/204828221-d0d5f8e9-c078-40d5-96ca-62adf7b4bdea.png)



-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## The files
- ECI_model.m = main script for the model's analysis;
- ECI_equations.m = code of model's equations (the generatrix functions of the state flow);
- hill_function.m = code of the shifted hill function equation;
- p_hill_function.m = code of the pointwise version of the shifted hill function equation;
- parameters.mat = model parameters;
- intersections.m = code of the function adopted to find out the quasi steady state equilibrium points;
- stability.m = code of the function used to verify the points stability in the bifurcation diagrams.
