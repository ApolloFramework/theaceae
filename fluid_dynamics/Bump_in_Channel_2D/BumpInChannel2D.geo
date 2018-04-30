/*********************************************************************
 *
 *  ehowell 4/30/2018
 *
 *  This is an attempt to create a gmsh mesh consistent with the 
 *  2D Bump-in-channel verification case on NASA's turbulence
 *  modeling webpage. A description of the test problem can be found
 *  at the url: https://turbmodels.larc.nasa.gov/bump.html
 *
 *  This problem has a plate located at 0.0 < x < 1.5. The plate has
 *  a bump with the following parameterization:
 *  y = 0 for x in [0,0.3) + (1.2,1.5] 
 *  y = 0.05*(sin(pi* x/0.9 - pi/3))**4 for x in [0.3,1.2]
 *
 *  The computational domain extended beyond the plate
 *  the size of the domain is not specified in the problem
 *  but the rough dimensions shown in the diagram are [-25,25]x[0,5]
 *  some experimentation may be needed to find the right domain 
 *********************************************************************/

// Define the dimensions of the domain

x0 = -25; //left side of the domain
xf = 25;  //right side of the domain

y0 = 0.0 ; //bottom of domain
yf = 5.0 ; // top of the domain

z = 0.0 ; // this is a 2d problem

xp0 = 0.0; // start of the plate
xp1 = 0.3; // stat of the bump

xp2 = 1.2; // end of the bump
xp3 = 1.5; // end of the plate

npts = 100; //number of points to represent the bump (includes xp1 and xp2)

// todo
// tune mesh scales when I get a working eddy-viscosity model
lc1 = 5e-1 ; // characteristic mesh scale length far from the plate
lc2 = 1e-2 ; // characteristic mesh scale near bump


//These points define the boundary of the domain
Point(1) = {x0, y0, z, lc1} ;
Point(2) = {xf, y0, z, lc1} ;
Point(3) = {xf, yf, z, lc1} ;
Point(4) = {x0, yf, z, lc1} ;

// These points define the start of the plate, end of the plate
// start of the bump and end of the bump
Point(5) = {xp0, y0, z, lc2} ;
Point(6) = {xp1, y0, z, lc2} ;
Point(7) = {xp2, y0, z, lc2} ;
Point(8) = {xp3, y0, z, lc2} ;

//Define the points the lie on the bump
For i In {1:(npts)}
    xi = xp1+i*(xp2-xp1)/(1.0*npts+1.0) ;
    yi = 0.05*(Sin(Pi* xi/0.9 - Pi/3))^4 ;
    Point (8+i) = {xi, yi, z, lc2} ;
EndFor 


Line(1) = {1,5} ;
Line(2) = {5,6} ;
Line(3) = {7,8} ;
Line(4) = {8,2} ;
Line(5) = {2,3} ;
Line(6) = {3,4} ;
Line(7) = {4,1} ;


Spline(8) = {6,(8+1) :(8+npts),7} ;

// Define a line loop that encircles the domain
Line Loop(1) = {1,2,8,3,4,5,6,7} ;

// Define the plane that represents the 2D computational domain
Plane Surface(1) = {1} ;