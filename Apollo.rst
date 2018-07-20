.. Apollo project master documentation follow
   This file must at least contain the root `toctree` directive.

Apollo project
~~~~~~~~~~~~~~

The Apollo project's goal is to develop open-source software for rapidly
prototyping and evaluation methods relevant to the mission of the United States
Air Force.  The United States Air Force Office of Scientific Research (AFOSR)
has invested heavily in the development of higher order discretization methods,
scalable solver approaches, adaptive mesh refinement (AMR) and Aribtrary
Langrangian Eulerian (ALE) methods.  However, it is not always clear how well
these methods can be applied to the production level problems required for
use in the engineering design.  The Apollo project aims to fulfill that need.

This large task is ameliorated by leveraging previous open-source efforts,
especially the Camellia toolkit for discontinuous Petrov-Galerkin (DPG) methods,
and the Trilinos mathematics library, and the MOAB library for meshing.  By
leveraging these libraries, the Apollo project can focus on just the specifics
of defining the appropriate test problems and applying a range of methods to
these problems.

The disadvantage of leveraging these libraries is that using and building these
libraries can be onerous.  To solve this problem, we use a `git` repository
named *apolloall* to act as a container for our software libraries, and the new
*spack package manager* to facilitate this development.   

We discuss using *apolloall* next.  To go to the evaluation problems, one can
jump directly to the *theaceae* documentation.


.. toctree::
   :maxdepth: 3

   theadocs/apollo_readme.rst
   Theaceae.rst



Indices and tables
~~~~~~~~~~~~~~~~~~

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
