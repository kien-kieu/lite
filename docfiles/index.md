Introduction {#index}
============

<!-- Line Tessellation (LiTe) library
     |||Development version
     Authors: Katarzyna Adamczyk and Kiên Kiêu.
     |||Copyright INRA 2006-yyyy.
     Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
     License: GPL v3. -->

\authors Katarzyna Adamczyk-Chauvat and Kiên Kiêu, UR1404 MaIAGE, INRA, Jouy-en-Josas, France
\copyright INRA
\version development

LiTe is a C++ library devoted to planar Gibbsian T-tessellations. It is a software companion of the paper "A completely random T-tessellation model and Gibbsian extensions" published in Spatial Statistics 2013, vol. 6 (preliminary version freely available on [HAL](http://hal.archives-ouvertes.fr/hal-00785980) and [arXiv](http://fr.arxiv.org/abs/1302.1809)). Note that LiTe is also available to R users, see the [RLiTe](rlite.html) page.

Using LiTe, one can 
- Perform some basic operations on T-tessellations.
- Define a Gibbs model of T-tessellations.
- Simulate a given Gibbs model of T-tessellations.
- Estimate the parameters of a given model. 

What are T-tessellations? All their internal vertices are T-vertices. A vertex is a T-vertex if it is of degree 3 and 2 of its incident edges are aligned. An  example of T-tessellation is shown below.
\image html static_ttessellation.png "T-tessellation in a square domain"
\image latex static_ttessellation.pdf "T-tessellation in a square domain"

LiTe stands for _Line Tessellation_. By line tessellation, we mean tessellations that are built upon a line pattern rather than a point pattern. Examples of point tessellations are Voronoi and Laguerre tessellations. A T-tessellation is only a special case of line tessellation. Other types of line tessellations may be supported in LiTe in the future.

In addition to the reference manual documenting classes and their members, specific documentation topics are provided: check the [Related Pages tab](pages.html).

LiTe relies heavily on [CGAL](http://www.cgal.org), the Computational Geometry Algorithms Library. By using CGAL data structures and algorithms, LiTe developers saved a lot of time and effort. The slight drawback of using CGAL is that LiTe depends on a third-party library that must be installed separately by LiTe users. Fortunately, CGAL has many users working in diverse computing environments (hardware, OS, compiler). Therefore CGAL installation is largely tested and documented.

Development status
==================
The development of LiTe started many years ago. So many years that we are not able to find the exact date! It was at least before 2006. At that time, the core classes and functions of LiTe were already written. And they do not have changed too much since. In that sense, LiTe seems rather stable. But this may be due to the fact that LiTe developers were also almost the only LiTe users (although LiTe has been used occasionaly by close collaborators in the past years).

Due to its rather small user community, LiTe has not been extensively tested. And some unexpected bugs are likely.

Also it should be noticed that we are not professional developers. Therefore we are prone to coding misbehaviour that can make LiTe not a very sustainable library. An example of that is the scheme used for naming class method. At the beginning of LiTe development, we used lower-case names with underscores as word separators. Recently, we turned to a scheme without separators but with upper-case letters for starting each word. We may have to unify the naming scheme in the future.

No major changes are planned in the near future. In the medium term, LiTe and [EBSpatCGAL](https://github.com/rcqls/EBSpatCGAL) may merge. 

Licence
=======
LiTe is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LiTe is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LiTe.  If not, see <http://www.gnu.org/licenses/>.
