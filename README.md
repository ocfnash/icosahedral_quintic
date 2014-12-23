The icosahedral solution of the quintic in Python
=================================================

Although the general quintic equation cannot be solved 'in
radicals', it can be solved using some additional special
functions associated to the icosahedron.

The scripts here implement the icosahedral solution
for quintics that take the form:
 * y^5 + y + c = 0 in `bring-jerrard_solution_demo.py`
 * y^5 + 5ay^2 + 5by + c = 0 in `quintic_icosahedral_solution_demo.py`

Note that it is possible to reduce a general quintic equation to 
either of these forms and it requires solving only one auxiliary
quadratic equation to reduce it to the latter form.

I waffled on about this in a bit more detail in the notes
published [here](http://www.sciencedirect.com/science/article/pii/S0723086913000571)
available for free at [this link](http://arxiv.org/abs/1308.0955).
