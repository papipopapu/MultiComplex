# multi-complex

Multi complex numbers allow for fast, machine limit precision differenciation ([e.g.](https://dl.acm.org/doi/abs/10.1145/3378538)), while also being
on their own a pretty interesting subject.

This is a header only implementation of multicomplex numbers through recursive templating. Dirty template metaprogramming makes it faster than any other implementation I have seen if well optmised. Also allows inverse functions like the logarithm and inverse trigonometric functions, although the brach cut chosen makes them 
only bijective for the first complex level (C0).

Some stuff could be improved but overall pretty decent.
