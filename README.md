# multi-complex

Multi complex numbers allow for fast, machine limit precision differenciation ([e.g.](https://dl.acm.org/doi/abs/10.1145/3378538)), while also being
on their own a pretty interesting subject.

This is a fast, concise, no external dependencies implementation of multicomplex numbers through recursive templating. Compiles slow, runs like
lightning. Also allows inverse functions like the logarithm and inverse trigonometric functions, although the brach cut chosen makes them 
only bijective for the first complex level (C0).

It took some brain squeezing to make it actually work, and I have seen papers failing to implement inverse functions or make programs fast enough
to be used in differentiation applicatoins, so please reference the existance of this repo if used :). 
