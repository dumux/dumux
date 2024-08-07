@addtogroup PQ1BubbleDiscretization

An extension of the box method is obtained when enriching the trial space by a bubble function.
While the box method is a control-volume finite element scheme with piecewise constant basis functions,
the bubble scheme uses basis functions enriched by a bubble function. The scheme is introduced in @cite Schneider2024.

![Control volume partitioning for the box with bubble method (PQ1Bubble).](pq1bubble.svg){html: width=50%}

The additional degree of freedom at the element center is used to ensure stability, for example in the case
of velocity spaces for the Stokes problem. An additional overlapping control volume is introduced for
the element degree of freedom.