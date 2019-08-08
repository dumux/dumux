import dune.grid
import dumux.porousmediumflow as pmf

gridView = dune.grid.structuredGrid([0,0], [1,1], [10,10])
print(gridView._typeName)

ttag = pmf.OnePTpfa( gridView )
print(ttag._typeName)

problem = pmf.OnePProblem(gridView, ttag)
