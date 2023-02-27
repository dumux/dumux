## gridview

A gridView is a representation of the grid given a hierarchical view and provides read-only access to certain parts of the grid. A LeafGridView is a view on all elements of the grid without descendants in the hierarchy (which would be the leaves of a hierarchy) while a LevelGridView is a view on all elements of a given level of a refinement hierarchy. The grid view is given by the grid. The leaf grid view can be accessed by `gridManager.grid().leafGridView()`.
