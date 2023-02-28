## gridmanager

A grid manager handles the grid data from an input (file) and constructs a grid from that information. Only works for supported grid types: YASP, OneD, Alberta, UG, ALU, FOAM, SP, MMESH, Sub. A grid can be constructed either from parameters or from files.

### Functionalities

1. init()
2. grid()
3. loadBalance()
4. makeGridFromFile()
5. makeGridFromDgfFile()
6. makeStructuredGrid()

### Overview
```{mermaid}
flowchart LR
    A(Grid) -->|"template parameter for declaration"| B(gridManager)
    B -->|"init()"| C(grid)
    click A "./grid.html"
    click B "./gridmanager.html"
    click C "./grid.html"
```
