#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <stdio.h>

template<class G, class V>
void vtkout (const G& grid, const V& Sat,const V& press,const V& perm, const char* name, int k)
{
    Dune::VTKWriter<G> vtkwriter(grid);
    char fname[128];
    sprintf(fname,"%s-%05d",name,k);
    vtkwriter.addCellData(Sat,"saturation Sw");
    vtkwriter.addCellData(press,"total pressure p~");
    vtkwriter.addCellData(perm,"absolute permeability k");
    vtkwriter.write(fname,Dune::VTKOptions::ascii);
}

template<class G, class V>
void vtkout_compressible (const G& grid, const V& press, const V& velocX, const V& velocY, const V& perm, const char* name, int k)
{
    Dune::VTKWriter<G> vtkwriter(grid);
    char fname[128];
    sprintf(fname,"%s-%05d",name,k);
    vtkwriter.addCellData(press,"pressure p");
    vtkwriter.addCellData(velocX,"vellocityX");
    vtkwriter.addCellData(velocY,"vellocityY");
    vtkwriter.addCellData(perm,"absolute permeability k");
    vtkwriter.write(fname,Dune::VTKOptions::ascii);
}

template<class G, class V>
void vtkout_pipeflow (const G& grid, const V& press, const char* name, int k)
{
    Dune::VTKWriter<typename G::LeafGridView> vtkwriter(grid.leafView());
    char fname[128];
    sprintf(fname,"%s-%05d",name,k);
    vtkwriter.addVertexData(press,"pressure");
    //  vtkwriter.vertexData(velocity,"vellocity");
    vtkwriter.write(fname,Dune::VTKOptions::ascii);
}

template<class G, class V>
void vtkout2P2C (const G& grid, const V& Sat, const V& totC, const V& Cw, const V& Cn,  const V& press,const V& perm, const char* name, int k)
{
    Dune::VTKWriter<G> vtkwriter(grid);
    char fname[128];
    sprintf(fname,"%s-%05d",name,k);
    vtkwriter.addCellData(Sat,"saturation Sw");
    vtkwriter.addCellData(totC,"total concetration C1");
    vtkwriter.addCellData(Cw,"wetting phase concentration C1w");
    vtkwriter.addCellData(Cn,"nonwetting phase concentration C1n");
    vtkwriter.addCellData(press,"total pressure p~");
    vtkwriter.addCellData(perm,"absolute permeability k");
    vtkwriter.write(fname,Dune::VTKOptions::ascii);
}
