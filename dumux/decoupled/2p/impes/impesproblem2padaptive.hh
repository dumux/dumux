/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Base class for all 2-phase problems which use an impes algorithm
 * @author Markus Wolff
 */
#ifndef DUMUX_IMPESPROBLEM_2P_ADAPTIVE_HH
#define DUMUX_IMPESPROBLEM_2P_ADAPTIVE_HH

#include <dumux/decoupled/2p/impes/impesproblem2p.hh>


namespace Dumux
{

struct RestrictedValue
{
	double sat;
	double press;
	double volCorr;
	int count;
	RestrictedValue()
	{
		sat = 0;
		press = 0;
		volCorr = 0;
		count = 0;
	}
};
/*!
 * \ingroup IMPESproblem
 * \ingroup IMPES
 * \brief  Base class for all 2-phase problems which use an impes algorithm
 *
 * @tparam TypeTag The Type Tag
 * @tparam Implementation The Problem implementation
 */
template<class TypeTag>
class IMPESProblem2Padaptive : public IMPESProblem2P<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Implementation;
    typedef IMPESProblem2P<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::Grid Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    // material properties
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;


    enum {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<Scalar, dimWorld>      GlobalPosition;

    //typedefs für die Gitteradaption
    //*******************************
    typedef typename Grid::LeafGridView LeafGridView;
    typedef typename Grid::LevelGridView LevelGridView;
    typedef typename LeafGridView::template Codim<0>::Iterator LeafIterator;
    typedef typename LevelGridView::template Codim<0>::Iterator LevelIterator;
    typedef typename Grid::HierarchicIterator SonIterator;
    typedef typename LeafGridView::IntersectionIterator LeafIntersectionIterator;
//    typedef typename LeafIntersectionIterator::Intersection Intersection;
    typedef typename Grid::template Codim<0>::Entity Entity;
    typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;
    typedef typename Grid::LocalIdSet IdSet;
    typedef typename IdSet::IdType IdType;
    //*******************************

    typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;

    IMPESProblem2Padaptive(const IMPESProblem2Padaptive &)
    {}

public:
    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     * \param verbose Output flag for the time manager.
     */
    IMPESProblem2Padaptive(TimeManager &timeManager, Grid &grid)
        : ParentType(timeManager, grid.leafView()), grid_(grid)
    {
    	initialize();
    }
    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     * \param spatialParameters SpatialParameters instantiation
     * \param verbose Output flag for the time manager.
     */
    IMPESProblem2Padaptive(TimeManager &timeManager, Grid &grid, SpatialParameters &spatialParameters)
        : ParentType(timeManager, grid.leafView(), spatialParameters), grid_(grid)
    {
    	initialize();
    }

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     * \param verbose Output flag for the time manager.
     */
    IMPESProblem2Padaptive(Grid &grid, bool verbose = true)
    DUNE_DEPRECATED // use IMPESProblem2P(TimeManager &, const GridView &)
        : ParentType(grid.leafView(), verbose), grid_(grid)
    {
    	initialize();
    }
    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     * \param spatialParameters SpatialParameters instantiation
     * \param verbose Output flag for the time manager.
     */
    IMPESProblem2Padaptive(Grid &grid, SpatialParameters &spatialParameters, bool verbose = true)
    DUNE_DEPRECATED // use IMPESProblem2Padaptive(TimeManager &, const GridView &)
        : ParentType(grid.leafView(), spatialParameters ,verbose), grid_(grid)
    {
    	initialize();
    }

    virtual ~IMPESProblem2Padaptive()
    {
    }

    void initialize()
    {
    	levelmin_ = Params::tree().template get<int>("levelMin");
    	levelmax_ = Params::tree().template get<int>("levelMax");
    	refinetol_ = Params::tree().template get<Scalar>("refineTol");
    	coarsentol_ = Params::tree().template get<Scalar>("coarsenTol");
    }

    void postTimeStep()
    {
    	LeafGridView leafView = this->gridView();
    	ScalarSolutionType indicator;
    	indicator.resize(this->variables().saturation().size());
    	indicator = -1e10;

    	Scalar globalmax = -1e100;
    	Scalar globalmin = 1e100;

    	// 1)
		// Schleife über alle Leaf-Elemente
		for (LeafIterator it = leafView.template begin<0>();
			it!=leafView.template end<0>(); ++it)
		{
			// Find maximal und minimal Saturation
			int indexi = this->variables().elementMapper().map(*it);
			globalmin = std::min(this->variables().saturation()[indexi][0], globalmin);
			globalmax = std::max(this->variables().saturation()[indexi][0], globalmax);

			// Calculate Indicator in all cells
			LeafIntersectionIterator isItEnd = leafView.iend(*it);
			for (LeafIntersectionIterator isIt = leafView.ibegin(*it); isIt!= isItEnd; ++isIt)
			{
				if (!isIt->neighbor())
					continue ;

				const EntityPointer pOutside = isIt->outside();
				const Entity &outside =*pOutside;
				int indexj = this->variables().elementMapper().map( outside );

				if ( it.level() > outside.level() ||
						(it.level() == outside.level() && indexi<indexj) )
				{
					Scalar localdelta = std::abs(this->variables().saturation()[indexi][0] - this->variables().saturation()[indexj][0]);
					indicator[indexi] = std::max(indicator[indexi][0], localdelta);
					indicator[indexj] = std::max(indicator[indexj][0], localdelta);
				}
			}
		}

		Scalar globaldelta = globalmax- globalmin;
   		globaldelta = std::max(globaldelta,0.1);
//		globaldelta = 1;
		int marked = 0;

		for (LeafIterator it = leafView.template begin<0>();
					it!=leafView.template end<0>(); ++it)
		{
			// refine?
			if (indicator[this->variables().elementMapper().map(*it)] > refinetol_*globaldelta
					&& it.level()<levelmax_)
			{
				const Entity &entity =*it;
				grid_.mark( 1, entity );
				++marked;

//				// Taken from grid-howto. What is it for?
//				LeafIntersectionIterator isend = leafView.iend(entity);
//				for(LeafIntersectionIterator is = leafView.ibegin(entity); is != isend; ++is)
//				{
//					const typename LeafIntersectionIterator::Intersection intersection = *is;
//					if(!intersection.neighbor())
//						continue;
//
//					const EntityPointer pOutside = intersection.outside();
//					const Entity &outside =*pOutside;
//					if ((outside.level()<levelmax_))
//						grid_.mark(1, outside);
//				}

			}
			// Coarsen?
			if (indicator[this->variables().elementMapper().map(*it)] < coarsentol_*globaldelta
					&& it.level()>levelmin_)
			{
				grid_.mark( -1, *it );
				++marked;
			}
		}

		if (marked==0)
			return;

		// 3)
		// Sättigung und Druck für alle father-elemente von leaf-elementen berechnen

		std::map<IdType ,RestrictedValue> restrictionmap;
		const IdSet& idset = grid_.localIdSet();

		for (LeafIterator it = leafView.template begin<0>();
					it!=leafView.template end<0>(); ++it)
		{
			// get your map entry
			IdType idi = idset.id(*it);
			RestrictedValue& rv = restrictionmap[idi];

			// put your value in the map
			int indexi = this->variables().elementMapper().map(*it);
			rv.sat = this->variables().saturation()[indexi];
			rv.press = this->variables().pressure()[indexi];
			rv.volCorr = this->variables().volumecorrection(indexi);
			rv.count = 1;

			//Average in father
			if (it.level()>0)
			{
				EntityPointer ep = it->father();
				IdType idf = idset.id(*ep);
				RestrictedValue& rvf = restrictionmap[idf];
				rvf.sat += rv.sat/rv.count;
				rvf.press += rv.press/rv.count;
				rvf.volCorr += rv.volCorr/rv.count;
				rvf.count += 1;
			}
		}

		//4) Adapt und Größe der Vektoren ändern
		//*********************

		// Check for maximum refinement ratio
		for (int level=grid_.maxLevel(); level>=levelmin_; --level)
			{
				LevelGridView levelView = grid_.levelView(level);
				for (LevelIterator it = levelView.template begin<0>();
						it!=levelView.template end<0>(); ++it)
				{
					if (it->isLeaf())
					{
						// compute minimal allowed level for this element
						int maxNeighborLevel=0;
						LeafIntersectionIterator isItEnd = leafView.iend(*it);
						for (LeafIntersectionIterator isIt = leafView.ibegin(*it); isIt!= isItEnd; ++isIt)
						{
							if(!isIt->neighbor())
								continue;
							const EntityPointer pOutside = isIt->outside();
							const Entity &outside =*pOutside;
							int levelOutside=outside.level()+grid_.getMark(*pOutside);
							maxNeighborLevel = (maxNeighborLevel < levelOutside) ? levelOutside : maxNeighborLevel;
						}

						int mark = grid_.getMark(*it);
						// The following "-1" defines the maximum refinement ratio
						int minAllowedLevel = maxNeighborLevel - 1;
						// set refinement mark and remove coarsening mark if neccessary
						if (level+mark < minAllowedLevel)
						{
							mark = minAllowedLevel - level;
							assert((-1 <= mark) and (mark <= 1));
//    		                    std::cout<<"Vorher..."<<mark<<this->variables().index(*it);
							grid_.mark(mark, *it);
//    		                    std::cout<<"...Hinterher"<<std::endl;
						}
					}
				}
				for (LevelIterator it = levelView.template begin<0>();
												it!=levelView.template end<0>(); ++it)
				{
					if((it->isLeaf())
						&&(it->hasFather())
						&&(-1 == grid_.getMark(*it)))
					{
						int maxNeighborLevel=0;
						LeafIntersectionIterator isItEnd = leafView.iend(*it);
						for (LeafIntersectionIterator isIt = leafView.ibegin(*it); isIt!= isItEnd; ++isIt)
						{
							if(!isIt->neighbor())
								continue;
							const EntityPointer pOutside = isIt->outside();
							const Entity &outside =*pOutside;
							int levelOutside=outside.level()+grid_.getMark(*pOutside);
							maxNeighborLevel = (maxNeighborLevel < levelOutside) ? levelOutside : maxNeighborLevel;
						}
						if (maxNeighborLevel>level)
						{
							grid_.mark(0,*it);
						}
					}
				}

				for (LevelIterator it = levelView.template begin<0>();
							it!=levelView.template end<0>(); ++it)
				{
					if((it->isLeaf())
						&&(it->hasFather())
						&&(-1 == grid_.getMark(*it)))
					{
						int maxMark=-1;
						const EntityPointer fatherPointer = it->father();
						const SonIterator sonEndIt(fatherPointer->hend(level+1));
						for (SonIterator sonIt(fatherPointer->hbegin(level+1)); sonIt != sonEndIt; ++sonIt)
						{
							if (sonIt->isLeaf())
							{
								maxMark=std::max(maxMark,
										sonIt->level()-level+grid_.getMark(*sonIt));
							}
						}
//							if (maxMark>=0)
//								grid_.mark(0,*it);
						grid_.mark(std::min(0,maxMark),*it);
					}
				}
			}

		grid_.preAdapt();
		grid_.adapt();
		grid_.postAdapt();

		this->variables().elementMapper().update();

		this->variables().adaptVariableSize2p(this->variables().elementMapper().size());

		//5) Write saturation into new primary variable vector
		for (LeafIterator it = leafView.template begin<0>();
					it!=leafView.template end<0>(); ++it)
		{
			// get your id
			IdType idi = idset.id(*it);

			// check map entry
			typename std::map<IdType , RestrictedValue>::iterator rit
				= restrictionmap.find(idi);
			if (rit!=restrictionmap.end())
			{
				//entry is in map, write in leaf
				int indexi = this->variables().elementMapper().map(*it);
				this->variables().saturation()[indexi] = rit->second.sat/rit->second.count;
				this->variables().pressure()[indexi] = rit->second.press/rit->second.count;
				this->variables().volumecorrection(indexi) = rit->second.volCorr/rit->second.count;
			}
			else
			{
				// value is not in map, interpolate from father element
				EntityPointer ep = it->father();
				IdType idf = idset.id(*ep);
				RestrictedValue& rvf = restrictionmap[idf];
				int indexi = this->variables().elementMapper().map(*it);
				this->variables().saturation()[indexi] = rvf.sat/rvf.count;
				this->variables().pressure()[indexi] = rvf.press/rvf.count;
				this->variables().volumecorrection(indexi) = rvf.volCorr/rvf.count;
			}
		}

//    		this->pressureModel().updateMaterialLaws();
//            this->transportModel().updateMaterialLaws();
		return;
    }


private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc Dumux::IMPESProblem2Padaptive::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

int levelmin_,levelmax_;
Scalar refinetol_,coarsentol_;
Grid& grid_;
};

}

#endif
