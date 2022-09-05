#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/defaultvariables.hh>
#include <dumux/experimental/new_assembly/dumux/timestepping/multistagevariables.hh>
#include <dumux/experimental/new_assembly/dumux/timestepping/concepts.hh>

int main(int argc, char* argv[])
{
    using Dofs = double;
    using Variables = Dumux::DefaultVariables<Dofs>;
    using MSVariables = Dumux::MultiStageVariables<Variables>;
    using MSVariablesView = Dumux::MultiStageVariablesView<Variables>;

    static_assert(Dumux::Concepts::MultiStageVariables<MSVariables>);
    static_assert(Dumux::Concepts::MultiStageVariables<MSVariablesView>);

    static_assert(!Dumux::Concepts::View<MSVariables>);
    static_assert(Dumux::Concepts::View<MSVariablesView>);

    {  // constructible from temporary
        MSVariables msVars{Variables{Dofs{42.0}}};
        if (msVars.dofs() != 42.0)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected dof value");
    }

    {  // constructible from ctor arguments
        MSVariables msVars{Dofs{42.0}};
        if (msVars.dofs() != 42.0)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected dof value");
    }

    {  // test a view on variables
        Variables vars{10.0};
        MSVariablesView msVarsView{vars};
        if (msVarsView.dofs() != 10.0)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected dof value");

        vars.update(Dofs{42.0});
        if (msVarsView.dofs() != 42.0)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected dof value");
    }

    return 0;
}
