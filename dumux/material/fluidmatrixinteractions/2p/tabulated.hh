#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include <dumux/material/fluidmatrixinteractions/2p/materiallaw.hh>

namespace Dumux::FluidMatrix {

template<int steps>
class TabulatedProperties
{

public:
    template<class Scalar>
    struct Params
    {
        Params (std::string basename)
        {
            init(basename);
        }

        void init(std::string basename)
        {
            basename_ = basename;
            init();
        }

        void init()
        {
            tables_[0] = readTable_(basename_ + "-" + relationFileNames_[0] + ".dat",0);
            for (int i = 0; i < relationFileNames_.size(); ++i)
            {
                auto relation = relationFileNames_[i];
                tables_[i+1] = readTable_(basename_ + "-" + relation + ".dat",1);
            }
        }

        std::array<Scalar, steps> table(int property) const
        {
            return tables_[property];
        }

        Scalar pcEntry() const
        {
            return tables_[1][steps-1]/factor_;
        }

        bool operator== (const Params& p) const
        {
            return p.basename_ == basename_;
        }

    private:
        std::array<Scalar, steps> readTable_(std::string filename, int column)
        {
            std::ifstream tableFile(filename);
            std::string line;
            std::array<Scalar, steps> table;
            int i = 0;
            while (std::getline(tableFile, line) && i < steps)
            {
                size_t pos = 0;
                size_t col = 0;
                Scalar token;
                while (col < column && (pos = line.find(" ")) != std::string::npos)
                    line.erase(0, pos+1);
                token = stod(line.substr(0, line.find(" ")));
                // TODO: flexible direction
                table[steps-i-1] = token;
                ++i;
            }
            return table;
        }

        std::string basename_;
        std::array<std::string, 3> relationFileNames_ = {"CapillaryPressureVsAverageSaturation",
            "Wetting-EffPermeabilityVsAverageSaturation", "Non-wetting-EffPermeabilityVsAverageSaturation"};
        // sw, pc, krw, krn
        std::array<std::array<Scalar, steps>, 4> tables_;
    };

    template<class Scalar = double>
    static Params<Scalar> makeParams(const std::string& paramGroup)
    {
        std::string basename = getParamFromGroup<std::string>(paramGroup, "TabulatedPropertiesName");
        return Params<Scalar>(basename);
    }

    template<class Scalar>
    static Scalar endPointPc(const Params<Scalar>& params)
    { return params.pcEntry(); }

    template<class Scalar>
    static Scalar pc (Scalar sw, const Params<Scalar>& basicParams)
    {
        return interpolate_ (0, 1, sw, basicParams);
    }

    template<class Scalar>
    static Scalar swe (Scalar pc, const Params<Scalar>& basicParams)
    {
        return interpolate_ (1, 0, pc, basicParams, -1);
    }

    template<class Scalar>
    static Scalar krw (Scalar sw, const Params<Scalar>& basicParams)
    {
        return interpolate_ (0, 2, sw, basicParams);
    }

    template<class Scalar>
    static Scalar krn (Scalar sw, const Params<Scalar>& basicParams)
    {
        return interpolate_ (0, 3, sw, basicParams);
    }

private:
    template<class Scalar>
    static int index_ (int property, Scalar value, const Params<Scalar>& params, int order = 1)
    {
        const auto& table = params.table(property);
        const auto comp = [order](Scalar a, Scalar b){return order*a < order*b;};
        auto it = std::lower_bound(table.begin(), table.end(), value, comp);
        return std::clamp(int(std::distance(table.begin(), it)), 1, steps-1)-1;
    }

    template<class Scalar>
    static Scalar interpolate_ (int in, int out, Scalar value, const Params<Scalar>& params, int order = 1)
    {
        auto index = index_(in, value, params, order);
        const auto& tablei = params.table(in);
        const auto& tableo = params.table(out);
        Scalar offset = std::clamp((value-tablei[index])/(tablei[index+1]-tablei[index]),0.0,1.0);
        return offset * (tableo[index+1]-tableo[index]) + tableo[index];
    }
};

template<int steps, class Scalar = double>
using TabulatedPropertiesDefault = TwoPMaterialLaw<Scalar, TabulatedProperties<steps>,
      NoRegularization, TwoPEffToAbsDefaultPolicy>;

} // end namespace Dumux::FluidMatrix
