from dune.generator.generator import SimpleGenerator

def SpatialParams(gridGeometry, scalar):
    moduleName = 'spatialparams'
    typeName = "Dumux::FVSpatialParamsOnePConstant<{}, {}>".format(gridGeometry._typeName, scalar._typeName)
    includes = ["dumux/material/spatialparams/fv1pconstant.hh"]
    includes += ["dumux/python/material/spatialparams/spatialparams.hh"]
    includes += gridGeometry._includes
    holderType = "std::shared_ptr<{}>".format(typeName)

    generator = SimpleGenerator("SpatialParams", "Dumux::Python")
    module = generator.load(includes, typeName, moduleName, options=[holderType])

    return module.SpatialParams(gridGeometry)
