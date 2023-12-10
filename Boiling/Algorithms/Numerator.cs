using Boiling.FemContext.BasisInfo;
using Boiling.FemContext.BasisInfo.Interfaces;
using Boiling.Meshing;

namespace Boiling.Algorithms;

public static class Numerator
{
    public static BasisInfoCollection NumerateBasisFunctions(Mesh mesh, IBasis basis)
    {
        var basisInfo = new BasisInfoCollection(mesh.Elements.Length, basis.BasisSize);

        for (var ielem = 0; ielem < mesh.Elements.Length; ielem++)
        {
            var nodes = mesh.Elements[ielem].Nodes;
            basisInfo[ielem, 0] = new BasisInfoItem(nodes[0], BasisFunctionType.ByGeometricNode, nodes[0]);
            basisInfo[ielem, 1] = new BasisInfoItem(nodes[1], BasisFunctionType.ByGeometricNode, nodes[1]);
            basisInfo[ielem, 2] = new BasisInfoItem(nodes[2], BasisFunctionType.ByGeometricNode, nodes[2]);
            basisInfo[ielem, 3] = new BasisInfoItem(nodes[3], BasisFunctionType.ByGeometricNode, nodes[3]);
        }

        return basisInfo;
    }
}