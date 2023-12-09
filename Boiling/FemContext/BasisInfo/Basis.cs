using Boiling.FemContext.BasisInfo.Interfaces;
using Boiling.Meshing.Geometry;

namespace Boiling.FemContext.BasisInfo;

public readonly record struct LinearBasis : IBasis
{
    public int BasisSize => 4;

    public double Phi(int function, double x, double y)
        => function switch
        {
            0 => (1.0 - x) * (1.0 - y),
            1 => x * (1.0 - y),
            2 => (1.0 - x) * y,
            3 => x * y,
            _ => throw new ArgumentOutOfRangeException(nameof(function), function, "Not expected function number!")
        };

    public double DPhi(int function, int variable, double x, double y)
        => variable switch
        {
            0 => function switch
            {
                0 => y - 1.0,
                1 => 1.0 - y,
                2 => -y,
                3 => y,
                _ => throw new ArgumentOutOfRangeException(nameof(function), function, "Not expected function number!")
            },
            1 => function switch
            {
                0 => x - 1.0,
                1 => -x,
                2 => 1.0 - x,
                3 => x,
                _ => throw new ArgumentOutOfRangeException(nameof(function), function, "Not expected function number!")
            },
            _ => throw new ArgumentOutOfRangeException(nameof(variable), variable, "Not expected var number!")
        };
}

// public struct BiQuadraticBasis : IBasis
// {
//     public int BasisSize => 9;
//     
//     public double Phi(int function, double x, double y)
//     {
//         return function switch
//         {
//             0 => 4.0 * (x - 0.5) * (x - 1.0) * (y - 0.5) * (y - 1.0),
//             1 => -8.0 * x * (x - 1.0) * (y - 0.5) * (y - 1.0),
//             2 => 4.0 * x * (x - 0.5) * (y - 0.5) * (y - 1.0),
//             3 => -8.0 * (x - 0.5) * (x - 1.0) * y * (y - 1.0),
//             4 => 16.0 * x * (x - 1.0) * y * (y - 1.0),
//             5 => -8.0 * x * (x - 0.5) * y * (y - 1.0),
//             6 => 4.0 * (x - 0.5) * (x - 1.0) * y * (y - 0.5),
//             7 => -8.0 * x * (x - 1.0) * y * (y - 0.5),
//             8 => 4.0 * x * (x - 0.5) * y * (y - 0.5),
//             _ => throw new ArgumentOutOfRangeException(nameof(function), function,
//                 $"Function {function} doesn't match interval [0, {BasisSize - 1}]")
//         };
//     }
//
//     public double DPhi(int function, int variable, double x, double y)
//     {
//         return variable switch
//         {
//             0 => function switch
//             {
//                 0 => 4.0 * (x - 1.0 + (x - 0.5)) * (y - 0.5) * (y - 1.0),
//                 1 => -8.0 * (x - 1.0 + x) * (y - 0.5) * (y - 1.0),
//                 2 => 4.0 * (x - 0.5 + x) * (y - 0.5) * (y - 1.0),
//                 3 => -8.0 * (x - 1.0 + (x - 0.5)) * y * (y - 1.0),
//                 4 => 16.0 * (x - 1.0 + x) * y * (y - 1.0),
//                 5 => -8.0 * (x - 0.5 + x) * y * (y - 1.0),
//                 6 => 4.0 * (x - 1.0 + (x - 0.5)) * y * (y - 0.5),
//                 7 => -8.0 * (x - 1.0 + x) * y * (y - 0.5),
//                 8 => 4.0 * (x - 0.5 + x) * y * (y - 0.5),
//                 _ => throw new ArgumentOutOfRangeException(nameof(function), function,
//                     $"Function {function} doesn't match interval [0, {BasisSize - 1}]")
//             },
//             1 => function switch
//             {
//                 0 => 4.0 * (x - 0.5) * (x - 1.0) * (y - 1.0 + (y - 0.5)),
//                 1 => -8.0 * x * (x - 1.0) * (y - 1.0 + (y - 0.5)),
//                 2 => 4.0 * x * (x - 0.5) * (y - 1.0 + (y - 0.5)),
//                 3 => -8.0 * (x - 0.5) * (x - 1.0) * (y - 1.0 + y),
//                 4 => 16.0 * x * (x - 1.0) * (y - 1.0 + y),
//                 5 => -8.0 * x * (x - 0.5) * (y - 1.0 + y),
//                 6 => 4.0 * (x - 0.5) * (x - 1.0) * (y - 0.5 + y),
//                 7 => -8.0 * x * (x - 1.0) * (y - 0.5 + y),
//                 8 => 4.0 * x * (x - 0.5) * (y - 0.5 + y),
//                 _ => throw new ArgumentOutOfRangeException(nameof(function), function,
//                     $"Function {function} doesn't match interval [0, {BasisSize - 1}]")
//             },
//             _ => throw new ArgumentOutOfRangeException(nameof(variable), variable,
//                 $"Function {function} doesn't match interval [0, 1]")
//         };
//     }
// }