using Boiling.FemContext.BasisInfo.Interfaces;

namespace Boiling.FemContext.BasisInfo;

public readonly record struct BiLinearBasis : IBasis
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