namespace Boiling.FemContext.BasisInfo.Interfaces;

public interface IBasis
{
    int BasisSize { get; }
    double Phi(int function, double x, double y);
    double DPhi(int function, int variable, double x, double y);
}