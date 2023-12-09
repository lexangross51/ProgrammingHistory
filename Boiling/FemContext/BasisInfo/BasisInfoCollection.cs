using System.Collections;

namespace Boiling.FemContext.BasisInfo;

public readonly struct BasisInfoItem(int functionNumber, BasisFunctionType type, int index)
{
    public int FunctionNumber { get; } = functionNumber;
    public BasisFunctionType Type { get; } = type;
    public int Index { get; } = index;
}

public class BasisInfoCollection : IEnumerable<KeyValuePair<int, BasisInfoItem[]>>
{
    private readonly Dictionary<int, BasisInfoItem[]> _basisByElements;
    private int _functionsCount = -1;

    public int FunctionsCount
    {
        get
        {
            if (_functionsCount == -1 && _basisByElements.Count != 0)
            {
                _functionsCount = _basisByElements.Values
                    .SelectMany(pair => pair)
                    .Max(item => item.FunctionNumber) + 1;
            }

            return _functionsCount;
        }
    }
    
    public BasisInfoCollection(int elementCount, int functionsPerElement)
    {
        _basisByElements = new Dictionary<int, BasisInfoItem[]>(elementCount * functionsPerElement);

        for (int i = 0; i < elementCount; i++)
        {
            _basisByElements.Add(i, new BasisInfoItem[functionsPerElement]);
        }
    }
    
    public BasisInfoItem this[int elementIndex, int functionIndex]
    {
        get => _basisByElements[elementIndex][functionIndex];
        set => _basisByElements[elementIndex][functionIndex] = value;
    }

    public BasisInfoItem GetFunctionAtNode(int ielem, int node)
    {
        return _basisByElements[ielem]
            .FirstOrDefault(b => b.Type == BasisFunctionType.ByGeometricNode && b.Index == node);
    }
    
    public BasisInfoItem GetFunctionAtEdge(int ielem, int edge)
    {
        return _basisByElements[ielem]
            .FirstOrDefault(b => b.Type == BasisFunctionType.ByEdgeNode && b.Index == edge);
    }
    
    public IEnumerator<KeyValuePair<int, BasisInfoItem[]>> GetEnumerator()
        => _basisByElements.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
}