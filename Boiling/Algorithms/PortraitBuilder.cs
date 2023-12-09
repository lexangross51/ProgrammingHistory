using Boiling.FemContext.BasisInfo;

namespace Boiling.Algorithms;

public static class PortraitBuilder
{
    public static (int[] Ig, int[] Jg) GeneratePortrait(BasisInfoCollection basisInfo)
    {
        var funcCount = basisInfo.Last().Value[^1].FunctionNumber + 1;
        var connectivityList = new List<SortedSet<int>>();

        for (int i = 0; i < funcCount; i++)
        {
            connectivityList.Add(new SortedSet<int>());
        }

        foreach (var nodes in basisInfo)
        {
            foreach (var nodeToInsert in nodes.Value)
            {
                foreach (var posToInsert in nodes.Value)
                {
                    if (nodeToInsert.FunctionNumber < posToInsert.FunctionNumber)
                    {
                        connectivityList[posToInsert.FunctionNumber].Add(nodeToInsert.FunctionNumber);
                    }
                }
            }
        }

        var ig = new int[connectivityList.Count + 1];

        ig[0] = 0;
        ig[1] = 0;

        for (int i = 1; i < connectivityList.Count; i++)
        {
            ig[i + 1] = ig[i] + connectivityList[i].Count;
        }

        var jg = new int[ig[^1]];

        for (int i = 1, j = 0; i < connectivityList.Count; i++)
        {
            foreach (var it in connectivityList[i])
            {
                jg[j++] = it;
            }
        }

        return (ig, jg);
    }
}