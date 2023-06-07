# Visuals used in final project presentation

| Tool | Method | Run time | Space |
|------|--------|----------| ----------|
|BWA MEM| Maximal Exact Matches | Linear | Linear |
|Bowtie2| Dual Direction Backttracking | Almost Linear | Linear |
|STAR | Sequential Maximum Mappable Seed | Close to Linear | Linear |

| Tool | Method | Run time | Space |
|------|--------|----------| ----------|
|wf-align| BWA Backtracking | O(Reads Length * Reference Length) | Linear |

| Structure | Method | Run time | Space |
|------|--------|----------| ----------|
|Suffix Array| Prefix Doubling | O(nlogn) | O(n) |
|Burrows Wheeler Transformation | Directly Read from SA | O(n) | O(n) |
|Last to First Alignment| Radix Sort | O(n) | O(n) |
|Count Vector| Iterating through First Column of BWT Matrix | O(n) | O(1) |


| Method | Run time | Space |
|--------|----------| ----------|
| Exact Backtracking | O(m * n / (size of alphabet) ) | O(n / (size of alphabet)) |