"""
An implementation of the Todd-Coxeter Algorithm
author: J. Simon Richard (jsimonrichard@gmail.com)
Explanation and Interactive examples: https://jsimonrichard.com/blog/todd-coxeter

A reorganized version of https://math.berkeley.edu/~kmill/notes/tc.py
Explanation: https://math.berkeley.edu/~kmill/notes/todd_coxeter.html
"""


import argparse
from typing import Callable, List


class Permutation:
  def __init__(self, perm: List[int]):
    self.cycles = []

    for i in range(len(perm)):
      cycle = [i]
      k = perm[i]
      while k != i:
        if k < i: break # We've already covered this cycle, so don't add
        cycle.append(k)
        k = perm[k]
      else:
        self.cycles.append(cycle)
  
  def __str__(self) -> str:
    return "(" + ")(".join(
      [" ".join([str(i+1) for i in cycle]) for cycle in self.cycles]
    ) + ")"


class UnionFind:
  def __init__(self):
    self.parent = []

  def newRoot(self) -> int:
    c = len(self.parent)
    self.parent.append(c)
    return c
  
  def getRoots(self) -> List[int]:
    return [c for i,c in enumerate(self.parent) if i == c]
  
  def find(self, a: int) -> int:
    if self.parent[a] == a:
      return a
    else:
      return self.find(self.parent[a])
  
  def union(self, a: int, b: int):
    # Get roots of trees
    a = self.find(a)
    b = self.find(b)

    # Re order
    a, b = min(a,b), max(a,b)

    # Unify (setting the lower value as the root)
    self.parent[b] = a
  
  def dynLoopThroughRoots(self, action: Callable[[int], None]):
    # While loop used since labels may be added to self.parents during iteration
    i = 0
    while i < len(self.parent):
      if i == self.find(i): # It's a root
        action(i)
      i += 1
  

class PartialSchreierGraph:
  def __init__(self, ngens: int):
    self.labels = UnionFind()

    self.adjList = []
    self.gens = list(range(1, ngens+1))
  
  def addNode(self) -> int:
    d = {}
    for g in self.gens:
      d[g] = None
      d[-g] = None

    self.adjList.append(d)
    return self.labels.newRoot()
  
  def followWord(self, node: int, word: List[int]) -> int:
    # Get lowest representative of this node
    node = self.labels.find(node)

    for g in reversed(word): # Reverse since we're doing the left cosets
      # Create new node if it doesn't exist
      if self.adjList[node][g] == None:
        self.adjList[node][g] = self.addNode()
      
      # Get lowest representative
      node = self.labels.find(self.adjList[node][g])
    
    return node

  def unifyNodes(self, a: int, b: int):
    # Find lowest representatives
    a = self.labels.find(a)
    b = self.labels.find(b)

    # Return if they're already the same
    if a == b:
      return
    
    # Order
    a, b = min(a,b), max(a,b)

    # Unify labels
    self.labels.union(a,b)

    # Unify graph connections
    for g in self.gens + [-g for g in self.gens]: # Include inverses
      # Just move the connection
      if self.adjList[a][g] == None:
        self.adjList[a][g] = self.adjList[b][g]
      
      # Two more nodes need to be unified
      elif self.adjList[b][g] != None:
        self.unifyNodes(
          self.adjList[a][g],
          self.adjList[b][g]
        )
    
  def dynLoopThroughNodes(self, action: Callable[[int], None]):
    self.labels.dynLoopThroughRoots(action)

  
class SchreierGraph(PartialSchreierGraph):
  """
  Builds a schreier graph for G/H using the Todd-Coxeter algorithm
  """
  def __init__(self, ngens, rels, hgens):
    super().__init__(ngens)
    self.rels = rels
    self.hgens = hgens

    # Step one
    basepoint = self.addNode()

    # Step two (unify using the generators of H)
    for hgen in hgens:
      endpoint = self.followWord(basepoint, hgen)
      self.unifyNodes(basepoint, endpoint)

    # Step three (follow the relations from each node and unify accordingly)
    trivial_rels = [[g,-g] for g in self.gens] # Include the trivial relation g*g^-1
    def stepThree(node: int):
      for rel in trivial_rels + rels:
        endpoint = self.followWord(node, rel)
        self.unifyNodes(node, endpoint)
    # Nodes may change during the iteration, so this function is used
    self.dynLoopThroughNodes(stepThree)
  
  def getPermutationRepresentation(self) -> List[Permutation]:
    """
    In this function we essentially realize Cayley's theorem for the G-set G/H (although
    it's not necessarily a group)
    """
    cosets = self.labels.getRoots()

    perms = []
    for g in self.gens:
      # See how g changes each coset
      perm = []
      for coset in cosets:
        # We're using cosets.index to essentially rename each coset to use the lowest
        # label possible
        perm.append(cosets.index(self.followWord(coset, [g])))
      perms.append(Permutation(perm))
    
    return perms


def parseWord(str) -> List[int]:
  """
  Format input word string into tuple
  """
  return [int(s) for s in str.split(",")]


if __name__=="__main__":
  # Parse command line arguments
  parser = argparse.ArgumentParser(description='Enumerate the left cosets of H in G using the Todd-Coxeter algorithm')
  required = parser.add_argument_group("required arguments")
  required.add_argument('--ngens', metavar='N', type=int, required=True,
                      help="""number of generators of G (generators will be labeled 1...n
                      and their inverses will be labeled -1...-n)""")
  required.add_argument('--rels', metavar="a,b,c", type=parseWord, nargs="+", required=True,
                      help="""relations of the form a,b,c where a*b*c = e""")
  required.add_argument('--hgens', metavar="a,b,c", type=parseWord, nargs="+", required=True,
                      help="""generators of H (of the form a,b,c meaning a*b*c)""")

  args = parser.parse_args()


  graph = SchreierGraph(args.ngens, args.rels, args.hgens)

  print("g = permutation representation")
  print("------------------------------")
  for i, perm in enumerate(graph.getPermutationRepresentation()):
    print(f"{i+1} = {perm}")