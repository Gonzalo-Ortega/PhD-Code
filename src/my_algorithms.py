
from sage.all import *

def my_fundamental_group(complex: SimplicialComplex, base_point=None, simplify=True):
    # Check for connectivity:
    if not complex.is_connected():
        if base_point is None:
            raise ValueError("this complex is not connected, so you must specify a base point")
        return complex.connected_component(Simplex([base_point])).fundamental_group(simplify=simplify)
    
    from sage.groups.free_group import FreeGroup
    from sage.libs.gap.libgap import libgap as gap

    # Get a spannig tree:
    G = complex.graph()
    spanning_tree = {frozenset((u, v)) for u, v, _ in G.min_spanning_tree()}

    # Get generators (edges not in spanning tree):
    gens = [e for e in G.edge_iterator(labels=False)
                if frozenset(e) not in spanning_tree]
    if not gens:
        return gap.TrivialGroup()
    
    # Generate generators dictionary for later consult:
    gens_dict = {frozenset(g): i for i, g in enumerate(gens)}

    # Make free group out from generators:
    FG = FreeGroup(len(gens), 'e')

    # Get relators:
    rels = []
    for cell in complex._n_cells_sorted(2): # Iterate over 2-simplices.
        boundary = [tuple(e) for e in cell.faces()]
        z = {}
        for i in range(3): # Iterate over simplex edges.
            x = frozenset(boundary[i])
            if x in spanning_tree:
                z[i] = FG.one() # Add a one if edge is in spannig tree.
            else:
                z[i] = FG.gen(gens_dict[x]) # Add generator if edge is not.
        rels.append(z[0]*z[1].inverse()*z[2]) # Add combination to relators.

    # Symplify if desired
    if simplify:
        return FG.quotient(rels).simplified()
    else:
        return FG.quotient(rels)