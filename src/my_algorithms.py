import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection



from sage.all import *
import gudhi as gd


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


def plot_points(points, n):
    """
    Plot points based on their dimension (n).
    - n = 1: Scatter plot on x-axis
    - n = 2: 2D scatter plot
    - n = 3: 3D scatter plot
    """
    if n == 1:
        x = [p[0] for p in points]
        plt.scatter(x, [0]*len(points), c='blue')
        plt.xlabel('X')
        plt.title('1D Points')
        plt.yticks([])
        plt.show()

    elif n == 2:
        x = [p[0] for p in points]
        y = [p[1] for p in points]
        plt.scatter(x, y, c='green')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('2D Points')
        plt.show()

    elif n == 3:
        x = [p[0] for p in points]
        y = [p[1] for p in points]
        z = [p[2] for p in points]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c='red')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('3D Points')
        plt.show()
    else:
        raise ValueError("Plotting is only supported for n = 1, 2, or 3")





def plot_filtration(simplices, points, dim=2, steps=None, labels=None):
    """
    Plot a simplicial complex filtration.

    Parameters:
    simplices: list of tuples (list of vertex indices, filtration value)
    points: np.array of shape (n_points, dim)
    dim: 2 or 3 for 2D or 3D plot
    steps: None (full complex), True (all steps), or int (limit steps)
    labels: None, 'step', 'value', or 'both' for labeling simplices
    """
    # Group simplices by filtration value
    filtration_dict = {}
    for simplex, value in simplices:
        filtration_dict.setdefault(value, []).append(simplex)

    # Sort filtration values
    filtration_values = sorted(filtration_dict.keys())

    # Prepare steps
    grouped_steps = []
    accumulated = []
    for val in filtration_values:
        accumulated.extend([(s, val) for s in filtration_dict[val]])
        grouped_steps.append(list(accumulated))

    if steps is None:
        steps_to_plot = [grouped_steps[-1]]
    else:
        max_steps = len(grouped_steps) if steps is True else min(steps, len(grouped_steps))
        steps_to_plot = grouped_steps[:max_steps]

    for step_idx, step in enumerate(steps_to_plot):
        fig = plt.figure(figsize=(8, 6))
        if dim == 2:
            ax = fig.add_subplot(111)
            ax.scatter(points[:, 0], points[:, 1], color='black')
            for simplex, val in step:
                coords = points[simplex]
                if len(simplex) == 1:
                    ax.scatter(coords[0, 0], coords[0, 1], color='blue')
                    if labels:
                        label_text = _get_label(labels, step_idx+1, val)
                        ax.text(coords[0, 0], coords[0, 1], label_text, fontsize=8)
                elif len(simplex) == 2:
                    p1, p2 = coords
                    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color='green')
                    if labels:
                        center = (p1 + p2) / 2
                        label_text = _get_label(labels, step_idx+1, val)
                        ax.text(center[0], center[1], label_text, fontsize=8)
                elif len(simplex) == 3:
                    ax.fill(coords[:, 0], coords[:, 1], alpha=0.3, color='orange')
                    if labels:
                        center = coords.mean(axis=0)
                        label_text = _get_label(labels, step_idx+1, val)
                        ax.text(center[0], center[1], label_text, fontsize=8)
        elif dim == 3:
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(points[:, 0], points[:, 1], points[:, 2], color='black')
            for simplex, val in step:
                coords = points[simplex]
                if len(simplex) == 1:
                    ax.scatter(coords[0, 0], coords[0, 1], coords[0, 2], color='blue')
                    if labels:
                        label_text = _get_label(labels, step_idx+1, val)
                        ax.text(coords[0, 0], coords[0, 1], coords[0, 2], label_text, fontsize=8)
                elif len(simplex) == 2:
                    p1, p2 = coords
                    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color='green')
                    if labels:
                        center = (p1 + p2) / 2
                        label_text = _get_label(labels, step_idx+1, val)
                        ax.text(center[0], center[1], center[2], label_text, fontsize=8)
                elif len(simplex) == 3:
                    poly = Poly3DCollection([coords], alpha=0.3, facecolor='orange')
                    ax.add_collection3d(poly)
                    if labels:
                        center = coords.mean(axis=0)
                        label_text = _get_label(labels, step_idx+1, val)
                        ax.text(center[0], center[1], center[2], label_text, fontsize=8)
                elif len(simplex) == 4:
                    faces = [[coords[i] for i in face] for face in [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]]
                    poly = Poly3DCollection(faces, alpha=0.2, facecolor='purple')
                    ax.add_collection3d(poly)
                    if labels:
                        center = coords.mean(axis=0)
                        label_text = _get_label(labels, step_idx+1, val)
                        ax.text(center[0], center[1], center[2], label_text, fontsize=8)
        ax.set_title(f"Filtration Step {step_idx+1}" if steps is not None else "Full Complex")
        plt.show()

def _get_label(option, step, value):
    if option == 'step':
        return f"Step {step}"
    elif option == 'value':
        return f"Val {value:.3f}"
    elif option == 'both':
        return f"Step {step}\nVal {value:.3f}"
    return ""


def get_random_filtration():

    return

def manual_persitent_fundamental_group():
    
    return

#filtration = gd.SimplexTree()


