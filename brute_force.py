from ppl import Constraint_System, C_Polyhedron, Variable
from vertex import Vertex
import itertools


def solve_les(hyperplanes):
    dim = len(hyperplanes[0]) - 1
    x = [Variable(i) for i in range(dim)]
    constraints = Constraint_System()
    for hyp in hyperplanes:
        constraints.insert(sum(hyp[i + 1] * x[i] for i in range(dim)) + hyp[0] == 0)
    poly = C_Polyhedron(constraints)
    ppl_points = [pt for pt in poly.minimized_generators() if pt.is_point()]
    if len(ppl_points) != len(poly.minimized_generators()):
        # Not uniquely determined.
        return None
    if len(ppl_points) == 1:
        vertex = Vertex(tuple(c / ppl_points[0].divisor() for c in ppl_points[0].coefficients()), None)
        return vertex
    elif len(ppl_points) == 0:
        return None
    else:
        raise ValueError


def inside_box(coordinates, box):
    for box_constraint in box:
        if sum(coordinates[i] * box_constraint[i+1] for i in range(len(coordinates))) + box_constraint[0] < 0:
            return False
    return True


def brute_force_vertices(hyperplanes, box_constraints):
    vertices = set()
    d = len(hyperplanes[0]) - 1
    for combination in itertools.combinations(range(len(hyperplanes)), d):
        hyperplane_combination = [hyperplanes[i] for i in combination]
        vertex = solve_les(hyperplane_combination)
        if vertex is not None:
            vertex.cobasis = combination
            vertices.add(vertex)

    vertices_inside_box = set()
    for v in vertices:
        if inside_box(v.coordinates, box_constraints):
            vertices_inside_box.add(v)
    return vertices, vertices_inside_box
