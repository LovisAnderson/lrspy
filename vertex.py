class Vertex:
    def __init__(self, coordinates, cobasis=None):
        self.coordinates = tuple(coordinates)
        self.cobasis = cobasis
        self.position_vector = None

    def compute_position_vector(self, matrix, B, nr_hyperplanes, d):
        self.position_vector = [0] * nr_hyperplanes
        for i, b in enumerate(B[d:]):
            if b.slack_variable and not b.box_variable:
                if matrix[B.order[d+i]][0] > 0:
                    self.position_vector[b.hyperplane_index] = 1
                elif matrix[B.order[d+i]][0] < 0:
                    self.position_vector[b.hyperplane_index] = -1
                else:
                    self.position_vector[b.hyperplane_index] = 0
                    self.cobasis.append(b)

    def cobasis_lexicographic_minimal(self):
        if len(self.cobasis) == len(self.coordinates):
            return True
        for i in range(len(self.cobasis) - 1):
            if self.cobasis[i+1] < self.cobasis[i]:
                return False
        return True

    def __str__(self):
        return str(self.coordinates)

    def __eq__(self, other):
        if self.coordinates == other.coordinates:
            return True
        else:
            return False

    def __hash__(self):
        return hash(self.coordinates)
