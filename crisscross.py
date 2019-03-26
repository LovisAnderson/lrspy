import lrs


class CrissCross(lrs.Lrs):
    def __init__(self, hyperplane_matrix, m, d):
        super().__init__(hyperplane_matrix, m, d)

    def select_pivot(self):
        basis_index = self.d
        cobasis_index = 0
        for i in range(self.d, self.d + self.m):
            if basis_index <= self.m and self.B[basis_index] == i:
                if self.matrix[self.B.order[basis_index]][0] < 0:
                    # Bi is primal infeasible
                    print('B[{}] = {} primal infeasible'.format(basis_index, self.B[basis_index]))
                    for cobasis_index, c in enumerate(self.C):
                        if self.boxed and not self.pivot_stays_in_box(basis_index, cobasis_index):
                            continue
                        if self.matrix[self.B.order[basis_index]][self.C.order[cobasis_index]] > 0:
                            return basis_index, cobasis_index
                    raise ValueError
                basis_index += 1
            elif cobasis_index < self.d and self.C[cobasis_index] == i:
                if self.matrix[0][self.C.order[cobasis_index]] > 0:
                    # Ci is dual infeasible
                    print('C[{}] = {} dual infeasible'.format(cobasis_index, self.C[cobasis_index]))
                    for basis_index, b in enumerate(self.B):
                        if basis_index < self.d:
                            continue
                        if self.boxed and not self.pivot_stays_in_box(basis_index, cobasis_index):
                            continue
                        if self.matrix[self.B.order[basis_index]][self.C.order[cobasis_index]] < 0:
                            return basis_index, cobasis_index
                    raise ValueError
                cobasis_index += 1

    def necessaryConditionForReverse(self):
        if self.boxed:
            if not self.pivot_stays_in_box(self.i, self.j):
                print('Reverse Does not stay in box!')
                return False
            return self.matrix[self.B.order[self.i]][self.C.order[self.j]] != 0
        if self.matrix[self.B.order[self.i]][0] > 0:
            if (
                    self.matrix[self.B.order[self.i]][self.C.order[self.j]] > 0 and
                    all(
                        self.matrix[self.B.order[self.i]][self.C.order[k]] >= 0
                        for k in range(0, self.maxIndexSmallerNumber(self.C, self.B[self.i]) + 1)
                    )
            ):
                return  True
        if self.matrix[0][self.C.order[self.j]] < 0:
            if (self.matrix[self.B.order[self.i]][self.C.order[self.j]] < 0 and
                all(
                    self.matrix[self.B.order[k]][self.C.order[self.j]] <= 0
                    for k in range(1, self.maxIndexSmallerNumber(self.C, self.C[self.j]) + 1)
                )
            ):
                return True
        return False