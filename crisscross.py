import lrs


class CrissCross(lrs.Lrs):
    def __init__(self, hyperplane_matrix, m, d, bounding_box=None):
        super().__init__(hyperplane_matrix, m, d, bounding_box=bounding_box)

    def select_pivot(self):
        basis_index = self.d
        cobasis_index = 0
        for i in range(self.d, self.d + self.m):
            if basis_index <= self.m and self.B[basis_index] == i:
                if self.matrix[self.B.order[basis_index]][0] < 0:
                    # Bi is primal infeasible
                    for cobasis_index, c in enumerate(self.C[:-1]):
                        if self.boxed and not self.pivot_stays_in_box(basis_index, cobasis_index):
                            continue
                        if self.matrix[self.B.order[basis_index]][self.C.order[cobasis_index]] > 0:
                            print('Primal infeasible! pivot i={}, j={}'.format(basis_index,
                                                                               cobasis_index))
                            return basis_index, cobasis_index
                    if not self.boxed:
                        raise ValueError
                basis_index += 1
            elif cobasis_index < self.d and self.C[cobasis_index] == i:
                if self.matrix[0][self.C.order[cobasis_index]] > 0:
                    # Ci is dual infeasible
                    for basis_index, b in enumerate(self.B):
                        if b < self.d:
                            continue
                        if self.boxed and not self.pivot_stays_in_box(basis_index, cobasis_index):
                            continue
                        if self.matrix[self.B.order[basis_index]][self.C.order[cobasis_index]] < 0:
                            print('Dual infeasible! pivot i={}, j={}'.format(basis_index,
                                                                               cobasis_index))
                            return basis_index, cobasis_index
                    if not self.boxed:
                        # todo Is it possible to have dual infeasible solutions with no valid pivot?
                        raise ValueError
                cobasis_index += 1
        return 0, 0

    def necessary_condition_for_reverse(self):
        def lower_index_pivot(k, k_in_cobasis=False):
            i = self.i if k_in_cobasis else k
            j = k if k_in_cobasis else self.j
            m_i = self.B.order[self.i] if k_in_cobasis else self.B.order[k]
            m_j = self.C.order[k] if k_in_cobasis else self.C.order[self.j]
            pivot_element = self.matrix[m_i][m_j]
            if (k_in_cobasis and pivot_element < 0) or (not k_in_cobasis and pivot_element > 0):
                if not self.boxed or self.pivot_stays_in_box(i, j):
                    return True
            return False

        if self.boxed and not self.pivot_stays_in_box(self.i, self.j):
            return False

        if self.matrix[self.B.order[self.i]][0] > 0:
            if (
                    self.matrix[self.B.order[self.i]][self.C.order[self.j]] > 0 and
                    not any(lower_index_pivot(k, k_in_cobasis=True)
                            for k in
                            range(0, self.max_index_of_smaller_number(self.C, self.B[self.i]) + 1)
                            )
            ):
                return True
        if self.matrix[0][self.C.order[self.j]] < 0:
            if (self.matrix[self.B.order[self.i]][self.C.order[self.j]] < 0 and
                not any(lower_index_pivot(k, k_in_cobasis=False)
                        for k in
                        range(self.d, self.max_index_of_smaller_number(self.C, self.C[self.j]) + 1)
                        )):
                return True
        print('Necessary Condition for reverse not fulfilled!')
        return False