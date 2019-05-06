import lrs


class Bland(lrs.Lrs):
    def __init__(self, hyperplane_matrix, m, d, bounding_box=None):
        super().__init__(hyperplane_matrix, m, d, bounding_box=bounding_box)

    def first_basis(self):
        self.variables_into_basis()
        if self.boxed:
            self.resort_inequalities() # todo necessary?
            self.move_to_box_corner()

        self.make_feasible()
        self.set_objective()
        self.resort_inequalities()
        self.append_solution()

    def select_pivot(self):
        i = 0
        j = self.d
        for k, c in enumerate(self.C[:-1]):
            if self.matrix[0][self.C.order[k]] > 0:
                j = k
                break
        if j == self.d:
            return i, j
        min = 1e+130
        min_sign = 'plus'
        for k, b in enumerate(self.B[self.d:]):
            m_k_j = self.matrix[self.B.order[k+self.d]][self.C.order[j]]
            if m_k_j == 0:
                continue
            if self.boxed and not self.pivot_stays_in_box(k + self.d, j):
                continue
            ratio = self.matrix[self.B.order[k+self.d]][0] / m_k_j
            if ratio > 0 and min_sign == 'minus':
                continue
            elif ratio < 0 and min_sign == 'plus':
                min = abs(ratio)
                min_sign = 'minus'
                i = k + self.d
            elif abs(ratio) < min:
                min = abs(ratio)
                i = k + self.d
        return i, j

    def necessary_condition_for_reverse(self):
        if self.boxed and not self.pivot_stays_in_box(self.i, self.j):
            return False
        if self.matrix[self.B.order[self.i]][self.C.order[self.j]] == 0:
            return False
        return True
