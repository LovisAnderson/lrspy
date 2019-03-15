import lrs


class Bland(lrs.Lrs):
    def __init__(self, inequality_matrix, m, d):
        super().__init__(inequality_matrix, m, d)


    def select_pivot(self):
        i = 0
        j = self.d
        for k, c in enumerate(self.C[:-1]):
            if self.matrix[0][self.Column[k]] > 0:
                j = k
                break
        if j == self.d:
            return i, j
        min = 1e+10
        for k, b in enumerate(self.B[self.d:]):
            if self.matrix[self.Row[k+self.d]][self.Column[j]] == 0:
                continue
            ratio = self.matrix[self.Row[k+self.d]][0] / self.matrix[self.Row[k+self.d]][self.Column[j]]
            if abs(ratio) < min:
                min = abs(ratio)
                i = k + self.d
        return i, j

    def necessaryConditionForReverse(self):
        if self.matrix[self.Row[self.i]][self.Column[self.j]] == 0:
            return False
        return True
