import itertools

import lrs

import logging, logging_config
logger = logging.getLogger(__name__)


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
                            logger.detailed_debug('Primal infeasible! pivot i={}, j={}'.format(
                                basis_index, cobasis_index)
                            )
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
                            logger.detailed_debug('Dual infeasible! pivot i={}, j={}'.format(
                                basis_index, cobasis_index)
                            )
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
        logger.detailed_debug('Necessary Condition for reverse not fulfilled!')
        return False

    def forest_search(self):
        degenerated_basis_vars = self.get_degenericies()
        search_status = lrs.SearchStatus.NONE
        if len(degenerated_basis_vars) == 0:
            search = self.search()
            while search_status != lrs.SearchStatus.DONE:
                search_status = search.__next__()
                yield search_status
        else:
            logging.info('Degenerated start basis: Start forest search!')
            degeneracy = self.C[:-1] + degenerated_basis_vars
            degeneracy_hyperplanes = [v.hyperplane_index for v in degeneracy]
            start_bases_by_hyperplanes = itertools.combinations(degeneracy_hyperplanes, self.d - 1)
            self.append_solution(check_lexicographic_order=False)
            for basis_hyperplanes in start_bases_by_hyperplanes:
                logging.info(f'Start tree search with start hyperplanes {basis_hyperplanes}')
                self.move_to_hyperplanes(basis_hyperplanes)
                search_status = lrs.SearchStatus.NEWTREE
                yield search_status
                search = self.search()
                while True:
                    search_status = search.__next__()
                    if search_status == lrs.SearchStatus.DONE:
                        break
                    yield search_status
            yield lrs.SearchStatus.DONE
