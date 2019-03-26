from lrs_datastructures import LrsDict, Variable


def test_dict_sort():
    a = LrsDict()
    a.append(1)
    a.append(4)
    a.append(2)
    a.order = [0, 1, 2]
    a.sort_respecting_order()
    assert a == [1, 2, 4]
    assert a.order == [0, 2, 1]