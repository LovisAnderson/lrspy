class LrsDict(list):
    def __init__(self, iterable=[]):
        list.__init__(self, iterable)
        self.order = list(range(len(iterable))) # If Basis this gives rows in matrix / cobasis columns

    def sort_respecting_order(self):
        newDic, newLoc = zip(*sorted(zip(self, self.order)))
        self.sort()
        self.order = list(newLoc)


class Variable(int):
    def __new__(cls, integer):
        v = int.__new__(cls, integer)
        v.box_variable = False
        v.slack_variable = False
        v.hyperplane_index = None
        return v

    def change_variable(self, index):
        v = Variable(index)
        v.box_variable = self.box_variable
        v.slack_variable = self.slack_variable
        v.hyperplane_index = self.hyperplane_index
        return v
