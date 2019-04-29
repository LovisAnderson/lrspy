import numpy as np
import matplotlib.pyplot as plt


def __intersect(rect, line):
    l = []
    xmin,xmax,ymin,ymax = rect
    a,b,c = line

    assert a!=0 or b!=0

    if a == 0:
        y = -c/b
        if y<=ymax and y>=ymin:
            l.append((xmin, y))
            l.append((xmax, y))
        return l
    if b == 0:
        x = -c/a
        if x<=xmax and x>=xmin:
            l.append((x, ymin))
            l.append((x, ymax))
        return l

    k = -a/b
    m = -c/b
    for x in (xmin, xmax):
        y = k*x+m
        if y<=ymax and y>= ymin:
            l.append((x,y))

    k = -b/a
    m = -c/a
    for y in (ymin, ymax):
        x = k*y+m
        if x<xmax and x> xmin:
            l.append((x,y))
    return l


def plotline(coef, *args, **kwargs):
    '''plot line: y=a*x+b or a*x+b*y+c=0'''
    coef = np.float64(coef[:])
    assert len(coef)==2 or len(coef)==3
    if len(coef) == 2:
        a, b, c = coef[0], -1., coef[1]
    elif len(coef) == 3:
        a, b, c = coef
    ax = plt.gca()

    limits = ax.axis()
    points = __intersect(limits, (a,b,c))
    if len(points) == 2:
        pts = np.array(points)
        ax.plot(pts[:,0], pts[:,1], *args, **kwargs)
        ax.axis(limits)


def circle_out(x, y, s=20, *args, **kwargs):
    '''Circle out points with size 's' and edgecolors'''
    ax = plt.gca()
    if 'edgecolors' not in kwargs:
        kwargs['edgecolors'] = 'g'
    ax.scatter(x, y, s, facecolors='none', *args, **kwargs)


def plot_arrangement(arrangement_matrix, bounding_box=[], ax=None, x_limits=(0, 4), y_limits=(0, 4), point=None):
    colors = well_distinguishable_colors(len(arrangement_matrix) + len(bounding_box) + 1)
    if ax is None:
        show = True
        ax = plt.gca()
    else:
        show = False
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    for i, hyperplane in enumerate(arrangement_matrix + bounding_box):
        line = hyperplane[1:] + hyperplane[0:1]
        if i < len(arrangement_matrix):
            plotline(line, color=colors[i])
        else:
            plotline(line, color=colors[i], linestyle='dashed')
    if point is not None:
        circle_out(point[0], point[1])
    if show:
        plt.show()


def well_distinguishable_colors(number_colors):
    hsv = plt.get_cmap('hsv')
    colors = hsv(np.linspace(0, 1.0, number_colors))
    return colors
