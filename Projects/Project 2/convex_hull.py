from which_pyqt import PYQT_VER


if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF, QObject
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF, QObject
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time


# Some global color constants that might be useful
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)

# Global variable that controls the speed of the recursion automation, in seconds
#
PAUSE = 0.25

#
# This is the class you have to complete.
#
class ConvexHullSolver(QObject):

    # Class constructor
    def __init__(self):
        super().__init__()
        self.pause = False

    # Some helper methods that make calls to the GUI, allowing us to send updates
    # to be displayed.

    def showTangent(self, line, color):
        self.view.addLines(line, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseTangent(self, line):
        self.view.clearLines(line)

    def blinkTangent(self, line, color):
        self.showTangent(line, color)
        self.eraseTangent(line)

    def showHull(self, polygon, color):
        self.view.addLines(polygon, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseHull(self, polygon):
        self.view.clearLines(polygon)

    def showText(self, text):
        self.view.displayStatusText(text)

    # This is the method that gets called by the GUI and actually executes
    # the finding of the hull
    def compute_hull(self, points, pause, view):
        self.pause = pause
        self.view = view
        assert (type(points) == list and type(points[0]) == QPointF)

        t1 = time.time()

        # Sorting the array of points, by horizontal position, done using TimSort
        # Time:O(nlog(n)) Space:O(n)
        def a(pt):
            return pt.x
        points.sort(key=lambda pt: pt.x())

        t2 = time.time()

        t3 = time.time()
        # Time:O(nlog(n)) Space: O(nlog(n))
        hull = ConvexHull.from_points(points)

        t4 = time.time()

        # when passing lines to the display, pass a list of QLineF objects.  Each QLineF
        # object can be created with two QPointF objects corresponding to the endpoints
        self.showHull(hull.to_polygon(), RED)
        print(t4 - t3)
        self.showText('Time Elapsed (Convex Hull): {:3.3f} sec'.format(t4 - t3))


'''
               Empyrical Analysis
       10:   0.001   0.000   0.001   0.000   0.001
      100:   0.007   0.006   0.006   0.006   0.002
    1,000:   0.066   0.050   0.058   0.066   0.051
   10,000:   0.228   0.223   0.222   0.206   0.211
  100,000:   2.261   2.235   2.217   2.226   2.233
  500,000:  11.586  11.345  11.391  11.350  11.536
1,000,000:  23.123  23.020  25.778  23.119  23.113
'''

class ConvexHull:
    """
    Class used to store the ConvexHull. Each Hull is made from a linked list of HullNodes [inner class within Convex Hull]
    storing both the left-most node and the right-most, both of which are the same linked list

    The init method takes in 2 Nodes, the left-most and the right-most node
    There are two class methods, the first from_points creates a hull from a list of points, ordered by x-values
    The other from_hull merges two separate hulls one on the left side and one on the right.

    Both of these class methods are recursive, and the from_points method calls from_hull on the way up the recursion

    Finally, the to_polygon method returns a list of QLineF's that form the edges of the convex hull
    """
    class HullNode:
        """
        Inner class, used to store each node in the hull. Each node consists of a QPointF point and a reference to
        the clockwise[cw] and counter-clockwise[ccw] neighbor
        """
        def __init__(self, pt: QPointF):
            """
            Assigns the cw and ccw neighbors to itself
                Time: O(1)
                Space: O(1)

            :param pt: location of the node
            :rtype: ConvexHull.HullNode
            """
            assert (type(pt) == QPointF)
            self.pt = QPointF(pt.x(), pt.y())
            self.cw = self
            self.ccw = self

        def slope_to(self, other) -> float:
            """
            Find the slope of the node relative to another
                Time: O(1)
                Space: O(1)

            :param other: Node to which to find the relative slope to
            :type other: ConvexHull.HullNode
            :return: the slope of the two nodes
            """
            assert (type(other) == ConvexHull.HullNode)
            return (self.pt.y() - other.pt.y()) / (self.pt.x() - other.pt.x())

    def __init__(self, lPt: HullNode, rPt: HullNode):
        """
        Create a ConvexHull, providing the HullNode for the left-most and right-most sides
            Time: O(1)
            Space: O(1)

        :param lPt: The leftmost HullNode, a linked-list connected to the rightmost HullNode
        :param rPt: The rightmost HullNode, a linked-list connected to the leftmost HullNode
        :rtype: ConvexHull
        """
        self.leftMostPt = lPt
        self.rightMostPt = rPt

    @classmethod
    def from_points(cls, points: [QPointF]):
        """
        This takes a list of QPointsF's, sorted by x value and recursively creates a ConvexHull by splitting the list
        of points into two equal portions, until the size of the list equals one. These lists of size one point are
        used to create small ConvexHulls, which are then merged back together using the from_hull class method.
            Time: O(nlog(n))
            Space: O(nlog(n))
        :param points: A list of sorted QPoints by x-value
        :return: A ConvexHull containing all the provided points :rtype: ConvexHull
        """
        n = len(points)
        if n > 1:
            return cls.from_hull(cls.from_points(points[:n // 2]), cls.from_points(points[n // 2:]))
        else:
            node = cls.HullNode(points[0])
            return cls(node, node)

    @classmethod
    def from_hull(cls, left, right):
        """
        Creates a ConvexHull created by merging two separate ConvexHulls (with size n and m)
            Time: O(n+m)
            Space: O(1)
        :param p: ConvexHull on the left-side :type: ConvexHull
        :param q: ConvexHull on the right-side :type: ConvexHull
        :return: The newly merged ConvexHUll :rtype: ConvexHull
        """
        assert (type(left) == ConvexHull and type(right) == ConvexHull)

        def right_side_upper_tangent(p, q):
            """
            Fings the upper tangent HullNode of the right side
                Time: O(n)
                Space: O(1)
            :param p: Current tangent HullNode of right-side ConvexHull
            :param q: Current tangent HullNode of left-side ConvexHull
            :return: The uppermost HullNode for the right-side, and a boolean value on whether the tangent was updated
            """
            updated = True
            while p.slope_to(q) > p.ccw.slope_to(q):
                p, updated = p.ccw, False
            return p, q, updated

        def left_side_upper_tangent(p, q):
            """
            Fings the upper tangent HullNode of the left side
                Time: O(m)
                Space: O(1)
            :param p: Current tangent HullNode of right-side ConvexHull
            :param q: Current tangent HullNode of left-side ConvexHull
            :return: The uppermost HullNode for the left-side, and a boolean value on whether the tangent was updated
            """
            updated = True
            while q.slope_to(p) < q.cw.slope_to(p):
                q, updated = q.cw, False
            return p, q, updated

        def right_side_lower_tangent(p, q):
            """
            Fings the lower tangent HullNode of the right side
                Time: O(n)
                Space: O(1)
            :param p: Current tangent HullNode of right-side ConvexHull
            :param q: Current tangent HullNode of left-side ConvexHull
            :return: The lowermost HullNode for the right-side, and a boolean value on whether the tangent was updated
            """
            updated = True
            while p.slope_to(q) < p.cw.slope_to(q):
                p, updated = p.cw, False
            return p, q, updated

        def left_side_lower_tangent(p, q):
            """
            Fings the upper tangent HullNode of the left side
                Time: O(m)
                Space: O(1)
            :param p: Current tangent HullNode of right-side ConvexHull
            :param q: Current tangent HullNode of left-side ConvexHull
            :return: The uppermost HullNode for the left-side, and a boolean value on whether the tangent was updated
            """
            updated = True
            while q.slope_to(p) > q.ccw.slope_to(p):
                q, updated = q.ccw, False
            return p, q, updated

        def merge_from_tangents():
            """
            This function merges the two ConvexHulls based on their upper and lower tangents between the two hulls.
                Time: O(n+m) - The time is most complex, despite this function calls four helper functions each with
                    O(n/m). This function will at the worst case visit every HullNode within the two hulls, that is
                    because the upper functions visit different nodes than the lower functions.
                Space: O(1)
            """
            # Find Upper Tangent
            upL, upR = left.rightMostPt, right.leftMostPt
            left_found, right_found = False, False
            while not left_found or not right_found:
                upL, upR, left_found = left_side_upper_tangent(upL, upR)
                upL, upR, right_found = right_side_upper_tangent(upL, upR)
            # Find Lower Tangent
            downL, downR = left.rightMostPt, right.leftMostPt
            left_found, right_found = False, False
            while not left_found or not right_found:
                downL, downR, left_found = left_side_lower_tangent(downL, downR)
                downL, downR, right_found = right_side_lower_tangent(downL, downR)



            upL.cw, upR.ccw = upR, upL
            downL.ccw, downR.cw = downR, downL

        merge_from_tangents()
        return cls(left.leftMostPt, right.rightMostPt)

    # Time: O(n) Space: O(n)
    def to_polygon(self):
        """
        Creates a list of QLineF's created by connecting the HullNodes clockwise, starting from the rightmost HullNode.
        This loops until we arrive at back the rightmost node forming a closed polygon
            Time: O(n)
            Space: O(n)
        :return: A list of QLineF's formed by polygon of the ConvexHull
        """
        start, curr = self.rightMostPt, self.rightMostPt
        poly = [QLineF(start.ccw.pt, start.pt)]
        while curr.cw.pt.x() != start.pt.x():
            poly.append(QLineF(curr.pt, curr.cw.pt))
            curr = curr.cw
        return poly
