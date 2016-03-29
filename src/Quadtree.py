# Credit to Steven Lambert for the tutorial

class Rectangle:
    """
    The Rectangle class defines the boundaries of a rectange used to bound an irregular shape
    """

    def __init__(self, x: float, y: float, width: float, height: float):
        """
        :param x: float representing the bottom left corner of the box
        :param y: float represting the bottom left corner of the box
        :param width: float representing width of the box
        :param height: float representing the width of the box
        """
        self.x = float(x)
        self.y = float(y)
        self.width = float(width)
        self.height = float(height)


class Quadtree:
    """
    The Quadtree class creates a tree for fast lookup of neighbors of a set of objects in a 2D
    space
    """

    # Adjust to increase total pool size or change maximum deepth
    max_objects = 10
    max_levels = 7

    def __init__(self, level: int, bounds: Rectangle):
        """
        :param level: indicating what level of the quadtree this particular node sits
        :param bounds: A Rectange defining the boundaries of this particular node
        """
        self.level = int(level)
        self.objects = []
        self.bounds = bounds
        self.members = 0
        self.nodes = [None, None, None, None]
        self.v_midpoint = self.bounds.x + self.bounds.width/2
        self.h_midpoint = self.bounds.y + self.bounds.height/2

    def clear(self):
        """
        Recursively remove object and nodes from the tree
        """
        self.members = 0
        self.objects.clear()
        for i in range(len(self.nodes)):
            if self.nodes[i] is not None:
                self.nodes[i].clear()
                self.nodes[i] = None

    def split(self):
        """
        Generates new nodes if the current node has more object than specified by max_objects
        """
        sub_width = self.bounds.width/2
        sub_height = self.bounds.height/2
        x = self.bounds.x
        y = self.bounds.y

        self.nodes[0] = Quadtree(self.level + 1, Rectangle(x + sub_width, y, sub_width, sub_height))
        self.nodes[1] = Quadtree(self.level + 1, Rectangle(x, y, sub_width, sub_height))
        self.nodes[2] = Quadtree(self.level + 1, Rectangle(x, y + sub_height, sub_width, sub_height))
        self.nodes[3] = Quadtree(self.level + 1, Rectangle(x + sub_width, y + sub_height, sub_width, sub_height))

    def get_index(self, obj: object):
        """
        Determine which node the current object belongs to if any
        :param obj: An object that as an attribute aabb or axis-aligned bounding box
        """
        bounds = obj.aabb

        if bounds.x > self.v_midpoint:
            if bounds.y > self.h_midpoint:
                return 3
            elif bounds.y + bounds.height < self.h_midpoint:
                return 0

        elif bounds.x + bounds.width < self.v_midpoint:
            if bounds.y > self.h_midpoint:
                return 2
            elif bounds.y + bounds.height < self.h_midpoint:
                return 3

        return -1

    def insert(self, obj: object):
        """
        Places an object into its correct node of the quadtree. This function will also
        split the tree if too many object are present at one node
        :param obj: an object that as the attribute aabb or axis-aligned bounding box
        """
        self.members += 1
        if self.nodes[0] is not None:
            index = self.get_index(obj)

            if index != -1:
                self.nodes[index].insert(obj)

                return

        self.objects.append(obj)

        if len(self.objects) > Quadtree.max_objects and self.level < Quadtree.max_levels:
            if self.nodes[0] is None:
                self.split()

            i = 0
            while i < len(self.objects):
                index = self.get_index(self.objects[i])
                if index != -1:
                    self.nodes[index].insert(self.objects.pop(i))
                else:
                    i += 1

    def set_leaf(self):
        """
        Marks all objects that are at a leaf as such
        """
        if self.nodes[0] is not None:
            self.nodes[0].set_leaf()
            self.nodes[1].set_leaf()
            self.nodes[2].set_leaf()
            self.nodes[3].set_leaf()
            return

        for obj in self.objects:
            obj.leaf = True

        return

    def retreive(self, return_obj: list, obj: object):
        """
        Return all objects that many intersect with an object of interest
        :param return_obj: An empty list that is filled by retreive
        :param obj: An object of interest that has attribute aabb or axis-aligned bounding box
        """
        index = self.get_index(obj)
        if index != -1 and self.nodes[0] is not None:
            self.nodes[index].retreive(return_obj, obj)

        return_obj.extend(self.objects)

        return return_obj


