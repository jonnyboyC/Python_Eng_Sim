# Credit to Steven Lambert for the tutorial

class Rectangle:

    def __init__(self, x: float, y: float, width: float, height: float, idx=-1):
        self.x = float(x)
        self.y = float(y)
        self.width = float(width)
        self.height = float(height)


class Quadtree:

    max_objects = 10
    max_levels = 7

    def __init__(self, p_level: int, p_bounds: Rectangle):
        self.level = int(p_level)
        self.objects = []
        self.bounds = p_bounds
        self.members = 0
        self.nodes = [None, None, None, None]
        self.v_midpoint = self.bounds.x + self.bounds.width/2
        self.h_midpoint = self.bounds.y + self.bounds.height/2

    def clear(self):
        self.members = 0
        self.objects.clear()
        for i in range(len(self.nodes)):
            if self.nodes[i] is not None:
                self.nodes[i].clear()
                self.nodes[i] = None

    def split(self):
        subwidth = self.bounds.width/2
        subheight = self.bounds.height/2
        x = self.bounds.x
        y = self.bounds.y

        self.nodes[0] = Quadtree(self.level + 1, Rectangle(x + subwidth, y, subwidth, subheight))
        self.nodes[1] = Quadtree(self.level + 1, Rectangle(x, y, subwidth, subheight))
        self.nodes[2] = Quadtree(self.level + 1, Rectangle(x, y + subheight, subwidth, subheight))
        self.nodes[3] = Quadtree(self.level + 1, Rectangle(x + subwidth, y + subheight, subwidth, subheight))

    def get_index(self, object: object):
        """
        TODO: Optimize

        """
        bounds = object.aabb

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

    def insert(self, object: object):
        self.members += 1
        if self.nodes[0] is not None:
            index = self.get_index(object)

            if index != -1:
                self.nodes[index].insert(object)

                return

        self.objects.append(object)

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
        if self.nodes[0] is not None:
            self.nodes[0].set_leaf()
            self.nodes[1].set_leaf()
            self.nodes[2].set_leaf()
            self.nodes[3].set_leaf()
            return

        for obj in self.objects:
            obj.leaf = True

        return

    def retreive(self, return_obj: list, object: object):
        index = self.get_index(object)
        if index != -1 and self.nodes[0] is not None:
            self.nodes[index].retreive(return_obj, object)

        return_obj.extend(self.objects)

        return return_obj

# quad = Quadtree(0, Rectangle(0, 0, 100, 100))


